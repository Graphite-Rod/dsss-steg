// Public API
#ifndef STEGO_AUDIO_H
#define STEGO_AUDIO_H

#include <vector>
#include <cstdint>
#include <cstddef>

struct StegoConfig {
    int chip_length = 512;
    float alpha = 0.04f;
    int low_bin = 200;
    int high_bin = 400;
	int seed = 42;
};

// Embeds bits into an MP3 memory buffer. Returns a new MP3 memory buffer.
template<std::size_t RS_CODE_LEN, std::size_t RS_DATA_LEN>
std::vector<uint8_t> stego_embed_memory(
    const std::vector<uint8_t>& container_data, 
    const std::vector<uint8_t>& secret_data, 
    const StegoConfig& config);

// Extracts bits from an MP3 memory buffer.
template<std::size_t RS_CODE_LEN, std::size_t RS_DATA_LEN>
std::vector<uint8_t> stego_extract_memory(
    const std::vector<uint8_t>& container_data, 
    const StegoConfig& config);

#endif // STEGO_AUDIO_H

// General implementation
#ifdef STEGO_IMPLEMENTATION

#include "dr_mp3.h"
#include <lame/lame.h>

#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <complex>
#include <utility>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "schifra/schifra_galois_field.hpp"
#include "schifra/schifra_galois_field_polynomial.hpp"
#include "schifra/schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra/schifra_reed_solomon_encoder.hpp"
#include "schifra/schifra_reed_solomon_decoder.hpp"
#include "schifra/schifra_reed_solomon_block.hpp"

namespace stego_internal {
    const std::vector<int> BARKER_13 = {1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1};

    inline std::vector<uint8_t> packBits(const std::vector<int>& bitstream) {
        std::vector<uint8_t> bytestream;
        uint8_t current_byte = 0;
        int bit_count = 0;
        for (int bit : bitstream) {
            current_byte = (current_byte << 1) | (bit & 1);
            if (++bit_count == 8) {
                bytestream.push_back(current_byte);
                current_byte = 0; bit_count = 0;
            }
        }
        if (bit_count > 0) bytestream.push_back(current_byte << (8 - bit_count));
        return bytestream;
    }

    inline std::vector<int> unpackBits(const std::vector<uint8_t>& bytestream, int original_bits) {
        std::vector<int> bitstream;
        bitstream.reserve(original_bits);
        int extracted = 0;
        for (uint8_t byte : bytestream) {
            for (int i = 7; i >= 0; --i) {
                if (extracted++ >= original_bits) break;
                bitstream.push_back((byte >> i) & 1);
            }
        }
        return bitstream;
    }

	typedef std::complex<float> Complex;

	// ---------------------------------------------------------
	// Helper: In-place Radix-2 FFT/IFFT
	// ---------------------------------------------------------
	inline void fft(std::vector<Complex>& a, bool invert) {
		int n = a.size();
		// Bit-reversal permutation
		for (int i = 1, j = 0; i < n; i++) {
			int bit = n >> 1;
			for (; j & bit; bit >>= 1) j ^= bit;
			j ^= bit;
			if (i < j) std::swap(a[i], a[j]);
		}
		// Cooley-Tukey decimation-in-time
		for (int len = 2; len <= n; len <<= 1) {
			float angle = 2.0f * M_PI / len * (invert ? 1.0f : -1.0f);
			Complex wlen(std::cos(angle), std::sin(angle));
			for (int i = 0; i < n; i += len) {
				Complex w(1.0f, 0.0f);
				for (int j = 0; j < len / 2; j++) {
					Complex u = a[i + j];
					Complex v = a[i + j + len / 2] * w;
					a[i + j] = u + v;
					a[i + j + len / 2] = u - v;
					w *= wlen;
				}
			}
		}
		if (invert) {
			for (Complex& x : a) x /= static_cast<float>(n);
		}
	}

	inline bool isPowerOfTwo(size_t n) {
		return n > 0 && (n & (n - 1)) == 0;
	}

	// ---------------------------------------------------------
	// Original Naive Fallbacks (for non-power-of-2 sizes)
	// ---------------------------------------------------------
	inline std::vector<float> applyDCT_Naive(const std::vector<float>& input) {
		size_t N = input.size();
		std::vector<float> output(N, 0.0f);
		for (size_t k = 0; k < N; ++k) {
			float sum = 0.0f;
			for (size_t n = 0; n < N; ++n) sum += input[n] * std::cos(M_PI / N * (n + 0.5f) * k);
			output[k] = sum;
		}
		return output;
	}

	inline std::vector<float> applyIDCT_Naive(const std::vector<float>& input) {
		size_t N = input.size();
		std::vector<float> output(N, 0.0f);
		for (size_t n = 0; n < N; ++n) {
			float sum = input[0] / 2.0f; 
			for (size_t k = 1; k < N; ++k) sum += input[k] * std::cos(M_PI / N * k * (n + 0.5f));
			output[n] = sum * (2.0f / N);
		}
		return output;
	}

	// ---------------------------------------------------------
	// Fast O(N log N) Implementations
	// ---------------------------------------------------------
	inline std::vector<float> applyDCT(const std::vector<float>& input) {
		size_t N = input.size();
		if (N == 0) return {};
		if (!isPowerOfTwo(N)) return applyDCT_Naive(input);

		// 1. Permute input (Makhoul's arrangement)
		std::vector<Complex> v(N);
		for (size_t n = 0; n < N / 2; ++n) {
			v[n] = input[2 * n];
			v[N - 1 - n] = input[2 * n + 1];
		}

		// 2. Compute N-point FFT
		fft(v, false);

		// 3. Phase rotation and extract real part
		std::vector<float> output(N);
		for (size_t k = 0; k < N; ++k) {
			float phase = -M_PI * k / (2.0f * N);
			Complex rot(std::cos(phase), std::sin(phase));
			output[k] = (v[k] * rot).real();
		}
		return output;
	}

	inline std::vector<float> applyIDCT(const std::vector<float>& input) {
		size_t N = input.size();
		if (N == 0) return {};
		if (!isPowerOfTwo(N)) return applyIDCT_Naive(input);

		// 1. Construct 2N-point conjugate symmetric array
		std::vector<Complex> W(2 * N, Complex(0.0f, 0.0f));
		W[0] = input[0];
		for (size_t k = 1; k < N; ++k) {
			float phase = M_PI * k / (2.0f * N);
			W[k] = input[k] * Complex(std::cos(phase), std::sin(phase));
			W[2 * N - k] = std::conj(W[k]);
		}

		// 2. Compute 2N-point IFFT
		fft(W, true); 

		// 3. Extract first N points and scale
		std::vector<float> output(N);
		for (size_t n = 0; n < N; ++n) {
			output[n] = 2.0f * W[n].real();
		}
		return output;
	}

    inline std::vector<float> genPN(size_t len, int min_b, int max_b, int pn_seed) {
        std::vector<float> pn(len);
        std::mt19937 gen(pn_seed); 
        std::uniform_int_distribution<> dis(0, 1);
        for (size_t i = 0; i < len; ++i) pn[i] = dis(gen) ? 1.0f : -1.0f;

        std::vector<float> pn_freq = applyDCT(pn);
        for (size_t k = 0; k < len; ++k) { if (k < min_b || k > max_b) pn_freq[k] = 0.0f; }
        
        std::vector<float> pn_filt = applyIDCT(pn_freq);
        float m = 0.0f;
        for (float v : pn_filt) m = std::max(m, std::abs(v));
        if (m > 0.0f) for (size_t i = 0; i < len; ++i) pn_filt[i] /= m;
        return pn_filt;
    }
}

template<std::size_t RS_CODE_LEN, std::size_t RS_DATA_LEN>
std::vector<uint8_t> stego_embed_memory(const std::vector<uint8_t>& container_data, const std::vector<uint8_t>& secret_data, const StegoConfig& config) {
    constexpr std::size_t RS_FEC_LEN = RS_CODE_LEN - RS_DATA_LEN;
    
    // Padding & RS Encoding
    std::vector<uint8_t> data_padded = secret_data;
    while (data_padded.size() % RS_DATA_LEN != 0) data_padded.push_back(0x00);
    
    schifra::galois::field field(8, schifra::galois::primitive_polynomial_size06, schifra::galois::primitive_polynomial06);
    schifra::galois::field_polynomial gen_poly(field);
    if (!schifra::make_sequential_root_generator_polynomial(field, 0, RS_FEC_LEN, gen_poly)) return {};
    schifra::reed_solomon::encoder<RS_CODE_LEN, RS_FEC_LEN> encoder(field, gen_poly);

    std::vector<uint8_t> encoded_bytes;
    encoded_bytes.reserve((data_padded.size() / RS_DATA_LEN) * RS_CODE_LEN);
    for (size_t b = 0; b < data_padded.size() / RS_DATA_LEN; ++b) {
        schifra::reed_solomon::block<RS_CODE_LEN, RS_FEC_LEN> block;
        for (size_t i = 0; i < RS_DATA_LEN; ++i) block[i] = data_padded[b * RS_DATA_LEN + i];
        encoder.encode(block);
        for (size_t i = 0; i < RS_CODE_LEN; ++i) encoded_bytes.push_back(block[i]);
    }

    std::vector<int> payload_bits = stego_internal::unpackBits(encoded_bytes, encoded_bytes.size() * 8);
    std::vector<int> tx_stream = stego_internal::BARKER_13;
    uint32_t p_len = payload_bits.size();
    for (int i = 31; i >= 0; --i) tx_stream.push_back((p_len >> i) & 1);
    tx_stream.insert(tx_stream.end(), payload_bits.begin(), payload_bits.end());

    // Decode Audio
    drmp3 mp3;
    if (!drmp3_init_memory(&mp3, container_data.data(), container_data.size(), NULL)) return {};
    std::vector<float> pcm;
    std::vector<float> chunk(4096 * mp3.channels);
    drmp3_uint64 frames;
    while ((frames = drmp3_read_pcm_frames_f32(&mp3, 4096, chunk.data())) > 0) {
        pcm.insert(pcm.end(), chunk.begin(), chunk.begin() + frames * mp3.channels);
    }
    
    size_t total_frames = pcm.size() / mp3.channels;
	std::cout << tx_stream.size() * config.chip_length << "/" << total_frames << std::endl;
    if (tx_stream.size() * config.chip_length > total_frames) {
        std::cerr << "ERROR: Carrier audio too small.\n";
        drmp3_uninit(&mp3);
        return {};
    }

    std::vector<float> left(total_frames), right(total_frames);
    for (size_t i = 0; i < total_frames; ++i) {
        left[i] = pcm[i * mp3.channels];
        right[i] = (mp3.channels > 1) ? pcm[i * mp3.channels + 1] : pcm[i * mp3.channels];
    }
    drmp3_uninit(&mp3);

    // Embedding
    std::vector<float> pn_seq = stego_internal::genPN(config.chip_length, config.low_bin, config.high_bin, config.seed);
    for (size_t b = 0; b < tx_stream.size(); ++b) {
        float bit_val = (tx_stream[b] == 1) ? 1.0f : -1.0f;
        size_t off = b * config.chip_length;
        
        std::vector<float> host_chip(config.chip_length);
        for (size_t i = 0; i < config.chip_length; ++i) {
            host_chip[i] = left[off + i];
        }
        std::vector<float> host_freq = stego_internal::applyDCT(host_chip);
        
        // Calculate energy ONLY in the target band
        float band_energy = 0.0f;
        int band_count = config.high_bin - config.low_bin + 1;
        for (int k = config.low_bin; k <= config.high_bin; ++k) {
            band_energy += std::abs(host_freq[k]);
        }
        
        // 0.0005f since high-freq bands naturally have less energy
        float energy = std::max(band_energy / band_count, 0.0005f); 
        
        for (size_t i = 0; i < config.chip_length; ++i) {
            float alpha = config.alpha * energy;
            left[off + i] += alpha * bit_val * pn_seq[i];
            right[off + i] += alpha * bit_val * pn_seq[i];
        }
    }

    // MP3 Encoding
    lame_t lame = lame_init();
    lame_set_in_samplerate(lame, mp3.sampleRate);
    lame_set_num_channels(lame, 2);
    lame_set_VBR(lame, vbr_default); 
    lame_set_VBR_q(lame, 0); // Force High Quality
    lame_init_params(lame);

    std::vector<uint8_t> out_mp3(left.size() * 1.25 + 7200); // whyyyy?
    int w1 = lame_encode_buffer_ieee_float(lame, left.data(), right.data(), left.size(), out_mp3.data(), out_mp3.size());
    int w2 = lame_encode_flush(lame, out_mp3.data() + w1, out_mp3.size() - w1);
    out_mp3.resize(w1 + w2);
    lame_close(lame);

    return out_mp3;
}

template<std::size_t RS_CODE_LEN, std::size_t RS_DATA_LEN>
std::vector<uint8_t> stego_extract_memory(const std::vector<uint8_t>& container_data, const StegoConfig& config) {
    constexpr std::size_t RS_FEC_LEN = RS_CODE_LEN - RS_DATA_LEN;

    // Decode Audio
    drmp3 mp3;
    if (!drmp3_init_memory(&mp3, container_data.data(), container_data.size(), NULL)) return {};
    std::vector<float> target(drmp3_get_pcm_frame_count(&mp3));
    std::vector<float> chunk(4096 * mp3.channels);
    drmp3_uint64 frames, idx = 0;
    while ((frames = drmp3_read_pcm_frames_f32(&mp3, 4096, chunk.data())) > 0) {
        for(drmp3_uint64 i = 0; i < frames; ++i) target[idx++] = chunk[i * mp3.channels];
    }
    drmp3_uninit(&mp3);

    // Sync & Extraction
    std::vector<float> pn_seq = stego_internal::genPN(config.chip_length, config.low_bin, config.high_bin, config.seed);
    std::vector<float> barker_wf(stego_internal::BARKER_13.size() * config.chip_length);
    for (size_t b = 0; b < stego_internal::BARKER_13.size(); ++b) {
        float bit_val = (stego_internal::BARKER_13[b] == 1) ? 1.0f : -1.0f;
        for (size_t i = 0; i < config.chip_length; ++i) barker_wf[b * config.chip_length + i] = bit_val * pn_seq[i];
    }

    int best_offset = 0; float max_corr = 0.0f;
    for (int off = 0; off < 8000 && off + barker_wf.size() < target.size(); ++off) {
        float c = 0.0f;
        for (size_t i = 0; i < barker_wf.size(); ++i) c += target[off + i] * barker_wf[i];
        if (c > max_corr) { max_corr = c; best_offset = off; }
    }

    size_t cur = best_offset + barker_wf.size();
    uint32_t p_len = 0;
    for (int b = 0; b < 32; ++b) {
        float c = 0.0f;
        for (size_t i = 0; i < config.chip_length; ++i) c += target[cur + i] * pn_seq[i];
        p_len = (p_len << 1) | ((c > 0) ? 1 : 0);
        cur += config.chip_length;
    }

    if (p_len == 0 || p_len > (target.size() / config.chip_length)) return {};

    std::vector<int> extracted; extracted.reserve(p_len);
    for (size_t b = 0; b < p_len; ++b) {
        float c = 0.0f;
        for (size_t i = 0; i < config.chip_length; ++i) c += target[cur + i] * pn_seq[i];
        extracted.push_back((c > 0) ? 1 : 0);
        cur += config.chip_length;
    }

    // RS Decode
    std::vector<uint8_t> ex_bytes = stego_internal::packBits(extracted);
    schifra::galois::field field(8, schifra::galois::primitive_polynomial_size06, schifra::galois::primitive_polynomial06);
    schifra::reed_solomon::decoder<RS_CODE_LEN, RS_FEC_LEN> decoder(field, 0);
    
    std::vector<uint8_t> result;
    size_t blocks = ex_bytes.size() / RS_CODE_LEN;
    for (size_t b = 0; b < blocks; ++b) {
        schifra::reed_solomon::block<RS_CODE_LEN, RS_FEC_LEN> block;
        for (size_t i = 0; i < RS_CODE_LEN; ++i) block[i] = ex_bytes[b * RS_CODE_LEN + i];
        decoder.decode(block);
        for (size_t i = 0; i < RS_DATA_LEN; ++i) result.push_back(block[i]);
    }

    while(!result.empty() && result.back() == 0x00) result.pop_back();
    return result;
}

#endif // STEGO_IMPLEMENTATION

// Command-Line-Interface
#ifdef STEGO_CLI
#include <fstream>
#include <string>

std::vector<uint8_t> readFile(const std::string& path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) { std::cerr << "Failed to open: " << path << "\n"; exit(1); }
    std::vector<uint8_t> buf(f.tellg());
    f.seekg(0, std::ios::beg);
    f.read(reinterpret_cast<char*>(buf.data()), buf.size());
    return buf;
}

void writeFile(const std::string& path, const std::vector<uint8_t>& data) {
    std::ofstream f(path, std::ios::binary);
    f.write(reinterpret_cast<const char*>(data.data()), data.size());
}

int main(int argc, char* argv[]) {
    if (argc < 7) {
        std::cerr << "Usage: \n"
                  << "  Embed:   dsss_cli embed --container <in.mp3> --seed <int> --chip <int> --lo_bin <int> --hi_bin <int> --secret <data.txt> --out <out.mp3>\n"
                  << "  Extract: dsss_cli extract --container <stego.mp3> --seed <int> --chip <int> --lo_bin <int> --hi_bin <int> --out <data.txt>\n";
        return 1;
    }

    std::string mode = argv[1];
    std::string container, secret, out;
	int seed, chip, hi_bin, lo_bin;

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--container" && i+1 < argc) container = argv[++i];
        else if (arg == "--secret" && i+1 < argc) secret = argv[++i];
        else if (arg == "--out" && i+1 < argc) out = argv[++i];
        else if (arg == "--seed" && i+1 < argc) seed = atoi(argv[++i]);
        else if (arg == "--chip" && i+1 < argc) chip = atoi(argv[++i]);
        else if (arg == "--lo_bin" && i+1 < argc) lo_bin = atoi(argv[++i]);
        else if (arg == "--hi_bin" && i+1 < argc) hi_bin = atoi(argv[++i]);
    }

    StegoConfig config;
	config.chip_length = chip;
	config.low_bin = lo_bin;
	config.high_bin = hi_bin;
	config.seed = seed;
	
	if (lo_bin>chip | hi_bin>chip) {
		std::cerr << "ERROR: Chosen frequency band overcomes the range defined by the chip size.\n";
	}
    if (mode == "embed") {
        std::cout << "Loading container and secret...\n";
        auto c_data = readFile(container);
        auto s_data = readFile(secret);
        std::cout << "Embedding (this may take a moment)...\n";
        auto result = stego_embed_memory<255, 191>(c_data, s_data, config);
        if (result.empty()) return 1;
        writeFile(out, result);
        std::cout << "Success! Saved to " << out << "\n";
    } 
    else if (mode == "extract") {
        std::cout << "Loading container...\n";
        auto c_data = readFile(container);
        std::cout << "Extracting...\n";
        auto result = stego_extract_memory<255, 191>(c_data, config);
        if (result.empty()) { std::cerr << "Extraction failed.\n"; return 1; }
        writeFile(out, result);
        std::cout << "Success! Extracted " << result.size() << " bytes to " << out << "\n";
    }

    return 0;
}
#endif // STEGO_CLI