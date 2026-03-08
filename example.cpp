// 1. Define the implementation macro exactly ONCE before including your library
#define STEGO_IMPLEMENTATION

// Include dr_mp3 implementation here if you are using the STB version of it as well
#define DR_MP3_IMPLEMENTATION 
#include "dr_mp3.h"

#include "stego_audio.hpp"

#include <iostream>
#include <fstream>
#include <string>

// Helper to load a real MP3 file into memory for the test
std::vector<uint8_t> read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) return {};
    std::vector<uint8_t> buf(f.tellg());
    f.seekg(0, std::ios::beg);
    f.read(reinterpret_cast<char*>(buf.data()), buf.size());
    return buf;
}

int main() {
    std::cout << "--- Steganography Library Example ---\n";

    // 1. Get our raw data
    std::vector<uint8_t> clean_mp3_buffer = read_file("yana.mp3");
    if (clean_mp3_buffer.empty()) {
        std::cerr << "Could not find yana.mp3 to run the test.\n";
        return 1;
    }

    std::string my_secret_string = "Hello from the pure-memory API! This leaves no trace on disk.";
    std::vector<uint8_t> secret_data(my_secret_string.begin(), my_secret_string.end());

    // 2. Configure the DSP parameters
    StegoConfig my_config;
    my_config.chip_length = 512;
    my_config.alpha = 0.05f;      // Slight bump to be sure
    my_config.low_bin = 150;      // Pushed lower into the safe midrange
    my_config.high_bin = 350;

    // 3. Embed (Using template parameters for RS Code Length and RS Data Length)
    std::cout << "Embedding data in memory...\n";
    std::vector<uint8_t> new_mp3_buffer = stego_embed_memory<255, 127>(clean_mp3_buffer, secret_data, my_config);

    if (new_mp3_buffer.empty()) {
        std::cerr << "Embedding failed (carrier likely too small).\n";
        return 1;
    }
    std::cout << "Successfully generated a " << new_mp3_buffer.size() << " byte stego MP3 in RAM.\n";

    // 4. Extract directly from the RAM buffer
    std::cout << "Extracting data from memory...\n";
    
    // We must use the exact same template parameters to decode properly!
    std::vector<uint8_t> extracted_data = stego_extract_memory<255, 127>(new_mp3_buffer, my_config);

    // 5. Verify
    std::string recovered_string(extracted_data.begin(), extracted_data.end());
    std::cout << "\n--- Recovered Secret ---\n";
    std::cout << recovered_string << "\n";

    return 0;
}