@echo off
@set OUT_EXE=./dsss_cli.exe
@set INCLUDES=-Iinclude -I.
@set LDFLAGS=-L. -lmp3lame
@set SOURCES=-x c++ stego_audio.hpp
@set DEFINES=-DSTEGO_IMPLEMENTATION -DSTEGO_CLI -DDR_MP3_IMPLEMENTATION

g++ -g -O3 %INCLUDES% %DEFINES% %SOURCES% -o %OUT_EXE% %LDFLAGS%