@echo off
@set OUT_EXE=./example.exe
@set INCLUDES=-Iinclude -I.
@set LDFLAGS=-L. -lmp3lame
@set SOURCES=example.cpp

g++ -g -O3 %INCLUDES% %SOURCES% -o %OUT_EXE% %LDFLAGS%