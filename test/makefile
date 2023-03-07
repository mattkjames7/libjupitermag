# CFLAGS for CC
CFLAGS=-lm -std=c++17 -Wextra -ggdb3 -O3  -no-pie  #-fsanitize=address

# Compiler
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)

ifeq ($(OS),Windows_NT)
	LDFLAGS=-L ../lib -ljupitermag
else 
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		LDFLAGS=-L ../lib -ljupitermag
	else
		LDFLAGS=-Wl,-rpath='$$ORIGIN/../lib' -L ../lib -ljupitermag  -static-libasan
	endif
endif

all:
	$(CC)  -fsanitize=address  -I ../include cpptest.cc -o cpptest	$(LDFLAGS)
	$(CC) timing.cc -o timing	$(LDFLAGS)
	gcc -I ../include ctest.c -o ctest -lm -g $(LDFLAGS) 

	./ctest
	./cpptest
	./timing

clean:
	-rm -v ctest
	-rm -v cpptest
	-rm -v timing