# CFLAGS for CC
CFLAGS=-lm -std=c++17 -Wextra -g -O3

# Compiler
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)

ifeq ($(OS),Windows_NT)
	LDFLAGS=-L ../lib/libjupitermag -ljupitermag
else 
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		LDFLAGS=-L ../lib/libjupitermag -ljupitermag
	else
		LDFLAGS=-Wl,-rpath='$$ORIGIN/../lib/libjupitermag' -L ../lib/libjupitermag -ljupitermag
	endif
endif

all:
	$(CC) test.cc -o test	$(LDFLAGS)
	$(CC) timing.cc -o timing	$(LDFLAGS)

	./test
	./timing

clean:
	-rm -v test
	-rm -v timing