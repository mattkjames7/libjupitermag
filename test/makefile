# CFLAGS for CC
CFLAGS=-lm -std=c++17 -Wextra

# Compiler
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)

all:
	echo $(ORIGIN)
	$(CC) test.cc -o test	-Wl,-rpath='$$ORIGIN/../lib/libjupitermag' -L ../lib/libjupitermag -ljupitermag
	./test

clean:
	-rm -v test