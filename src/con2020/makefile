# CFLAGS for CC
CFLAGS=-lm -fPIC -std=c++17 -Wextra

# Compiler
CCo=g++ -c $(CFLAGS)
CC=g++ $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)

ifndef BUILDDIR 
	BUILDDIR=.
endif


all: obj lib 

obj:
	$(CCo) con2020.cc -o $(BUILDDIR)/con2020.o 
	$(CCo) bessel.cc -o $(BUILDDIR)/bessel.o 
	$(CCo) polyeval.cc -o $(BUILDDIR)/polyeval.o 
	$(CCo) trap.cc -o $(BUILDDIR)/trap.o 
	$(CCo) smoothd.cc -o $(BUILDDIR)/smoothd.o 
	$(CCo) libcon2020.cc -o $(BUILDDIR)/libcon2020.o

lib:
	$(CC) *.o -shared -o $(BUILDDIR)/libcon2020.so
	
windows:
	$(CCWo) con2020.cc -o $(BUILDDIR)/con2020.o 
	$(CCWo) bessel.cc -o $(BUILDDIR)/bessel.o 
	$(CCWo) polyeval.cc -o $(BUILDDIR)/polyeval.o 
	$(CCWo) trap.cc -o $(BUILDDIR)/trap.o 
	$(CCWo) smoothd.cc -o $(BUILDDIR)/smoothd.o 
	$(CCWo) libcon2020.cc -o $(BUILDDIR)/libcon2020.o

test:
	$(CC) test.cc *.o -o test

clean:
	rm -v *.o 
	rm -v *.so
	#rm test
