# CFLAGS for CC
CFLAGS=-lm -std=c++17 -Wextra

# Compiler
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)

ifndef BUILDDIR 
	BUILDDIR=.
endif

all: obj lib

obj:
	$(CCo) spline.cc -o $(BUILDDIR)/spline.o
	$(CCo) libspline.cc -o $(BUILDDIR)/libspline.o

windows:
	$(CCWo) spline.cc -o $(BUILDDIR)/spline.o
	$(CCWo) libspline.cc -o $(BUILDDIR)/libspline.o
		
lib: 
	$(CC) -shared -o libspline.so libspline.o spline.o
	
clean:
	rm -v *.o *.so
