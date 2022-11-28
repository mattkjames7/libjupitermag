# CFLAGS for CC
CFLAGS=-lm -std=c++17 -Wextra  -O3

# Compiler
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)

ifndef BUILDDIR 
	BUILDDIR=.
endif

OS=$(shell uname -s)
ifeq ($(OS),Linux)
	LIBFILE=libjupitermag.so
else
	LIBFILE=libjupitermag.dylib
endif

all: obj lib


.PHONY: all lib clean

obj:
	$(CCo) model.cc -o $(BUILDDIR)/model.o
	$(CCo) trace.cc -o $(BUILDDIR)/trace.o
	$(CCo) interptraceclosestpos.cc -o $(BUILDDIR)/interptraceclosestpos.o

	
lib:
	
	$(CC) -shared libjupitermag.cc $(BUILDDIR)/*.o \
			-o ../lib/libjupitermag/$(LIBFILE)

test:
	$(CC) test.cc *.o \
			con2020/*.o \
			spline/*.o \
			-o test	-Wl,-rpath='$$ORIGIN/libinternalfield' -L ./libinternalfield -linternalfield


winobj:
#	cd con2020;	make windows
#	cd spline; make windows
	$(CCWo) model.cc -o $(BUILDDIR)/model.o
	$(CCWo) trace.cc -o $(BUILDDIR)/trace.o
	$(CCWo) interptraceclosestpos.cc -o $(BUILDDIR)/interptraceclosestpos.o	

winlib:
	$(CCW) -shared libjupitermag.cc $(BUILDDIR)/*.o -o ../lib/libjupitermag/libjupitermag.dll



