# CFLAGS for CC
ifdef DEBUG
	CFLAGS=-lm -std=c++17 -Wextra -ggdb3  -no-pie -static-libasan
else 
	CFLAGS=-lm -std=c++17 -Wextra 
endif

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
	$(CCo) footprint.cc -o $(BUILDDIR)/footprint.o
	$(CCo) coordconv.cc -o $(BUILDDIR)/coordconv.o 
	
lib:
	
	$(CC) -shared libjupitermag.cc $(BUILDDIR)/*.o \
			-o ../lib/$(LIBFILE) 

winobj:
#	cd con2020;	make windows
#	cd spline; make windows
	$(CCWo) model.cc -o $(BUILDDIR)/model.o
	$(CCWo) trace.cc -o $(BUILDDIR)/trace.o
	$(CCWo) interptraceclosestpos.cc -o $(BUILDDIR)/interptraceclosestpos.o	
	$(CCWo) footprint.cc -o $(BUILDDIR)/footprint.o
	$(CCWo) coordconv.cc -o $(BUILDDIR)/coordconv.o 

winlib:
	$(CCW) -shared libjupitermag.cc $(BUILDDIR)/*.o -o ../lib/libjupitermag.dll



