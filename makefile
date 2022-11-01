# CFLAGS for CC
CFLAGS=-lm -std=c++17 -Wextra

# Compilers - the mingw ones allow us to compile for Windows using Linux
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)

#chose compiler based on OS
ifeq ($(OS),Windows_NT)
#windows stuff here
	export BUILDDIR=$(shell pwd)/build/windows
	MD=mkdir
	RBUILD=rm $(BUILDDIR)/*; rmdir $(BUILDDIR); rmdir build
	RLIB=rm lib/libinternalfield/libinternalfield.dll; rmdir lib/libinternalfield
else
#linux and mac here
	export BUILDDIR=$(shell pwd)/build/unix
	MD=mkdir -p
	RBUILD=rm -vfr $(BUILDIR)
	RLIB=rm lib/libinternalfield/libinternalfield.so; rmdir lib/libinternalfield
endif


.PHONY: all lib clean

all: internal obj lib

internal:
	cd lib/libinternalfield; make all

obj:
	$(MD) $(BUILDDIR)
	cd src; make obj

lib:
	$(MD) lib/libjupitermag
	cd src; make lib


windows:
	BUILDDIR=build/windows

clean:
	cd lib/libinternalfield; make clean
	$(RBUILD)
	$(RLIB)
	