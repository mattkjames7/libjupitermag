# CFLAGS for CC
CFLAGS=-lm -std=c++17 -Wextra

# Compilers - the mingw ones allow us to compile for Windows using Linux
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)


export BUILDDIR=$(shell pwd)/build
ifeq ($(OS),Windows_NT)
#windows stuff here
	MD=mkdir
else
#linux and mac here
	MD=mkdir -p
endif


.PHONY: all lib clean header

all: internal obj lib header

internal:
	cd lib/libinternalfield; make all

obj:
	$(MD) $(BUILDDIR)
	cd src; make obj

lib:
	$(MD) lib/libjupitermag
	cd src; make lib

header:
ifneq (,$(shell which python3))
	python3 generateheader.py
else
	@echo "python3 command doesn't appear to exist - skipping header regeneration..."
endif

windows:
#this should build the library for Windows using mingw
	$(MD) $(BUILDDIR)
	cd lib/libinternalfield; make windows
	cd src; make winobj
	$(MD) lib/libjupitermag
	cd src; make winlib


clean:
	cd lib/libinternalfield; make clean
	-rm -v lib/libjupitermag/libjupitermag.so
	-rm -v lib/libjupitermag/libjupitermag.dll
	-rm -vfr lib/libjupitermag
	-rm -v build/*.o
	-rmdir -v build