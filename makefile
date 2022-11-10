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

ifeq ($(PREFIX),)
#install path
	PREFIX=/usr/local
endif

.PHONY: all lib clean header test

all: obj internal con2020 spline lib header

internal:
	cd lib/libinternalfield; make obj

con2020:
	cd lib/libcon2020; make obj

spline:
	cd lib/libspline; make obj

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
	cd lib/libinternalfield; make winobj
	cd lib/libcon2020; make winobj
	cd lib/libspline; make winobj
	cd src; make winobj
	$(MD) lib/libjupitermag
	cd src; make winlib

test:
	cd test; make all

install:
	cp -v include/jupitermag.h $(PREFIX)/include
	cp -v lib/libjupitermag/libjupitermag.so $(PREFIX)/lib
	chmod 0775 $(PREFIX)/lib/libjupitermag.so
	ldconfig

uninstall:
	rm -v $(PREFIX)/include/jupitermag.h
	rm -v $(PREFIX)/lib/libjupitermag.so
	ldconfig

clean:
	cd lib/libinternalfield; make clean
	cd lib/libcon2020; make clean
	cd lib/libspline; make clean
	-rm -v build/*.o
	-rmdir -v build
	cd test; make clean