# CFLAGS for CC
ifdef DEBUG
	CFLAGS=-lm -std=c++17 -Wextra -ggdb3 -no-pie 
else 
	CFLAGS=-lm -std=c++17 -Wextra
endif 


# Compilers - the mingw ones allow us to compile for Windows using Linux
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)


export BUILDDIR=$(shell pwd)/build
ifeq ($(OS),Windows_NT)
#windows stuff here
	MD=mkdir
	LIBFILE=libjupitermag.dll
	ARCFILE=jupitermag.lib
else
#linux and mac here
	OS=$(shell uname -s)
	ifeq ($(OS),Linux)
		LIBFILE=libjupitermag.so
	else
		LIBFILE=libjupitermag.dylib
	endif
	ARCFILE=libjupitermag.a
	MD=mkdir -p
endif

#install paths
ifeq ($(PREFIX),)
	PREFIX=/usr/local
endif
ifeq ($(PREFIX_INC),)
	PREFIX_INC=$(PREFIX)/include
endif
ifeq ($(PREFIX_LIB),)
	PREFIX_LIB=$(PREFIX)/lib
endif


.PHONY: all lib clean header test internal con2020 spline obj install uninstall  windows

all: internal con2020 spline
	$(MD) $(BUILDDIR)
	$(MD) lib
	cd src; make all

internal:
	cd lib/libinternalfield; make header obj

con2020:
	cd lib/libcon2020; make header obj

spline:
	cd lib/libspline; make header obj

obj:
	$(MD) $(BUILDDIR)
	cd src; make obj

lib:
	cd src; make lib

header:
	cd src; make header
# ifneq (,$(shell which python3))
# 	python3 generateheader.py
# else
# 	@echo "python3 command doesn't appear to exist - skipping header regeneration..."
# endif

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
	cp -v include/jupitermag.h $(PREFIX_INC)
	cp -v lib/$(LIBFILE) $(PREFIX_LIB)
	cp -v lib/$(ARCFILE) $(PREFIX_LIB)
	chmod 0775 $(PREFIX_LIB)/$(LIBFILE)
ifeq ($(OS),Linux)
	-ldconfig
endif

uninstall:
	rm -v $(PREFIX_INC)/jupitermag.h
	rm -v $(PREFIX_LIB)/$(LIBFILE)
ifeq ($(OS),Linux)
	-ldconfig
endif

clean:
	cd lib/libinternalfield; make clean
	cd lib/libcon2020; make clean
	cd lib/libspline; make clean
	-rm -vfr build
	cd test; make clean
	rm -v lib/$(LIBFILE)
