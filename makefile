# CFLAGS for CC
CFLAGS=-lm -std=c++17 -Wextra

# Compilers
defCCo=g++ -fPIC -c $(CFLAGS)
defCC=g++ -fPIC $(CFLAGS)
mingwCCo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
mingwCC=x86_64-w64-mingw32-g++-win32 $(CFLAGS)

#chose compiler based on OS
ifeq ($(OS),Windows_NT)
#windows stuff here
	BUILDDIR=build/windows
else
#linux and mac here
	BUILDDIR=build/unix
endif



all: obj lib


obj:
	echo $(BUILDDIR)
	mkdir -pv $(BUILDDIR)
	CC="$(defCC)"
	CCo="$(defCCo)"
	cd src; make all CC=$(CC) CCo="$(CCo)" BUILDDIR=$(BUILDDIR)

lib:
	echo "lib"


windows:
	BUILDDIR=build/windows
	CC=$(mingwCC)
	CCo=$(mingwCCo)
