# CFLAGS for CC
CFLAGS=-lm -std=c++17 -Wextra -ggdb3 -O3  -no-pie  #-fsanitize=address

# Compiler
CCo=g++ -fPIC -c $(CFLAGS)
CC=g++ -fPIC $(CFLAGS)
CCWo=x86_64-w64-mingw32-g++-win32 -c $(CFLAGS)
CCW=x86_64-w64-mingw32-g++-win32 $(CFLAGS)
GTEST_DIR=external/googletest
GTEST_REPO=https://github.com/google/googletest.git
GTEST_TAG=v1.14.0

ifeq ($(OS),Windows_NT)
	LDFLAGS=-L ../lib -ljupitermag
else 
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Darwin)
		LDFLAGS=-L ../lib -ljupitermag
	else
		LDFLAGS=-Wl,-rpath='$$ORIGIN/../lib' -L ../lib -ljupitermag  -static-libasan
	endif
endif

all:
	$(CC)  -fsanitize=address  -I ../include cpptest.cc -o cpptest	$(LDFLAGS)
	$(CC) timing.cc -o timing	$(LDFLAGS)
	gcc -I ../include ctest.c -o ctest -lm -g $(LDFLAGS) 

	./ctest
	./cpptest
#   ./timing

clean:
	-rm -v ctest
	-rm -v cpptest
	-rm -v timing
	-rm -v regression_gtest

setup-gtest:
	@mkdir -p external
	@if [ ! -d $(GTEST_DIR)/.git ]; then \
		echo "Cloning googletest $(GTEST_TAG) into $(GTEST_DIR)"; \
		git clone --depth 1 --branch $(GTEST_TAG) $(GTEST_REPO) $(GTEST_DIR); \
	else \
		echo "googletest already exists in $(GTEST_DIR)"; \
	fi

gtest: setup-gtest
	$(CC) -I ../include \
		-I $(GTEST_DIR)/googletest/include \
		-I $(GTEST_DIR)/googletest \
		regression_gtest.cc \
		$(GTEST_DIR)/googletest/src/gtest-all.cc \
		$(GTEST_DIR)/googletest/src/gtest_main.cc \
		-o regression_gtest $(LDFLAGS) -pthread
	./regression_gtest