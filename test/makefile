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
GTEST_SHARDS?=2
REGRESSION_COMMON_SRC=unit/regression_common.cc
REGRESSION_TEST_SRCS=$(wildcard unit/*_gtest.cc)
REGRESSION_TEST_BINS=$(patsubst unit/%.cc,%,$(REGRESSION_TEST_SRCS))

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
	-rm -v $(REGRESSION_TEST_BINS)

setup-gtest:
	@mkdir -p external
	@if [ ! -d $(GTEST_DIR)/.git ]; then \
		echo "Cloning googletest $(GTEST_TAG) into $(GTEST_DIR)"; \
		git clone --depth 1 --branch $(GTEST_TAG) $(GTEST_REPO) $(GTEST_DIR); \
	else \
		echo "googletest already exists in $(GTEST_DIR)"; \
	fi


regression_gtests: setup-gtest $(REGRESSION_TEST_BINS)

$(REGRESSION_TEST_BINS): %: unit/%.cc $(REGRESSION_COMMON_SRC)
	$(CC) -I ../include \
		-I $(GTEST_DIR)/googletest/include \
		-I $(GTEST_DIR)/googletest \
		$^ \
		$(GTEST_DIR)/googletest/src/gtest-all.cc \
		$(GTEST_DIR)/googletest/src/gtest_main.cc \
		-o $@ $(LDFLAGS) -pthread

regression_gtest: regression_gtests

gtest: regression_gtests
	@for t in $(REGRESSION_TEST_BINS); do \
		echo "Running $$t"; \
		./$$t; \
	done

gtest-parallel: regression_gtests
	@echo "Running gtest binaries in parallel"
	@for t in $(REGRESSION_TEST_BINS); do \
		./$$t & \
	done; \
	wait