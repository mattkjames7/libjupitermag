CFLAGS=-std=c++17 -O3 -Wextra
IFLAGS=-I ../../include

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	LDFLAGS=-Wl,-rpath='$$ORIGIN/../../lib' -L ../../lib -ljupitermag
else
	LDFLAGS=-L ../../lib -ljupitermag
endif

CC=g++ $(CFLAGS)

all: generate_baseline
	./generate_baseline

generate_baseline: generate_baseline.cc
	$(CC) $(IFLAGS) generate_baseline.cc -o generate_baseline $(LDFLAGS)

clean:
	-rm -f generate_baseline

clean-data:
	-rm -rf ../data
