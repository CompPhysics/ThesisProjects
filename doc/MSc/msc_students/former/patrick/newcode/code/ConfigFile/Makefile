# Makefile for ConfigFile class

CC = g++
CFLAGS = -Wall -ansi

all: example tester

example: example.cpp ConfigFile.cpp
	$(CC) $(CFLAGS) -o example example.cpp ConfigFile.cpp

tester: tester.cpp ConfigFile.cpp
	$(CC) $(CFLAGS) -o tester tester.cpp ConfigFile.cpp

run: example
	./example | tee example.out

test: tester
	./tester | tee test.out

tidy:
	@ rm -f *.o

clean:  tidy
	@ rm -f example tester example.out test.out core*
