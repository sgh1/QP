CC=g++
CFLAGS=-c -g -O3 -Wall

all: quantum

quantum: main.o domain.o 
	$(CC) main.o domain.o -o quantum

main.o: main.cpp matrix.h examples.h
	$(CC) $(CFLAGS) main.cpp

domain.o: domain.cpp
	$(CC) $(CFLAGS) domain.cpp

clean:
	rm -rf *o quantum
