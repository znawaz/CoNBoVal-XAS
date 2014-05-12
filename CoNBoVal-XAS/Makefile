# the compiler to use.
CC=g++
#Compiler options
CFLAGS=-c -Wall

all: CoNBoVal-XAS

CoNBoVal-XAS: point.o Main.o
	$(CC) point.o Main.o -L$(ANN_HOME)/lib -lANN -o CoNBoVal-XAS 


Main.o: Main.cpp
	$(CC) $(CFLAGS) -I$(ANN_HOME)/include Main.cpp

point.o: point.cpp
	$(CC) $(CFLAGS) point.cpp



clean:
	rm *_bondVal.txt *.inp CoNBoVal-XAS
