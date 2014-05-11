# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
CC=g++
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.
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
	#rm files/* 
	#cp ../files/*.bin files
