all: conway
conway: main.c
	mpicc -o conway main.c -O2 -std=c99

