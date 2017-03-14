CC = mpicc
SRC = $(wildcard *.c)
EXE = $(basename $(SRC))

CFLAGS = -lrt -fopenmp

default: $(EXE)

run: $(EXE)
	/usr/lib64/openmpi/bin/mpiexec -N 3 $<
