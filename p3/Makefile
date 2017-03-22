SEQSRC=$(wildcard *.seq.cc)
MPISRC=$(wildcard *.mpi.cc)
CUDASRC=$(wildcard *.cu)

SEQEXE=$(basename $(SEQSRC))
MPIEXE=$(basename $(MPISRC))
CUDAEXE=$(basename $(CUDASRC))

CXXFLAGS=-Wall -pedantic
MPICXX=mpic++

default: $(SEQEXE) $(MPIEXE)

%.mpi: %.mpi.cc
	$(MPICXX) $(CXXFLAGS) $< -o $@
