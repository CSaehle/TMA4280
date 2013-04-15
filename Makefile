all: serial parallel

test: parallel
	./test.sh

kongullparallel:
	module load intel/compilers/11.1.059
	module load openmpi/1.4.3
	mpicc -openmp -o parallel parallel.c fst.o -g -lm

parallel: fst.o parallel.c
	mpicc -fopenmp -o parallel parallel.c fst.o -g -lm

serial: fst.o poisson.c
	gcc -o serial poisson.c fst.o -lm

fst.o: fst.f
	gfortran -o fst.o -c fst.f -g

clean:
	rm -rf tmp serial parallel fst.o
	
remake: | clean all
	
.PHONY: clean remake test
