all: serial parallel

test: parallel
	./test.sh

parallel: fst.o parallel.c
	mpicc -fopenmp -o parallel parallel.c fst.o -lm

serial: fst.o poisson.c
	gcc -o serial poisson.c fst.o -lm

fst.o: fst.f
	gfortran -o fst.o -c fst.f

clean:
	rm -rf tmp serial parallel fst.o
	
remake: | clean all
	
.PHONY: clean remake test
