all: parallel

parallel: parallel.c
	mpicc -o parallel parallel.c -lm

clean:
	rm -f parallel

remake: | clean all

.PHONY: clean remake
