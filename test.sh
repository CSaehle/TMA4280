#! /usr/bin/env bash

make
mkdir -p tmp
rm -f tmp/serial tmp/parallel
for i in $(seq 3 10); do
    ./serial $((2**i)) >> tmp/serial
    mpirun -n 4 parallel $((2**i)) >> tmp/parallel
done
diff tmp/serial tmp/parallel
