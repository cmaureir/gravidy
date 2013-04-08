#!/bin/bash
#PBS -q gpum
#PBS -N eta_65536
#PBS -l walltime=999:00:00


file=10-nbody-p65536_m1.in
etas=(0.001 0.003)

cd /user/c/cmaureir/GraviDy/

for e in "${etas[@]}"
do
    (time gpu/./gravidy -i input/$file -t 1 -e $e -o eta_65k)
done

