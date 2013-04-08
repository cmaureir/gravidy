#!/bin/bash
#PBS -q gpum
#PBS -N eta_new
#PBS -l walltime=999:00:00


file=07-nbody-p8192_m1.in
etas=(0.1)

cd /user/c/cmaureir/GraviDy/

for e in "${etas[@]}"
do
    (time gpu/./gravidy -i input/$file -t 10 -e $e -o 8k)
done

