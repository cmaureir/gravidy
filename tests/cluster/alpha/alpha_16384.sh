#!/bin/bash
#PBS -q gpum
#PBS -N alpha_16384
#PBS -l walltime=999:00:00


file=08-nbody-p16384_m1.in
alphas=(0.01 0.02 0.03 0.04 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)

cd /user/c/cmaureir/GraviDy/

for a in "${alphas[@]}
do
    (time gpu/./gravidy -i input/$file -t 1 -a $a)
done

