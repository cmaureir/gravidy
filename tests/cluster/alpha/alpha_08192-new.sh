#!/bin/bash
#PBS -q gpum
#PBS -N alpha_8192-new
#PBS -l walltime=999:00:00


file=07-nbody-p8192_m1.in
alphas=(0.6 0.65 0.7 0.75 0.8)

cd /user/c/cmaureir/GraviDy/

for a in "${alphas[@]}"
do
    (time gpu/./gravidy -i input/$file -t 1 -a $a -o new-alpha-8k)
done

