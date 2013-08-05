#!/bin/bash
#PBS -q gpum
#PBS -N eta_1024
#PBS -l walltime=999:00:00


file=04-nbody-p1024_m1.in
etas=(0.001 0.003 0.005 0.007 0.01 0.03 0.05 0.07 0.1)

cd /user/c/cmaureir/GraviDy/

for e in "${etas[@]}"
do
    (time gpu/./gravidy -i input/$file -t 1 -e $e -o eta_1k)
done

