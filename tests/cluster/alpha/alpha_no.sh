#!/bin/bash
#PBS -q gpum
#PBS -N alpha_no
#PBS -l walltime=999:00:00


cd /user/c/cmaureir/GraviDy/

for i in $(/bin/ls -1 input/ | grep "\.in$")
do
   time cpu/./gravidy -i input/$i -t 1 -o cpu_alpha0.01_$i
done

