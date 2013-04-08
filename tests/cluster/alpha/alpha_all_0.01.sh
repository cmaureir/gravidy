#!/bin/bash
#PBS -q gpum
#PBS -N alpha_all
#PBS -l walltime=999:00:00


cd /user/c/cmaureir/GraviDy/

for i in $(/bin/ls -1 input/ | grep "\.in$")
do
   time gpu/./gravidy -i input/$i -t 1 -a 0.01 -o alpha0.01_$i
done

