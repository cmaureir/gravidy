#!/bin/bash
#PBS -q gpum
#PBS -N lagrange_radii_4096
#PBS -l walltime=999:00:00


file=06-nbody-p4096_m1.in

cd /user/c/cmaureir/GraviDy/

(time gpu/./gravidy -i input/$file -t 1000 -o 4096_1000t)

