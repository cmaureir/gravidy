#!/bin/bash
#PBS -q gpum
#PBS -N lagrange_radii_1024
#PBS -l walltime=999:00:00


file=04-nbody-p1024_m1.in

cd /user/c/cmaureir/GraviDy/

(time gpu/./gravidy -i input/$file -t 1000 -o 1024_1000t)

