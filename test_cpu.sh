#! /bin/bash
#$ -r n
#$ -l h_rt=30:00:00
#$ -q gpu.q
#$ -pe gpu 1
#$
#$ -N gravidy_cpu_1k
#$ -M cmaurei@aei.mpg.de
#$ -m abe
#$ -o /home/cmaurei/repos/GraviDy/jobs/cpu_test-0001.out
#$ -e /home/cmaurei/repos/GraviDy/jobs/cpu_test-0001.err
#$ -cwd
#$ -R y

module add gpu/NVIDIA_GPU_Computing_SDK gpu/cudatoolkit4.1.28
export LD_LIBRARY_PATH=/home/cmaurei/boost/boost_build/lib:$LD_LIBRARY_PATH

cd /home/cmaurei/repos/GraviDy
./gravidy -i input/04-nbody-p1024_m1.in -t 21 -r cpu
