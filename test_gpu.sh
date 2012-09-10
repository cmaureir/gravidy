#! /bin/bash
#$ -r n
#$ -l h_rt=3:00:00
#$ -q gpu.q
#$ -pe gpu 1
#$
#$ -N gravidy_gpu_8k
#$ -M cmaurei@aei.mpg.de
#$ -m abe
#$ -o /home/cmaurei/repos/GraviDy/jobs/gpu_test-0001.out
#$ -e /home/cmaurei/repos/GraviDy/jobs/gpu_test-0001.err
#$ -cwd
#$ -R y

module add gpu/NVIDIA_GPU_Computing_SDK gpu/cudatoolkit4.1.28
export LD_LIBRARY_PATH=/home/cmaurei/boost/boost_build/lib:$LD_LIBRARY_PATH

cd /home/cmaurei/repos/GraviDy
./gravidy -i input/07-nbody-p8192_m8192 -t 0.01 -r gpu
