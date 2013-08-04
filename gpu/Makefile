CUDA_PATH=/opt/cuda
CC=         g++
NVCC=       $(CUDA_PATH)/bin/nvcc
NVFLAGS=    --pre-include gcc4.7_nvcc_fix.h -pg -O3 -m64 -arch=sm_20
CXXFLAGS=   -L$(CUDA_PATH)/lib64     \
            -I$(CUDA_PATH)/include   \
            -lcuda                   \
            -lcudart                 \
            -lboost_program_options  \
            -O3                      \
            -m64                     \
            -Wall                    \
            -pg                      \
            #-lgomp                   \
            #-fopenmp

OBJS_CPU=   include/options_parser.o \
			include/file_utils.o     \
			include/extra_utils.o    \
			include/equilibrium.o    \
			include/dynamics_cpu.o
#			include/memory_cpu.o         \

OBJS_GPU=   include/memory_gpu.cuo \
			include/memory_cpu.cuo \
			include/dynamics_gpu_kernels.cuo \
			include/dynamics_gpu.cuo \
            include/hermite.cuo

OBJS=       $(OBJS_GPU) $(OBJS_CPU)

all: gravidy

%.cuo : %.cu
	$(NVCC) $(NVFLAGS) -c $^ -o $@

gravidy: gravidy.cpp $(OBJS)
clean:
	rm -rf gravidy
distclean: clean
	rm -rf include/*.o
	rm -rf include/*.cuo
