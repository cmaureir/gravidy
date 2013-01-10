CUDA_PATH=/opt/cuda

CC=         g++
NVCC=       nvcc
NVFLAGS=    -O3 -m64 -G -g -arch=sm_20
#NVFLAGS=    -O3 -m64 -G -g -arch=sm_20 --ptxas-options="-v -dlcm=cg"
CXXFLAGS=   -O3 -m64 -Wall -lgomp -pg -fopenmp
CPPFLAGS=   -L$(CUDA_PATH)/lib64 -I$(CUDA_PATH)/include -lcuda -lcudart -lboost_program_options

OBJS_CPU=   include/options_parser.o \
			include/file_utils.o     \
			include/extra_utils.o    \
			include/equilibrium.o    \
			include/memory_cpu.o         \
			include/dynamics_cpu.o   \

OBJS_GPU=   include/memory_gpu.cuo \
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
