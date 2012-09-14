#CXXFLAGS=-O2 -pg -Wall -lboost_program_options -lgomp
CC=		g++
NVCC=		nvcc
NVFLAGS=	-O3 -m64 -arch=sm_20
CUDA_PATH=	/usr/local/cuda
BOOST_PATH=	/home/cmaurei/boost/boost_build
CXXFLAGS=	-O3 -m64 -Wall -lgomp -pg -fopenmp
CPPFLAGS=	-I $(BOOST_PATH)/include	\
		-I $(CUDA_PATH)/include        	\
		-L $(BOOST_PATH)/lib -lboost_program_options \
		-L $(CUDA_PATH)/lib64 -lcudart 	

OBJS_CPU= include/options_parser.o       \
	  include/file_utils.o           \
	  include/extra_utils.o          \
	  include/equilibrium.o          \
	  include/memory.o               \
	  include/kepler.o               \
          include/dynamics_cpu.o         \
	  include/hermite.o

OBJS_GPU= include/dynamics_gpu_kernels.o \
	  include/dynamics_gpu.o

OBJS = $(OBJS_GPU) $(OBJS_CPU)

include/dynamics_gpu_kernels.o: include/dynamics_gpu_kernels.cu
	$(NVCC) $(NVFLAGS) -c $^ -o $@
include/dynamics_gpu.o: include/dynamics_gpu.cu
	$(NVCC) $(NVFLAGS) -c $^ -o $@

gravidy: gravidy.cpp $(OBJS)
all: gravidy
clean:
	rm -rf gravidy include/dynamics_gpu.o
distclean: clean
	rm -rf include/*.o
