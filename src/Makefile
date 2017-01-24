## Makefile
HOST:=$(shell hostname | cut -d . -f 1)

# Cluster at AEI
ifeq ($(HOST), gpu-01)
BOOST_PATH=/cluster/boost/1.52.0
CUDA_PATH=/cluster/gpu/nvidia/5.5
MPI_LIB_PATH=/cluster/openmpi/SL6/1.5.4/gcc/lib
MPI_INC_PATH=/cluster/openmpi/SL6/1.5.4/gcc/include
endif
# Cluster at AEI
ifeq ($(HOST), pinatubo)
BOOST_PATH=/usr
CUDA_PATH=/usr
MPI_PATH=/usr/lib/openmpi
endif
ifeq ($(HOST), neXus)
BOOST_PATH=/usr
CUDA_PATH=/opt/cuda
MPI_LIB_PATH=/usr/lib/openmpi
MPI_INC_PATH=/usr/include/openmpi
endif

############################# Compilers
CXX=        g++-5
NVCC=       $(CUDA_PATH)/bin/nvcc

############################# FLAGS

# Boost FLAGS
BOOSTFLAGS= -L$(BOOST_PATH)/lib     \
            -I$(BOOST_PATH)/include \
            -lboost_program_options

# Links FLAGS
LINKFLAGS=  -L$(CUDA_PATH)/lib64     \
            -I$(CUDA_PATH)/include   \
            -lcuda                   \
            -lcudart                 \
            -DGPU

# C++ FLAGS
CXXFLAGS=   $(BOOSTFLAGS) \
            -O3           \
            -Wall         \
            -pg           \
            -lgomp        \
            -fopenmp      \
            -pipe         \
            -fstack-protector -Wl,-z,relro -Wl,-z,now -Wformat-security \
            -Wpointer-arith \
            -Wformat-nonliteral \
            -Werror \
            -Wl,-O1 \
            -Wl,--discard-all \
            -Wl,--no-undefined \
            -rdynamic \
            #-fno-exceptions \

# MPI FLAGS
MPIFLAGS=   -L$(MPI_LIB_PATH)             \
            -I$(MPI_INC_PATH)             \
            -pthread                      \
            -lmpi                         \
            -lmpi_cxx                     \
            -ldl                          \
            -D_MPI                        \
            -DOMPI_SKIP_MPICXX

# CUDA FLAGS
#NVFLAGS=   -DGPU -pg -O3 -std=c++11 -arch=sm_20 # -DBOOST_NOINLINE='__attribute__ ((noinline))'
#NVFLAGS=   -DGPU -pg -O3 -std=c++11 -arch=sm_20 -D_MWAITXINTRIN_H_INCLUDED # For gcc > 5.2.1
NVFLAGS=   -DGPU -pg -O3 -arch=sm_20

############################# Objects

OBJS=       include/utils/OptionsParser.o   \
            include/utils/NbodyUtils.o      \
            include/utils/Logger.o          \
            include/Hermite4.o              \
            include/NbodySystem.o

OBJS_CPU=   $(OBJS)                         \
            include/cpu/Hermite4CPU_integration.o \
            include/cpu/Hermite4CPU.o

OBJS_MPI=   $(OBJS)                         \
            include/mpi/Hermite4MPI_integration.o       \
            include/mpi/Hermite4MPI.o       \
            include/mpi/MPIUtils.o

OBJS_GPU=   $(OBJS)                               \
            include/gpu/Hermite4GPU_kernels.cuo   \
            include/gpu/Hermite4GPU_integration.o \
            include/gpu/Hermite4GPU.cuo

############################# Main rules
all: gpu # use `cpu`, `gpu` or `mpi` to set as default

cpu: gravidy-cpu

mpi: CXXFLAGS+=$(MPIFLAGS)
mpi: gravidy-mpi

gpu: CXXFLAGS+=$(LINKFLAGS)
gpu: gravidy-gpu

############################# Rules to generate CUDA object
%.cuo : %.cu
	$(NVCC) $(BOOSTFLAGS) $(NVFLAGS) -c $^ -o $@

############################# Rules to generate binaries

gravidy: gravidy.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

gravidy-cpu: gravidy.cpp $(OBJS_CPU)
	$(CXX) $(CXXFLAGS) $^ -o $@

gravidy-mpi: gravidy.cpp $(OBJS_MPI)
	$(CXX) $(CXXFLAGS) $^ -o $@

gravidy-gpu: gravidy.cpp $(OBJS_GPU)
	$(CXX) $(CXXFLAGS) $^ -o $@

doxy:
	@doxygen -s Doxyfile

clean:
	rm -f gravidy
	rm -f gravidy-gpu
	rm -f gravidy-cpu
	rm -f gravidy-mpi

distclean: clean
	rm -f include/*.o
	rm -f include/cpu/*.o
	rm -f include/utils/*.o
	rm -f include/gpu/*.{o,cuo}
	rm -f include/mpi/*.o
	rm -f gmon.out
