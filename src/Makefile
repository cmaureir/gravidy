## Makefile
HOST:=$(shell hostname | cut -d . -f 1)

# Cluster at AEI
ifeq ($(HOST), gpu-01)
BOOST_PATH=/cluster/boost/1.52.0
CUDA_PATH=/cluster/gpu/nvidia/5.5
MPI_PATH=/cluster/openmpi/SL6/1.5.4/gcc
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
MPI_PATH=/usr/lib/openmpi
endif

############################# Paths
#CUDA_PATH=/usr

############################# Compilers
#CXX=        vtc++
CXX=        g++
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
            -std=c++11    \
            -O3           \
            -Wall         \
            -pg           \
            -lgomp        \
            -fopenmp

# MPI FLAGS
MPIFLAGS=   -L$(MPI_PATH)                 \
            -L$(MPI_PATH)/lib             \
            -I$(MPI_PATH)/include         \
            -pthread                      \
            -lmpi                         \
            -lmpi_cxx                     \
            -ldl                          \
            -D_MPI                        \
            -DOMPI_SKIP_MPICXX

# CUDA FLAGS
NVFLAGS=   -DGPU -pg -O3 -arch=sm_20 -DBOOST_NOINLINE='__attribute__ ((noinline))'

############################# Objects

OBJS=       include/utils/OptionsParser.o   \
            include/NbodySystem.o           \
            include/utils/NbodyUtils.o      \
            include/utils/Logger.o          \
            include/Hermite4.o              \

OBJS_CPU=   $(OBJS)                         \
            include/cpu/Hermite4CPU.o       \

OBJS_MPI=   $(OBJS)                         \
            include/mpi/Hermite4MPI.o       \
            include/mpi/MPIUtils.o

OBJS_MPIGPU=   $(OBJS)                        \
            include/mpigpu/Hermite4MPIGPU.cuo \
            include/mpigpu/MPIUtils.o

OBJS_GPU=   $(OBJS)                         \
            include/gpu/Hermite4GPU.cuo

OBJS_KEPLER=$(OBJS)                         \
            include/kepler/Hermite4Kepler.o

############################# Main rules

all: gpu
cpu: gravidy-cpu
mpigpu: CXXFLAGS+=-DMPIGPU
mpigpu: NVFLAGS+=-DMPIGPU
mpigpu: NVFLAGS+=-L$(MPI_PATH)/lib -I$(MPI_PATH)/include
mpigpu: CXXFLAGS+=$(MPIFLAGS)
mpigpu: CXXFLAGS+=$(LINKFLAGS)
mpigpu: CXXFLAGS:=$(filter-out -D_MPI,$(CXXFLAGS))
mpigpu: CXXFLAGS:=$(filter-out -DGPU,$(CXXFLAGS))
mpigpu: gravidy

mpi: CXXFLAGS+=$(MPIFLAGS)
mpi: gravidy-mpi

gpu: CXXFLAGS+=$(LINKFLAGS)
gpu: gravidy-gpu

kepler: CXXFLAGS+=-DKEPLER
kepler: gravidy-kepler

############################# Rules to generate cuda object
%.cuo : %.cu
	$(NVCC) $(BOOSTFLAGS) $(NVFLAGS) -c $^ -o $@

############################# Rules to generate binaries

gravidy: gravidy.cpp $(OBJS_MPIGPU)
	$(CXX) $(CXXFLAGS) $^ -o $@

gravidy-cpu: gravidy.cpp $(OBJS_CPU)
	$(CXX) $(CXXFLAGS) $^ -o $@

gravidy-mpi: gravidy.cpp $(OBJS_MPI)
	$(CXX) $(CXXFLAGS) $^ -o $@

gravidy-gpu: gravidy.cpp $(OBJS_GPU)
	$(CXX) $(CXXFLAGS) $^ -o $@

gravidy-kepler: gravidy.cpp $(OBJS_KEPLER)
	$(CXX) $(CXXFLAGS) $^ -o $@



clean:
	rm -f gravidy
	rm -f gravidy-gpu
	rm -f gravidy-cpu
	rm -f gravidy-kepler
	rm -f gravidy-mpi

distclean: clean
	rm -f include/*.o
	rm -f include/utils/*.o
	rm -f include/gpu/*.{o,cuo}
	rm -f include/kepler/*.o
	rm -f include/mpi/*.o
	rm -f include/mpigpu/*.{o,cuo}
	rm -f gmon.out
