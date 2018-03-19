# Makefile
HOST:=$(shell hostname | cut -d . -f 1)

ifeq ($(HOST), gpu-01)
    BOOST_PATH=/cluster/boost/1.52.0
    CUDA_PATH=/cluster/gpu/nvidia/5.5
    MPI_LIB_PATH=/cluster/openmpi/SL6/1.5.4/gcc/lib
    MPI_INC_PATH=/cluster/openmpi/SL6/1.5.4/gcc/include
else ifeq ($(HOST), Gateway)
    BOOST_PATH=/usr
    CUDA_PATH=/opt/cuda
    MPI_LIB_PATH=/usr/lib/openmpi
    MPI_INC_PATH=/usr/include/openmpi
else
    BOOST_PATH=
    CUDA_PATH=
    MPI_LIB_PATH=
    MPI_INC_PATH=
endif

# Compilers
CXX=  g++
NVCC= $(CUDA_PATH)/bin/nvcc

# FLAGS

# Boost FLAGS
BOOSTFLAGS  = -L$(BOOST_PATH)/lib
BOOSTFLAGS += -I$(BOOST_PATH)/include
BOOSTFLAGS += -lboost_program_options

# Links FLAGS
LINKFLAGS  = -L$(CUDA_PATH)/lib64
LINKFLAGS += -I$(CUDA_PATH)/include
LINKFLAGS += -lcuda
LINKFLAGS += -lcudart
LINKFLAGS += -DGPU

# C++ FLAGS
CXXFLAGS  = $(BOOSTFLAGS)
CXXFLAGS += -O3 # Optimisation
CXXFLAGS += -Wall # Warning messages
CXXFLAGS += -fopenmp # Enable OpenMP
CXXFLAGS += -pipe # Pipes over temp files in the process
CXXFLAGS += -fstack-protector # Emit extra code to check buffer overflow
CXXFLAGS += -Wl,-z,relro # RO relocation table area in the final ELF
CXXFLAGS += -Wl,-z,now # Prevent GOT overwrite attacks
CXXFLAGS += -Wformat-security # Prevent formats vulnerabilities
CXXFLAGS += -Wpointer-arith # Prevent void size
CXXFLAGS += -Wformat-nonliteral # Prevent non literal format string
CXXFLAGS += -Wl,-O1 # Optimisation for the linker
CXXFLAGS += -Wl,--discard-all # Discard local symbols
CXXFLAGS += -Wl,--no-undefined # Prevent missing/unresolved symbols
CXXFLAGS += -rdynamic # Exported symbols in dynamics symbol table
#CXXFLAGS += -pg # For gprof
#CXXFLAGS += -Werror # Warnings will be errors

# MPI FLAGS
MPIFLAGS  = -L$(MPI_LIB_PATH)
MPIFLAGS += -I$(MPI_INC_PATH)
MPIFLAGS += -pthread
MPIFLAGS += -lmpi
MPIFLAGS += -lmpi_cxx
MPIFLAGS += -ldl
MPIFLAGS += -D_MPI
MPIFLAGS += -DOMPI_SKIP_MPICXX

# CUDA FLAGS
NVFLAGS  = -DGPU
NVFLAGS += -pg
NVFLAGS += -O3
#NVFLAGS=-std=c++11 -arch=sm_20 # -DBOOST_NOINLINE='__attribute__ ((noinline))'
#NVFLAGS=-std=c++11 -arch=sm_20 -D_MWAITXINTRIN_H_INCLUDED # For gcc > 5.2.1

# General Objects
OBJS  = include/utils/OptionsParser.o
OBJS += include/utils/NbodyUtils.o
OBJS += include/utils/Logger.o
OBJS += include/Hermite4.o
OBJS += include/NbodySystem.o

# CPU objects
OBJS_CPU  = $(OBJS)
OBJS_CPU += include/cpu/Hermite4CPU_integration.o
OBJS_CPU += include/cpu/Hermite4CPU.o

# MPI objects
OBJS_MPI  = $(OBJS)
OBJS_MPI += include/mpi/Hermite4MPI_integration.o
OBJS_MPI += include/mpi/Hermite4MPI.o
OBJS_MPI += include/mpi/MPIUtils.o

# GPU objects
OBJS_GPU  = $(OBJS)
OBJS_GPU += include/gpu/Hermite4GPU_kernels.cuo
OBJS_GPU += include/gpu/Hermite4GPU_integration.o
OBJS_GPU += include/gpu/Hermite4GPU.cuo

# Main rules
.PHONY: cpu pn mpi gpu doxy clean distclean

all: cpu

pn: print_pn
pn: gravidy-cpu

cpu: gravidy-cpu

mpi: CXXFLAGS+=$(MPIFLAGS)
mpi: gravidy-mpi

gpu: CXXFLAGS+=$(LINKFLAGS)
gpu: gravidy-gpu

print_pn:
	@echo "PN implementation is enabled in CPU version"

# Rules to generate CUDA object
%.cuo : %.cu
	$(NVCC) $(BOOSTFLAGS) $(NVFLAGS) -c $^ -o $@

# Rules to generate binaries
gravidy: gravidy.cpp
	$(CXX) $^ -o $@ $(CXXFLAGS)

gravidy-cpu: gravidy.cpp $(OBJS_CPU)
	$(CXX) $^ -o $@ $(CXXFLAGS)

gravidy-mpi: gravidy.cpp $(OBJS_MPI)
	$(CXX) $^ -o $@ $(CXXFLAGS)

gravidy-gpu: gravidy.cpp $(OBJS_GPU)
	$(CXX) $^ -o $@ $(CXXFLAGS)

# Generate documentation
doxy:
	@doxygen -s Doxyfile

# Clean
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
