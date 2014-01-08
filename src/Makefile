## Makefile

############################# Paths
#BOOST_PATH=/cluster/boost/1.52.0
#CUDA_PATH=/cluster/gpu/nvidia/5.5
CUDA_PATH=/opt/cuda

############################# Compilers
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
            -O3           \
            -Wall         \
            -pg           \
            -lgomp        \
            -fopenmp

# CUDA FLAGS
NVFLAGS=   -pg -O3 -arch=sm_20 -DGPU

############################# Objects

OBJS=       include/utils/OptionsParser.o   \
            include/NbodySystem.o           \
            include/utils/NbodyUtils.o      \
            include/utils/Logger.o          \
            include/Hermite4.o              \

OBJS_CPU=   $(OBJS)                         \
            include/cpu/Hermite4CPU.o       \

OBJS_GPU=   $(OBJS)                         \
            include/gpu/Hermite4GPU.cuo

OBJS_KEPLER=$(OBJS)                         \
            include/kepler/Hermite4Kepler.o

############################# Main rules

cpu: gravidy-cpu
gpu: CXXFLAGS+=$(LINKFLAGS)
gpu: gravidy-gpu
kepler: CXXFLAGS+=-DKEPLER
kepler: gravidy-kepler

############################# Rules to generate cuda object
%.cuo : %.cu
	$(NVCC) $(BOOSTFLAGS) $(NVFLAGS) -c $^ -o $@

############################# Rules to generate binaries

gravidy-cpu: gravidy.cpp $(OBJS_CPU)
	$(CXX) $(CXXFLAGS) $^ -o $@

gravidy-gpu: gravidy.cpp $(OBJS_GPU)
	$(CXX) $(CXXFLAGS) $^ -o $@

gravidy-kepler: gravidy.cpp $(OBJS_KEPLER)
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f gravidy-gpu
	rm -f gravidy-cpu
	rm -f gravidy-kepler

distclean: clean
	rm -f include/*.o
	rm -f include/utils/*.o
	rm -f include/gpu/*.{o,cuo}
	rm -f include/kepler/*.o
	rm -f gmon.out
