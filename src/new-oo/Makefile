# Makefile

CXX=        g++
CXXFLAGS=   -lboost_program_options  \
            -O3                      \
            -m64                     \
            -Wall                    \
            -pg                      \
            -lgomp                   \
            -fopenmp

CUDA_PATH=/opt/cuda
NVCC=       $(CUDA_PATH)/bin/nvcc

NVFLAGS=    --pre-include gcc4.7_nvcc_fix.h -pg -O3 -m64 -arch=sm_20 -DGPU

LINKFLAGS=  -L$(CUDA_PATH)/lib64     \
            -I$(CUDA_PATH)/include   \
            -lcuda                   \
            -lcudart                 \
            -DGPU

OBJS_UTILS= include/utils/extra_utils.o   \
            include/utils/OptionsParser.o \
            include/utils/NbodyUtils.o    \
            include/utils/Logger.o

OBJS= include/NbodySystem.o             \
      include/Hermite4.o                \
      include/NbodySystem_integration.o \
      include/NbodySystem_memory_host.o

OBJS_CPU= $(OBJS) $(OBJS_UTILS) \
          include/Hermite4CPU.o \

OBJS_GPU= $(OBJS_CPU)                               \
          include/gpu/NbodySystem_memory_device.cuo \
          include/gpu/NbodySystem_integration.cuo \
          include/gpu/Hermite4GPU.cuo



cpu: gravidy-cpu
gpu: CXXFLAGS+=$(LINKFLAGS)
gpu: gravidy-gpu

%.cuo : %.cu
	$(NVCC) $(NVFLAGS) -c $^ -o $@

gravidy-cpu: gravidy.cpp $(OBJS_CPU)
	$(CXX) $(CXXFLAGS) $^ -o $@

gravidy-gpu: gravidy.cpp $(OBJS_GPU)
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f gravidy
#	rm -f gravidy-gpu
distclean: clean
	rm -f include/*.{o,cuo}
	rm -f include/utils/*.{o,cuo}
	rm -f include/gpu/*.{o,cuo}
	rm -f gmon.out
