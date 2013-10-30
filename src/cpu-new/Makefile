CC=        g++
CXXFLAGS=  -O3 -Wall -lgomp -fopenmp -pg

CPPFLAGS= -lboost_program_options

OBJS_CPU= include/NbodySystem.o \
          include/NbodySystem_integration.o \
          include/NbodySystem_memory.o \
          include/OptionsParser.o \
          include/Hermite4.o \
          include/extra_utils.o
OBJS = $(OBJS_CPU)

all: gravidy

gravidy: gravidy.cpp $(OBJS)
clean:
	rm -rf gravidy
distclean: clean
	rm -rf include/*.o
	rm -rf gmon.out
