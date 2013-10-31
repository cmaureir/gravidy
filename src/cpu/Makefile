CC=        g++
CXXFLAGS=  -O3 -Wall -lgomp -fopenmp -pg -m64

CPPFLAGS= -lboost_program_options

OBJS_CPU= include/options_parser.o   \
	  include/file_utils.o           \
	  include/extra_utils.o          \
	  include/equilibrium.o          \
	  include/memory_cpu.o               \
	  include/dynamics_cpu.o         \
	  include/hermite.o

OBJS = $(OBJS_CPU)

all: gravidy

gravidy: gravidy.cpp $(OBJS)
clean:
	rm -rf gravidy
distclean: clean
	rm -rf include/*.o
	rm -rf gmon.out
