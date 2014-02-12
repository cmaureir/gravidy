#ifndef MPI_UTILS
#define MPI_UTILS
#include "Hermite4MPIGPU.cuh"

void define_forces_struct(MPI_Datatype *ts);
void forces_operation(void *in, void *inout, int *len, MPI_Datatype *type);

#endif // MPI_UTILS
