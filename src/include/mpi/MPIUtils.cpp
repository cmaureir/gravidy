#include "MPIUtils.hpp"

void define_forces_struct(MPI_Datatype *ts)
{
    const int count = 3;
    int          blocklens[count];
    MPI_Datatype types[count];
    MPI_Aint     disps[count];

    for (int i = 0; i < count; i++) {
        if (i == 2)
        {
            types[i] = MPI_INT;
            blocklens[i] = 1;
        }
        else
        {
            types[i] = MPI_DOUBLE;
            blocklens[i] = 3;
        }
    }

    disps[0] = offsetof(Forces,a);
    disps[1] = offsetof(Forces,a1);
    disps[2] = offsetof(Forces,nb);

    MPI_Type_create_struct(count, blocklens, disps, types, ts);
    MPI_Type_commit(ts);
}

void forces_operation(void *in, void *inout, int *len, MPI_Datatype *type){

    Forces *invals    = (Forces*)in;
    Forces *inoutvals = (Forces*)inout;

    for (int i = 0; i < *len; i++)
    {
        inoutvals[i].a[0]  += invals[i].a[0];
        inoutvals[i].a[1]  += invals[i].a[1];
        inoutvals[i].a[2]  += invals[i].a[2];

        inoutvals[i].a1[0]  += invals[i].a1[0];
        inoutvals[i].a1[1]  += invals[i].a1[1];
        inoutvals[i].a1[2]  += invals[i].a1[2];

        inoutvals[i].nb += invals[i].nb;
    }

    return;
}
