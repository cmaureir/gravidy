/*
 * Copyright (c) 2016
 *
 * Cristián Maureira-Fredes <cmaureirafredes@gmail.com>
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. The name of the author may not be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#include "MPIUtils.hpp"

void define_forces_struct(MPI_Datatype *ts)
{
    const int count = 2;
    int          blocklens[count];
    MPI_Datatype types[count];
    MPI_Aint     disps[count];

    for (int i = 0; i < count; i++) {
        types[i] = MPI_DOUBLE;
        blocklens[i] = 3;
    }

    disps[0] = offsetof(Forces,a);
    disps[1] = offsetof(Forces,a1);

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
    }

    return;
}
