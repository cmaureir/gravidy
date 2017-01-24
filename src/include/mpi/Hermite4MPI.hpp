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
#ifndef HERMITE4MPI_HPP
#define HERMITE4MPI_HPP
#include "../Hermite4.hpp"

#define MPI_NUM_SLAVES 600

/**
 * Class which implements on the MPI the structure of the Hermite4 scheme.
 *
 * This contains all the implementations of the requirements processes to perform
 * the integration, like the initialization of the forces, the prediction,
 * the correction, and the general integration of the system.
 *
 * Please note that this class uses the variable MPI_NUM_SLAVES to set an upper
 * limit for the amount of slaves that can be use.
 */
class Hermite4MPI : public Hermite4 {
    public:
        Hermite4MPI(NbodySystem *ns, Logger *logger, NbodyUtils *nu,
                    int rank, int nproc);
        ~Hermite4MPI();

        int rank;
        int nproc;
        int tag;
        MPI_Status   status;
        MPI_Datatype f_type;
        MPI_Op       f_op;

        unsigned int chunk_size;
        unsigned int chunk_begin;
        unsigned int chunk_end;
        Forces *tmp_f;

        void alloc_slaves_memory(int rank);
        void force_calculation(Predictor pi, Predictor pj, Forces &fi);
        void init_acc_jrk();
        void update_acc_jrk(unsigned int nact);
        void predicted_pos_vel(double ITIME);
        void correction_pos_vel(double ITIME, unsigned int nact);
        void integration();
};

#endif
