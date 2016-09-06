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

Hermite4MPI::Hermite4MPI(NbodySystem *ns, Logger *logger, NbodyUtils *nu,
                         int rank, int nproc )
                         : Hermite4(ns, logger, nu)
{
    this->rank = rank;
    this->nproc = nproc;
    this->tag = 42;

    // Configuration of the range that the nodes will use in different
    // calculations inside the integration method.
    //
    // Every node process will take care of a N/PROCS_NUMBER particles,
    // so after every node made some procedure with the particles
    // the root node (rank 0) will take care and perform the necesarry
    // synchronization and reduction process.
    this->chunk_size  = std::ceil((float)ns->n / this->nproc);
    this->chunk_begin = this->chunk_size * rank;
    this->chunk_end   = this->chunk_begin + this->chunk_size;
    //printf("Rank %d (%d, %d) %f\n", rank, this->chunk_begin, this->chunk_end, ns->h_r[5].y);

    if (rank < MPI_NUM_SLAVES)
    {
        alloc_slaves_memory(rank);
    }

}

Hermite4MPI::~Hermite4MPI()
{
    // free main and slaves
    free(tmp_f);
    free(ns->h_fout_tmp);

}

void Hermite4MPI::alloc_slaves_memory(int rank)
{

    define_forces_struct(&f_type);
    MPI_Op_create(forces_operation, 1, &f_op);

    tmp_f           = (Forces*)malloc(sizeof(Forces) * ns->n);
    ns->h_fout_tmp  = (Forces*)malloc(sizeof(Forces) * ns->n);

    //#pragma omp parallel for
    for (unsigned int i = 0; i < ns->n; i++)
    {
        tmp_f[i].a[0]  = 0.0;
        tmp_f[i].a[1]  = 0.0;
        tmp_f[i].a[2]  = 0.0;

        tmp_f[i].a1[0] = 0.0;
        tmp_f[i].a1[1] = 0.0;
        tmp_f[i].a1[2] = 0.0;

        ns->h_fout_tmp[i].a[0] = 0.0;
        ns->h_fout_tmp[i].a[1] = 0.0;
        ns->h_fout_tmp[i].a[2] = 0.0;

        ns->h_fout_tmp[i].a1[0] = 0.0;
        ns->h_fout_tmp[i].a1[1] = 0.0;
        ns->h_fout_tmp[i].a1[2] = 0.0;
    }
}

void Hermite4MPI::force_calculation(Predictor pi, Predictor pj, Forces &fi)
{
    double rx = pj.r[0] - pi.r[0];
    double ry = pj.r[1] - pi.r[1];
    double rz = pj.r[2] - pi.r[2];

    double vx = pj.v[0] - pi.v[0];
    double vy = pj.v[1] - pi.v[1];
    double vz = pj.v[2] - pi.v[2];

    double r2     = rx*rx + ry*ry + rz*rz + ns->e2;
    double rinv   = 1.0/sqrt(r2);
    double r2inv  = rinv  * rinv;
    double r3inv  = r2inv * rinv;
    double r5inv  = r2inv * r3inv;
    double mr3inv = r3inv * pj.m;
    double mr5inv = r5inv * pj.m;

    double rv = rx*vx + ry*vy + rz*vz;

    fi.a[0] += (rx * mr3inv);
    fi.a[1] += (ry * mr3inv);
    fi.a[2] += (rz * mr3inv);

    fi.a1[0] += (vx * mr3inv - (3 * rv ) * rx * mr5inv);
    fi.a1[1] += (vy * mr3inv - (3 * rv ) * ry * mr5inv);
    fi.a1[2] += (vz * mr3inv - (3 * rv ) * rz * mr5inv);

}

void Hermite4MPI::init_acc_jrk()
{

    if (rank < MPI_NUM_SLAVES)
    {
        for (unsigned int i = 0; i < ns->n; i++)
        {
            for (unsigned int j = chunk_begin; j < chunk_end; j++)
            {
                if(i == j) continue;
                force_calculation(ns->h_p[i], ns->h_p[j], tmp_f[i]);
            }
        }
        MPI_Allreduce(tmp_f, ns->h_f, ns->n, f_type, f_op, MPI_COMM_WORLD);
    }
}

void Hermite4MPI::update_acc_jrk(unsigned int nact)
{
    ns->gtime.update_ini = omp_get_wtime();
    if (rank < MPI_NUM_SLAVES)
    {

        // Temporal Predictor array with only the active particle information
        for (unsigned int i = 0; i < nact; i++)
        {
            int id = ns->h_move[i];
            ns->h_i[i] = ns->h_p[id];
        }
        MPI_Barrier(MPI_COMM_WORLD);

        // Forces loop with the temporal active particles array and the whole system
        for (unsigned int i = 0; i < nact; i++)
        {
            tmp_f[i].a[0]  = 0.0;
            tmp_f[i].a[1]  = 0.0;
            tmp_f[i].a[2]  = 0.0;
            tmp_f[i].a1[0] = 0.0;
            tmp_f[i].a1[1] = 0.0;
            tmp_f[i].a1[2] = 0.0;

            Predictor pi = ns->h_i[i];

            for (unsigned int j = chunk_begin; j < chunk_end; j++)
            {
                if(ns->h_move[i] == j) continue;
                force_calculation(pi, ns->h_p[j], tmp_f[i]);
            }
        }

        // All the nodes will reduce the forces, having the same results for the
        // new forces.
        MPI_Allreduce(tmp_f, ns->h_fout_tmp, nact, f_type, f_op, MPI_COMM_WORLD);

        for (unsigned int i = 0; i < nact; i++)
        {
            int id = ns->h_move[i];
            ns->h_f[id] = ns->h_fout_tmp[i];
        }
        MPI_Barrier(MPI_COMM_WORLD);

    }
    ns->gtime.update_end += omp_get_wtime() - ns->gtime.update_ini;
}

void Hermite4MPI::predicted_pos_vel(double ITIME)
{

    if (rank == 0)
    {
        ns->gtime.prediction_ini = omp_get_wtime();
        for (unsigned int i = 0; i < ns->n; i++)
        {
            double dt  = ITIME - ns->h_t[i];
            double dt2 = (dt  * dt);
            double dt3 = (dt2 * dt);

            ns->h_p[i].r[0] = (dt3/6 * ns->h_f[i].a1[0]) + (dt2/2 * ns->h_f[i].a[0]) + (dt * ns->h_v[i].x) + ns->h_r[i].x;
            ns->h_p[i].r[1] = (dt3/6 * ns->h_f[i].a1[1]) + (dt2/2 * ns->h_f[i].a[1]) + (dt * ns->h_v[i].y) + ns->h_r[i].y;
            ns->h_p[i].r[2] = (dt3/6 * ns->h_f[i].a1[2]) + (dt2/2 * ns->h_f[i].a[2]) + (dt * ns->h_v[i].z) + ns->h_r[i].z;

            ns->h_p[i].v[0] = (dt2/2 * ns->h_f[i].a1[0]) + (dt * ns->h_f[i].a[0]) + ns->h_v[i].x;
            ns->h_p[i].v[1] = (dt2/2 * ns->h_f[i].a1[1]) + (dt * ns->h_f[i].a[1]) + ns->h_v[i].y;
            ns->h_p[i].v[2] = (dt2/2 * ns->h_f[i].a1[2]) + (dt * ns->h_f[i].a[2]) + ns->h_v[i].z;

            ns->h_p[i].m = ns->h_r[i].w;

        }
        ns->gtime.prediction_end += omp_get_wtime() - ns->gtime.prediction_ini;
    }
    MPI_Bcast(ns->h_p, sizeof(Predictor) * ns->n, MPI_BYTE, 0, MPI_COMM_WORLD);
}

void Hermite4MPI::correction_pos_vel(double ITIME, unsigned int nact)
{
    ns->gtime.correction_ini = omp_get_wtime();
    for (unsigned int k = 0; k < nact; k++)
    {
        int i = ns->h_move[k];

        double dt1 = ns->h_dt[i];
        double dt2 = dt1 * dt1;
        double dt3 = dt2 * dt1;
        double dt4 = dt2 * dt2;
        double dt5 = dt4 * dt1;

        // Acceleration 2nd derivate
        ns->h_a2[i].x = (-6 * (ns->h_old[i].a[0] - ns->h_f[i].a[0] ) - dt1 * (4 * ns->h_old[i].a1[0] + 2 * ns->h_f[i].a1[0]) ) / dt2;
        ns->h_a2[i].y = (-6 * (ns->h_old[i].a[1] - ns->h_f[i].a[1] ) - dt1 * (4 * ns->h_old[i].a1[1] + 2 * ns->h_f[i].a1[1]) ) / dt2;
        ns->h_a2[i].z = (-6 * (ns->h_old[i].a[2] - ns->h_f[i].a[2] ) - dt1 * (4 * ns->h_old[i].a1[2] + 2 * ns->h_f[i].a1[2]) ) / dt2;

        // Acceleration 3rd derivate
        ns->h_a3[i].x = (12 * (ns->h_old[i].a[0] - ns->h_f[i].a[0] ) + 6 * dt1 * (ns->h_old[i].a1[0] + ns->h_f[i].a1[0]) ) / dt3;
        ns->h_a3[i].y = (12 * (ns->h_old[i].a[1] - ns->h_f[i].a[1] ) + 6 * dt1 * (ns->h_old[i].a1[1] + ns->h_f[i].a1[1]) ) / dt3;
        ns->h_a3[i].z = (12 * (ns->h_old[i].a[2] - ns->h_f[i].a[2] ) + 6 * dt1 * (ns->h_old[i].a1[2] + ns->h_f[i].a1[2]) ) / dt3;

        // Correcting position
        ns->h_r[i].x = ns->h_p[i].r[0] + (dt4/24)*ns->h_a2[i].x + (dt5/120)*ns->h_a3[i].x;
        ns->h_r[i].y = ns->h_p[i].r[1] + (dt4/24)*ns->h_a2[i].y + (dt5/120)*ns->h_a3[i].y;
        ns->h_r[i].z = ns->h_p[i].r[2] + (dt4/24)*ns->h_a2[i].z + (dt5/120)*ns->h_a3[i].z;

        // Correcting velocity
        ns->h_v[i].x = ns->h_p[i].v[0] + (dt3/6)*ns->h_a2[i].x + (dt4/24)*ns->h_a3[i].x;
        ns->h_v[i].y = ns->h_p[i].v[1] + (dt3/6)*ns->h_a2[i].y + (dt4/24)*ns->h_a3[i].y;
        ns->h_v[i].z = ns->h_p[i].v[2] + (dt3/6)*ns->h_a2[i].z + (dt4/24)*ns->h_a3[i].z;


        ns->h_t[i] = ITIME;
        double normal_dt  = nu->get_timestep_normal(i, ETA_N);
        normal_dt = nu->normalize_dt(normal_dt, ns->h_dt[i], ns->h_t[i], i);
        ns->h_dt[i] = normal_dt;

    }
    ns->gtime.correction_end += omp_get_wtime() - ns->gtime.correction_ini;
}
