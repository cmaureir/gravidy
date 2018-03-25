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
#include "Hermite4.hpp"

/** Constructor in charge of calling the allocation and
 * initialization methods */
Hermite4::Hermite4(NbodySystem *ns, Logger *logger, NbodyUtils *nu)
{
    this->ns     = ns;
    this->logger = logger;
    this->nu     = nu;

    alloc_arrays_host();
    init_data();
}

/** Destructor in charge of calling the method that free the memory
 */
Hermite4::~Hermite4()
{
    free_arrays_host();
}

/** Method in charge of finding all the particles that need to be
 * updated on the following integration step.
 * Since we are using DP, we base the comparison between times and
 * timesteps using the machine epsilon, to avoid overflows.
 *
 */
unsigned int Hermite4::find_particles_to_move(double ITIME)
{

    unsigned int j = 0;
    for (unsigned int i = 0; i < ns->n; i++)
    {
        ns->h_move[i] = -1;

        // Getting larger mass of the system
        if (ns->h_r[i].w > ns->max_mass)
            ns->max_mass = ns->h_r[i].w;

        double tmp_time = ns->h_t[i] + ns->h_dt[i];
        if(std::fabs(ITIME - tmp_time) <
            2*std::numeric_limits<double>::epsilon())
        {
            ns->h_move[j] = i;
            j++;
        }
    }
    return j;
}

/** Methods that look for the next integration time, looking for the
 * minimum combination of timestep and current time of all the
 * particles of the system
 */
void Hermite4::next_integration_time(double &CTIME)
{
    // Initial number as the maximum
    CTIME = ns->h_t[0] + ns->h_dt[0];
    for (unsigned int i = 0 + 1; i < ns->n; i++)
    {
        double time = ns->h_t[i] + ns->h_dt[i];
        if(time < CTIME)
        {
            CTIME = time;
        }
    }
}

/** Initialization of the timesteps of the system.
 * This method is based on the  Aarseth initial timestep definition
 * and follows a Gaussian distribution.
 */
void Hermite4::init_dt(double &CTIME, float ETA, double ITIME)
{
    // Aarseth initial timestep
    // dt_{i} = ETA_S * sqrt( (|a|) / (|j|) )
    double tmp_dt;
    for (unsigned int i = 0; i < ns->n; i++)
    {
        Forces f = ns->h_f[i];

        double a2 = nu->get_magnitude(f.a[0], f.a[1], f.a[2]);
        double j2 = nu->get_magnitude(f.a1[0], f.a1[1], f.a1[2]);

        tmp_dt = ETA * (a2/j2);

        // Adjusting to block timesteps
        // to the nearest-lower power of two
        int exp = (int)(std::ceil(log(tmp_dt)/log(2.0))-1);
        tmp_dt = pow(2,exp);
        //2 << (exp - 1);

        if (tmp_dt < D_TIME_MIN)
            tmp_dt = D_TIME_MIN;
        else if (tmp_dt > D_TIME_MAX)
            tmp_dt = D_TIME_MAX;

        ns->h_dt[i] = tmp_dt;
        ns->h_t[i] = ITIME;

        // Obtaining the first integration time
        if(tmp_dt < CTIME)
            CTIME = tmp_dt;
    }
}

/** Method in charge of saving the old values of the acceleration and
 * its first derivative to be use in the Corrector integration step
 */
void Hermite4::save_old_acc_jrk(unsigned int nact)
{
    #pragma omp parallel for
    for (unsigned int k = 0; k < nact; k++)
    {
        unsigned int i = ns->h_move[k];
        ns->h_old[i] = ns->h_f[i];
    }
}

/** Method in charge of the memory allocation of all the data
 * structures
 */
void Hermite4::alloc_arrays_host()
{
    ns->h_f    = new Forces[ns->n];
    ns->h_a2   = new double3[ns->n];
    ns->h_a3   = new double3[ns->n];
    ns->h_old  = new Forces[ns->n];
    ns->h_t    = new double[ns->n];
    ns->h_dt   = new double[ns->n];
    ns->h_move = new unsigned int[ns->n];
    ns->h_p    = new Predictor[ns->n];
    ns->h_i    = new Predictor[ns->n];
    ns->h_ekin = new double[ns->n];
    ns->h_epot = new double[ns->n];
}

/** Method in charge of the memory deallocation of all the data
 * structures
 */
void Hermite4::free_arrays_host()
{
    delete [] ns->h_f;
    delete [] ns->h_a2;
    delete [] ns->h_a3;
    delete [] ns->h_old;
    delete [] ns->h_t;
    delete [] ns->h_dt;
    delete [] ns->h_move;
    delete [] ns->h_p;
    delete [] ns->h_i;
    delete [] ns->h_ekin;
    delete [] ns->h_epot;
}

/** Method that initialize all the data structures that will be used
 * in the integration of hte system.
 */
void Hermite4::init_data()
{
    double3 empty = {0};
    ns->m_g = 0.0;
    for (unsigned int i = 0; i < ns->n; i++)
    {
        double mass = ns->h_r[i].w;

        ns->h_p[i].m    = mass;

        ns->h_p[i].r[0] = ns->h_r[i].x;
        ns->h_p[i].r[1] = ns->h_r[i].y;
        ns->h_p[i].r[2] = ns->h_r[i].z;

        ns->h_p[i].v[0] = ns->h_v[i].x;
        ns->h_p[i].v[1] = ns->h_v[i].y;
        ns->h_p[i].v[2] = ns->h_v[i].z;

        ns->h_f[i].a[0] = 0.0;
        ns->h_f[i].a[1] = 0.0;
        ns->h_f[i].a[2] = 0.0;

        ns->h_f[i].a1[0] = 0.0;
        ns->h_f[i].a1[1] = 0.0;
        ns->h_f[i].a1[2] = 0.0;

        ns->h_a2[i] = empty;
        ns->h_a3[i] = empty;

        ns->h_old[i].a[0] = 0.0;
        ns->h_old[i].a[1] = 0.0;
        ns->h_old[i].a[2] = 0.0;

        ns->h_old[i].a1[0] = 0.0;
        ns->h_old[i].a1[1] = 0.0;
        ns->h_old[i].a1[2] = 0.0;

        ns->h_t[i] = 0.0;
        ns->h_dt[i] = 0.0;

        ns->h_move[i] = 0;

        // Heaviest star
        if (mass > ns->m_g)
            ns->m_g = mass;
    }
}
