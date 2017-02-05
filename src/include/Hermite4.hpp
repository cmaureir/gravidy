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
#ifndef HERMITE4_HPP
#define HERMITE4_HPP
#include "utils/Logger.hpp"
#include "utils/NbodyUtils.hpp"
#include <limits>

/**
 * General class that define the structure to be follow by any implementation
 * of the integrator (CPU, MPI or GPU).
 * The methods on this class must be implemented to have a proper behaviour of the
 * integrator.
 */
class Hermite4 {
    public:
        Hermite4(NbodySystem *ns, Logger *logger, NbodyUtils *nu);
        ~Hermite4();

        /** NbodySystem object reference */
        NbodySystem *ns;
        /** Logger object reference */
        Logger      *logger;
        /** NbodyUtils object reference */
        NbodyUtils  *nu;

        unsigned int  find_particles_to_move(double ITIME);
        void next_integration_time(double &ATIME);
        void init_dt(double &ATIME, float ETA, double ITIME);
        void save_old_acc_jrk(unsigned int nact);
        void alloc_arrays_host();
        void free_arrays_host();
        void init_data();

        /* Virtual methods to be implemented by the different versions */
        /** Integration virtual method to be implemented */
        virtual void integration() {}
        /** Prediction virtual method to be implemented */
        virtual void predicted_pos_vel(double itime, double *t, double4 *r,
                                       double4 *v, Forces *f, Predictor *p) {}
        /** Correction virtual method to be implemented */
        virtual void correction_pos_vel(double itime, unsigned int nact, double *dt,
                                        double *t, unsigned int *move, Predictor *p,
                                        Forces *f, Forces *old, double3 *a2,
                                        double3 *a3, double4 *r, double4 *v) {}
        /** Force virtual method to be implemented */
        virtual void force_calculation(const Predictor &pi, const Predictor &pj, Forces &fi) {}
        /** Force initialization virtual method to be implemented */
        virtual void init_acc_jrk(Predictor *p, Forces *f) {}
        /** Force update virtual method to be implemented */
        virtual void update_acc_jrk(unsigned int nact, unsigned int *move,
                                    Predictor *p, Forces *f) {}

};

#endif // HERMITE4_HPP
