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

        void next_integration_time(double &ATIME);
        void save_old_acc_jrk(unsigned int nact);
        unsigned int  find_particles_to_move(double ITIME);
        void init_dt(double &ATIME, float ETA, double ITIME);

    private:
        void init_data();
        void alloc_arrays_host();
        void free_arrays_host();

        /* Virtual methods to be implemented by the different versions */
        /** Integration virtual method to be implemented */
        virtual void integration() {}
        /** Force virtual method to be implemented */
        virtual void force_calculation(const Predictor &,
            const Predictor &, Forces &) {}

        /* On top of the two methods above, every implementation provides its
         * own prediction (predicted_pos_vel), correction (correction_pos_vel),
         * force initialization (init_acc_jrk) and force update
         * (update_acc_jrk).  They are not declared here because each version
         * reads the particle arrays straight from the NbodySystem reference,
         * so their argument lists differ between the CPU, MPI and GPU
         * implementations. */
};

#endif // HERMITE4_HPP
