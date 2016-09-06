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
#ifndef NBODYUTILS_HPP
#define NBODYUTILS_HPP
#include "../common.hpp"
#include "../NbodySystem.hpp"


class NbodyUtils {
    public:
        NbodyUtils(NbodySystem *ns);
        ~NbodyUtils();

        NbodySystem *ns;
        // Variables
        float ratio;
        double3 cod;

        /** Radii array related to the center of density, and related to the
         * percentage distribution that we want to obtain */
        std::vector<Distance> radii;

        /** values of the radii of the different layers */
        std::vector<double> layers_radii;

        void nbody_attributes();
        void get_radii();
        void get_layers();
        void lagrange_radii();
        void core_radius_and_density();

        double get_virial_radius(double energy);
        double get_close_encounter_radius(double r_virial);
        double get_close_encounter_timestep(double r_cl);
        double get_core_radius();
        double get_half_mass_relaxation_time();
        double get_crossing_time(double r_virial);
        double get_halfmass_radius();
        double get_magnitude(double x, double y, double z);
        double get_timestep_normal(unsigned int i, float ETA);
        double normalize_dt(double new_dt, double old_dt, double t, unsigned int i);
        double get_timestep_central(unsigned int i);

        double get_potential();
        double get_kinetic();
        double get_energy(double ext = 0);

        double3 get_center_of_density();
};

#endif
