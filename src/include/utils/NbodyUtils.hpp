#ifndef NBODYUTILS_HPP
#define NBODYUTILS_HPP
#include "../common.hpp"
#include "../NbodySystem.hpp"


class NbodyUtils {
    public:
        NbodyUtils(NbodySystem *ns, float ratio);
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
        double get_timestep_normal(int i, float ETA);
        double normalize_dt(double new_dt, double old_dt, double t, int i);
        double get_timestep_central(int i);

        double get_potential();
        double get_kinetic();
        double get_energy(double ext = 0);

        double3 get_center_of_density();
};

#endif
