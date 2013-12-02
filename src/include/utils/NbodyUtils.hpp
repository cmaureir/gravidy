#ifndef NBODYUTILS_HPP
#define NBODYUTILS_HPP
#include "../common.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>

struct Distance
{
    int index;
    double value;
    bool operator<(const Distance& a) const
    {
        return value < a.value;
    }
};

class NbodyUtils {
    private:
        double energy;

    public:
        NbodyUtils(int n, double4 *r, float m, float ratio);
        ~NbodyUtils();

        // Variables
        int n;
        float total_mass;
        float ratio;
        double3 cod;

        double4 *r;
        /** Radii array related to the center of density, and related to the
         * percentage distribution that we want to obtain */
        std::vector<Distance> radii;

        /** values of the radii of the different layers */
        std::vector<double> layers_radii;

        void set_energy(double e);
        double get_energy();

        // Times
        double get_relaxation_time();
        double get_half_mass_relaxation_time();
        double get_crossing_time();

        /** Get the core radius according to a determinated ratio,
         * for example, if it is a Plummer sphere of 1024 particles
         * the core can be visible using a 0.2 percentage */
        double get_core_radius();
        double get_halfmass_radius();
        double3 get_center_of_density();

        void lagrange_radii();
        void core_radius_and_density();
        void get_radii();
        void get_layers();
};

#endif
