#ifndef NBODYUTILS_HPP
#define NBODYUTILS_HPP
#include "../common.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

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
    public:
        NbodyUtils(int n, double4 *r, float m);
        ~NbodyUtils();

        int n;
        float total_mass;
        double energy;
        double3 cod;
        double4 *r;

        void set_energy(double e);
        double get_core_radius(double3 pc);
        double get_half_mass_relaxation_time();
        double get_crossing_time();
        double get_halfmass_radius(double3 dc);
        double3 get_center_of_density();
};

#endif
