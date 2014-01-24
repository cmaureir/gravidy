#ifndef HERMITE4KEPLER_HPP
#define HERMITE4KEPLER_HPP
#include "../Hermite4.hpp"

#include <cmath>
#include <iostream>

static double
    _1_17_16 = 1./17./16.,
    _1_16_15 = 1./16./15.,
    _1_15_14 = 1./15./14.,
    _1_14_13 = 1./14./13.,
    _1_13_12 = 1./13./12.,
    _1_12_11 = 1./12./11.,
    _1_11_10 = 1./11./10.,
    _1_10_9  = 1./10./9.,
    _1_9_8   = 1./9./8.,
    _1_8_7   = 1./8./7.,
    _1_7_6   = 1./7./6.,
    _1_6_5   = 1./6./5.,
    _1_4_3   = 1./4./3.,
    _1_3_2   = 1./3./2.;

typedef struct orbital_elements
{
    double3 j; // Angular momentum vector
    double3 e; // Runge-Lenz vector
    double a;  // Semi-major axis
    double b;  // Semi-minor axis
    double w;  // Frequency of the orbit
    double3 a_vec; // Semi-major axis vector
    double3 b_vec; // Semi-minor axis vector
    double ecc;       // Eccentricity
    double e_anomaly; // Eccentric anomaly
    double e_anomaly0; // initial Eccentric anomaly
    double m_anomaly; // Mean anomaly
    double m_anomaly0; // initial Mean anomaly
} orbital_elements;

class Hermite4Kepler : public Hermite4 {
    public:
        Predictor pkepler;
        orbital_elements oe;

        double r_mag; // Relative position magnitude
        double rdotv; // r dot v product
        double v_mag; // Relative velocity magnitude

        Hermite4Kepler(NbodySystem *ns, Logger *logger, NbodyUtils *nu);
        ~Hermite4Kepler();

        void alloc_arrays_host_kepler();
        void free_arrays_host_kepler();
        void init_data_bh();

        // Hermite
        void force_calculation(int i, int j);
        void init_acc_jrk();
        void update_acc_jrk(int nact);
        void predicted_pos_vel(double ITIME);
        void correction_pos_vel(double ITIME, int nact);
        void integration();

        // Kepler
        void force_calculation_bh(int i);
        void init_acc_jrk_bh();

        void calculate_orbital_elements();
        void print_orbital_elements();


        void move_keplerian_orbit(double ITIME, int nact);
        void kepler_move(int i, double dt);

        double solve_kepler(double m_anomaly, double ecc);
        double kepler(const double ecc, double mean_anom);

        Predictor get_elliptical_pos_vel(double dt);
        Predictor get_hyperbolic_pos_vel(double dt);

};
#endif // HERMITE4KEPLER_HPP
