#ifndef MULTIPLESYSTEM_HPP
#define MULTIPLESYSTEM_HPP
#include "../utils/NbodyUtils.hpp"
#include <vector>

inline double4 operator+(const double4 &a, const double4 &b)
{
    double4 tmp = {a.x + b.x, a.y + b.y, a.z + b.z,  a.w + b.w};
    return tmp;
}

inline void operator+=(double4 &a, double4 &b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    a.w += b.w;
}

inline Forces operator+(Forces a, Forces b)
{
    Forces tmp;
    tmp.a[0] = a.a[0] + b.a[0];
    tmp.a[1] = a.a[1] + b.a[1];
    tmp.a[2] = a.a[2] + b.a[2];

    tmp.a1[0] = a.a1[0] + b.a1[0];
    tmp.a1[1] = a.a1[1] + b.a1[1];
    tmp.a1[2] = a.a1[2] + b.a1[2];

    return tmp;
}

inline Forces operator-(const Forces &a, const Forces &b)
{
    Forces tmp;
    tmp.a[0] = a.a[0] - b.a[0];
    tmp.a[1] = a.a[1] - b.a[1];
    tmp.a[2] = a.a[2] - b.a[2];

    tmp.a1[0] = a.a1[0] - b.a1[0];
    tmp.a1[1] = a.a1[1] - b.a1[1];
    tmp.a1[2] = a.a1[2] - b.a1[2];

    return tmp;
}

inline Forces operator*(Forces a, Forces b)
{
    Forces tmp;
    tmp.a[0] = a.a[0] * b.a[0];
    tmp.a[1] = a.a[1] * b.a[1];
    tmp.a[2] = a.a[2] * b.a[2];

    tmp.a1[0] = a.a1[0] * b.a1[0];
    tmp.a1[1] = a.a1[1] * b.a1[1];
    tmp.a1[2] = a.a1[2] * b.a1[2];

    return tmp;
}

inline void operator+=(Forces &a, Forces &b)
{
    a.a[0] += b.a[0];
    a.a[1] += b.a[1];
    a.a[2] += b.a[2];

    a.a1[0] += b.a1[0];
    a.a1[1] += b.a1[1];
    a.a1[2] += b.a1[2];
}

const short MAX_MEMBERS = 2;

typedef struct ShortParticle
{
    double4 r;
    double4 v;
    Forces f;
    Forces old;
} SParticle;

typedef struct MultipleSystemParticle
{
    int id;
    double4 r;
    double4 v;
    Forces f;
    Forces old;
    Forces pert;
    Predictor p;
    Predictor p0;
    double3 a2;
    double3 a3;
    double t;
    double dt;
    double old_dt;
    bool status;
} MParticle;

class MultipleSystem
{
    public:
        /// @brief Constructor of the class
        /// @param ns pointer to the NbodySystem object which contain
        //            all the information of the system
        MultipleSystem(NbodySystem *ns, NbodyUtils *nu);

        /// @brief Destructor of the class.
        ~MultipleSystem();

        /// Local reference to the NbodySystem pointer
        NbodySystem *ns;
        NbodyUtils  *nu;

        /// Center of Mass information
        SParticle com;
        SParticle com_evolving;

        /// Accuracy parameter for timestep calculation
        double ETA_B;

        /// Subsystem distance
        double dist;
        int *pert;
        int num_pert;
        double total_mass;

        /// Amount of member in the system
        short members;

        double ini_e;
        double bin_e;
        double semimajor_ini;
        double semimajor_end;
        double ecc_ini;
        double ecc_end;

        /// Particle members of the multiple system
        MParticle parts[MAX_MEMBERS];

        /// @brief Method for including a new particle member to the system
        /// @todo Not working yet
        void add_particle(int id);

        /// @brief Method for getting the center of mass between two MParticle
        ///        structures.
        /// @param p1 Mparticle member
        /// @param p2 Mparticle member
        /// @return SParticle with the center of mass information.
        //          (Also called virtual or ghost particle)
        SParticle get_center_of_mass(MParticle p1, MParticle p2);

        void print();
        void init_timestep();
        void next_itime(double &CTIME);
        void update_timestep(double ITIME);
        void update_information(double ITIME);
        void prediction(double ITIME);
        void force_calculation(Predictor pi, Predictor pj, Forces &fi);
        void evaluation(int *nb_list);
        void get_perturbers(int com_id);
        void predict_com(double CTIME, int com_id);
        void predict_perturbers(double CTIME);
        void perturbers(double CTIME, int ite);
        void correction(double ITIME, bool check);
        void save_old();
        double get_timestep_normal(MParticle p);
        void get_orbital_elements(bool t0);
        double get_energy();
        void adjust_particles(SParticle sp);
        void adjust_to_center_of_mass();
};
#endif
