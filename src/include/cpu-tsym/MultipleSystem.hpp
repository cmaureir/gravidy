#ifndef MULTIPLESYSTEM_HPP
#define MULTIPLESYSTEM_HPP
#include "../utils/NbodyUtils.hpp"
#include <vector>

const short MAX_MEMBERS = 2;

typedef struct ShortParticle
{
    double4 r;
    double4 v;
    Forces f;
} SParticle;

typedef struct MultipleSystemParticle
{
    int id;
    double4 r;
    double4 v;
    Forces f;
    Forces old_f;
    Predictor p;
    Predictor p0;
    double3 a2;
    double3 a3;
    double t;
    double dt;
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

        /// Accuracy parameter for timestep calculation
        double ETA_B;

        /// Amount of member in the system
        short members;

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

        /// @brief Forming
        MParticle get_new_particle(MParticle p1, MParticle p2);

        void print();
        void init_timestep(double &CTIME);
        void next_itime(double &CTIME);
        //void update_timestep(double ITIME);
        void update_information(double ITIME);
        void prediction(double ITIME);
        void force_calculation(Predictor pi, Predictor pj, Forces &fi);
        void evaluation();
        void correction(double ITIME);
        void save_old();
        double get_timestep_normal(MParticle p);
        double get_energy();
};
#endif
