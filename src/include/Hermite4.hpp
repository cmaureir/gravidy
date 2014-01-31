#ifndef HERMITE4_HPP
#define HERMITE4_HPP
#include "utils/Logger.hpp"
#include "utils/NbodyUtils.hpp"
#include <limits>

class Hermite4 {
    public:
        Hermite4(NbodySystem *ns, Logger *logger, NbodyUtils *nu);
        ~Hermite4();

        NbodySystem *ns;
        Logger      *logger;
        NbodyUtils  *nu;

        int  find_particles_to_move(double ITIME);
        void next_integration_time(double &ATIME);
        void init_dt(double &ATIME);
        void save_old_acc_jrk(int nact);
        void alloc_arrays_host();
        void free_arrays_host();
        void init_data();

        /** Virtual methods to be implemented by the different versions **/
        virtual void integration() {}
        virtual void predicted_pos_vel(double ITIME) {}
        virtual void correction_pos_vel(double ITIME, int nact) {}
        virtual void force_calculation(Predictor pi, Predictor pj, Forces &fi) {}
        virtual void init_acc_jrk() {}
        virtual void update_acc_jrk(int nact) {}

};

#endif // HERMITE4_HPP
