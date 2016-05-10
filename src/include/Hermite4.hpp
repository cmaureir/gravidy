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

        unsigned int  find_particles_to_move(double ITIME);
        void next_integration_time(double &ATIME);
        void init_dt(double &ATIME, float ETA, double ITIME);
        void save_old_acc_jrk(unsigned int nact);
        void alloc_arrays_host();
        void free_arrays_host();
        void init_data();

        /** Virtual methods to be implemented by the different versions **/
        virtual void integration() {}
        virtual void predicted_pos_vel(double itime, double *t, double4 *r,
                                       double4 *v, Forces *f, Predictor *p) {}
        virtual void correction_pos_vel(double itime, unsigned int nact, double *dt,
                                        double *t, unsigned int *move, Predictor *p,
                                        Forces *f, Forces *old, double3 *a2,
                                        double3 *a3, double4 *r, double4 *v) {}
        virtual void force_calculation(Predictor pi, Predictor pj, Forces &fi) {}
        virtual void init_acc_jrk(Predictor *p, Forces *f, double *r_sphere) {}
        virtual void update_acc_jrk(unsigned int nact, unsigned int *move, double *r_sphere,
                                    Predictor *p, Forces *f) {}

};

#endif // HERMITE4_HPP
