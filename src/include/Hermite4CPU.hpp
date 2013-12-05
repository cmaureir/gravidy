#ifndef HERMITE4CPU_HPP
#define HERMITE4CPU_HPP
#include "Hermite4.hpp"

class Hermite4CPU : public Hermite4 {
    public:
        Hermite4CPU(NbodySystem *ns, Logger *logger, NbodyUtils *nu)
            : Hermite4(ns, logger, nu) {}
        ~Hermite4CPU();

        void force_calculation(int i, int j);
        void init_acc_jrk();
        void update_acc_jrk(int nact);
        void predicted_pos_vel(double ITIME);
        void correction_pos_vel(double ITIME, int nact);
        void integration();
};

#endif
