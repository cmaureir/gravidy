#ifndef HERMITE4CPU_HPP
#define HERMITE4CPU_HPP
#include "Hermite4.hpp"

class Hermite4CPU : public Hermite4 {
    public:
        Hermite4CPU(int n, double e2, float eta) : Hermite4(n, e2, eta) { }

        void force_calculation(int i, int j, Predictor *p, Forces *f);
        void init_acc_jrk(Predictor *p, Forces* f);
        void update_acc_jrk(int nact, int *move, Predictor *p, Forces* f, Gtime &gtime);
};

#endif
