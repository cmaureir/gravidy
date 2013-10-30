#include "NbodySystem.hpp"

void NbodySystem::alloc_arrays_host()
{
    int d4_size = n * sizeof(double4);
    int d1_size = n * sizeof(double);
    int i1_size = n * sizeof(int);

    h_r      = new double4[d4_size];
    h_v      = new double4[d4_size];
    h_f      = new Forces[sizeof(Forces) * n];
    h_a2     = new double4[d4_size];
    h_a3     = new double4[d4_size];
    h_old_a  = new double4[d4_size];
    h_old_a1 = new double4[d4_size];
    h_p      = new Predictor[sizeof(Predictor) * n];
    h_ekin   = new double[d1_size];
    h_epot   = new double[d1_size];
    h_t      = new double[d1_size];
    h_dt     = new double[d1_size];
    //h_m      = new float[f1_size];
    h_move   = new int[i1_size];
}

void NbodySystem::free_arrays_host(){
    delete h_r;
    delete h_v;
    delete h_f;
    delete h_a2;
    delete h_a3;
    delete h_old_a;
    delete h_old_a1;
    delete h_p;
    delete h_ekin;
    delete h_epot;
    delete h_t;
    delete h_dt;
    //delete h_m;
    delete h_move;

}
