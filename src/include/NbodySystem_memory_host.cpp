#include "NbodySystem.hpp"

void NbodySystem::alloc_arrays_host()
{
    int d4_size = n * sizeof(double4);
    int d1_size = n * sizeof(double);
    int i1_size = n * sizeof(int);
    int ff_size = n * sizeof(Forces);
    int pp_size = n * sizeof(Predictor);

    h_r      = new double4[d4_size];
    h_v      = new double4[d4_size];
    h_f      = new Forces[ff_size];
    h_p      = new Predictor[pp_size];
    h_i      = new Predictor[pp_size];
    h_ekin   = new double[d1_size];
    h_epot   = new double[d1_size];
    h_t      = new double[d1_size];
    h_dt     = new double[d1_size];
    h_old    = new Forces[ff_size];
    h_fout_tmp= new Forces[ff_size*NJBLOCK];
    h_a2     = new double4[d4_size];
    h_a3     = new double4[d4_size];
    //h_m      = new float[f1_size];
    h_move   = new int[i1_size];
}

void NbodySystem::free_arrays_host(){
    delete h_r;
    delete h_v;
    delete h_f;
    delete h_a2;
    delete h_a3;
    delete h_old;
    delete h_p;
    delete h_i;
    delete h_ekin;
    delete h_epot;
    delete h_t;
    delete h_dt;
    delete h_move;
    //delete h_m;
    delete h_fout_tmp;

}
