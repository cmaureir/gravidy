#include "include/NbodySystem.hpp"

int main(int argc, char *argv[]) {

    NbodySystem ns;
    OptionsParser op(argc,argv);

    if(!op.check_options()) return 1;

    ns.get_parameters(op);

    ns.read_input_file();

    ns.alloc_arrays_host();
    ns.copy_initial_data();

    Hermite4 h(ns.n, ns.e2, ns.eta);
    ns.gtime.integration_ini = omp_get_wtime();
    ns.integration(h);
    ns.gtime.integration_end = omp_get_wtime();

    ns.free_arrays_host();

    return 0;
}
