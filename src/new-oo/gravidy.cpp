#include "include/NbodySystem.hpp"

int main(int argc, char *argv[]) {

    // Parsing the command-line parameters
    OptionsParser op(argc,argv);
    if(!op.check_options()) return 1;

    // Initialization of the system
    NbodySystem ns(op);

    ns.read_input_file();
    ns.alloc_arrays_host();
    #ifdef GPU
    ns.alloc_arrays_device();
    #endif
    ns.copy_initial_data();

    #ifdef GPU
    Hermite4GPU h4(ns.n, ns.e2, ns.eta);
    #else
    Hermite4CPU h4(ns.n, ns.e2, ns.eta);
    #endif

    Logger log(ns.print_log, ns.out_file);
    NbodyUtils nu(ns.n, ns.h_r, ns.total_mass);

    ns.gtime.integration_ini = omp_get_wtime();
    ns.integration(h4, log, nu);
    ns.gtime.integration_end = omp_get_wtime();

    #ifdef GPU
    ns.free_arrays_device();
    #endif
    ns.free_arrays_host();

    return 0;
}
