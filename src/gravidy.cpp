#include "include/NbodySystem.hpp"

/*
 * Main routing of the initialization and exeucutio of GraviDy
 *
 * The main objects based on the existing classes are created here
 * and passing to some methods to perform the interaction of them
 * to achieve an N-body system integration.
 */
int main(int argc, char *argv[]) {

    // Parsing the command-line parameters
    OptionsParser op(argc,argv);
    if(!op.check_options()) return 1;

    // Initialization of the system
    NbodySystem ns(op);

    // Reading the input file, storing the information in a temporal structure
    ns.read_input_file();

    // Memory allocation on the Host and Device
    ns.alloc_arrays_host();
    #ifdef GPU
    ns.alloc_arrays_device();
    #endif

    // Copy the data from the temporal structure to the Host variables
    ns.copy_initial_data();

    // Creating an Hermite object depending on the compiling rules
    #ifdef GPU
    Hermite4GPU h4(ns.n, ns.e2, ns.eta);
    #else
    Hermite4CPU h4(ns.n, ns.e2, ns.eta);
    #endif

    // Creating the object in charge of printing the program logs
    Logger log(ns.print_log, ns.out_file);
    // Creating the object with some extra functionallity for N-body codes
    NbodyUtils nu(ns.n, ns.h_r, ns.total_mass, RADIUS_RATIO);

    // Calling the integration process using the corresponding Hermite object
    ns.gtime.integration_ini = omp_get_wtime();
    ns.integration(h4, log, nu);
    ns.gtime.integration_end = omp_get_wtime();

    // Freeing memory on the Device and Host
    #ifdef GPU
    ns.free_arrays_device();
    #endif
    ns.free_arrays_host();

    return 0;
}
