#ifdef GPU
#include "include/gpu/Hermite4GPU.cuh"
#elif KEPLER
#include "include/kepler/Hermite4Kepler.hpp"
#else
#include "include/cpu/Hermite4CPU.hpp"
#endif

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
    ns.alloc_base_attributes();
    ns.copy_input_data();

    // Creating the object in charge of printing the program logs
    Logger logger(&ns);
    // Creating the object with some extra functionallity for N-body codes
    NbodyUtils nu(&ns,RADIUS_RATIO);

    // Creating an Hermite object depending on the compiling rules
    #ifdef GPU
    Hermite4GPU h4(&ns, &logger, &nu);
    #elif KEPLER
    Hermite4Kepler h4(&ns, &logger, &nu);
    #else
    Hermite4CPU h4(&ns, &logger, &nu);
    #endif

    // Calling the integration process using the corresponding Hermite object
    h4.integration();

    return 0;
}
