/*
 * Copyright (c) 2016
 *
 * Cristián Maureira-Fredes <cmaureirafredes@gmail.com>
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. The name of the author may not be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifdef GPU
#include "include/gpu/Hermite4GPU.cuh"
#elif _MPI
#include "include/mpi/Hermite4MPI.hpp"
#else // CPU
#include "include/cpu/Hermite4CPU.hpp"
#endif

/*
 * Main routing of the initialization and exeucution of GraviDy
 *
 * The main objects based on the existing classes are created here
 * and passing to some methods to perform the interaction of them
 * to achieve an N-body system integration.
 */

int main(int argc, char *argv[])
{
    #if defined(_MPI)
    int rank, nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    #endif

    // Parsing the command-line parameters
    OptionsParser op(argc,argv);
    if(!op.check_options()) return 1;

    // Initialization of the system
    NbodySystem ns(op);

    // Reading the input file, storing the information in a temporal structure
    #if defined(_MPI)
    if (rank == 0)
    {
        ns.read_input_file();
    }
    ns.alloc_base_attributes(rank);
    #else
    ns.read_input_file();
    ns.alloc_base_attributes(0);
    #endif

    ns.copy_input_data();

    // Creating the object in charge of printing the program logs
    Logger logger(&ns);

    // Creating the object with some extra functionallity for N-body codes
    NbodyUtils nu(&ns);

    // Creating an Hermite object depending on the compiling rules
    #ifdef GPU
    Hermite4GPU h4(&ns, &logger, &nu);
    #elif _MPI
    Hermite4MPI h4(&ns, &logger, &nu, rank, nprocs);
    #else
    Hermite4CPU h4(&ns, &logger, &nu);
    #endif

    // Calling the integration process using the corresponding Hermite object
    h4.integration();

    #if defined(_MPI)
    MPI_Finalize();
    #endif

    return 0;
}
