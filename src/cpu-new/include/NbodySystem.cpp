#include "NbodySystem.hpp"

NbodySystem::NbodySystem()
{
    std::cout << "NbodySystem::Creando" << std::endl;
}

NbodySystem::~NbodySystem()
{
    std::cout << "NbodySystem::Destruyendo" << std::endl;
}

void NbodySystem::get_parameters(OptionsParser op)
{
    input_filename   = op.input_filename;
    output_filename  = op.output_filename;
    integration_time = op.integration_time;
    e2               = op.softening * op.softening;
    eta              = op.eta;
    print_log        = op.print_log;

}

void NbodySystem::read_input_file()
{
    file_data tmp;
    total_mass = 0;
    std::string line;

    std::ifstream file(input_filename.c_str());

    if (file.is_open())
    {
        getline(file, line);
        while (file.good())
        {

            if(file.eof())
                break;
            if (line[0] == std::string("#")[0])
                continue;
            std::vector<std::string> tokens;
            std::istringstream iss(line);
            copy(std::istream_iterator<std::string>(iss),
                 std::istream_iterator<std::string>(),
                 std::back_inserter<std::vector<std::string> >(tokens));

            // We need only 7 parameters m, rx, ry, rz, vx, vy, vz
            if (tokens.size() != 7)
            {
                std::cerr << "gravidy: wrong line format in"
                          <<input_filename
                          << ": only seven numbers per line are allowed"
                          << std::endl;
                break;
            }
            else
            {
                tmp.m    = atof(tokens[0].c_str());
                tmp.r[0] = atof(tokens[1].c_str());
                tmp.r[1] = atof(tokens[2].c_str());
                tmp.r[2] = atof(tokens[3].c_str());
                tmp.v[0] = atof(tokens[4].c_str());
                tmp.v[1] = atof(tokens[5].c_str());
                tmp.v[2] = atof(tokens[6].c_str());
            }
            reader.push_back(tmp);
            total_mass += tmp.m;
            getline(file, line);
        }
        file.close();
    }
    n = (int)reader.size();
}

void NbodySystem::copy_initial_data()
{

    float4 empty = {0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < (int)reader.size(); i++) {

        h_r[i].w    = reader[i].m;

        h_r[i].x    = reader[i].r[0];
        h_r[i].y    = reader[i].r[1];
        h_r[i].z    = reader[i].r[2];

        h_v[i].x    = reader[i].v[0];
        h_v[i].y    = reader[i].v[1];
        h_v[i].z    = reader[i].v[2];

        h_p[i].r[0] = reader[i].r[0];
        h_p[i].r[1] = reader[i].r[1];
        h_p[i].r[2] = reader[i].r[2];

        h_p[i].v[0] = reader[i].v[0];
        h_p[i].v[1] = reader[i].v[1];
        h_p[i].v[2] = reader[i].v[2];

        h_f[i].a[0]  = 0.0;
        h_f[i].a[1]  = 0.0;
        h_f[i].a[2]  = 0.0;

        h_f[i].a1[0] = 0.0;
        h_f[i].a1[1] = 0.0;
        h_f[i].a1[2] = 0.0;

        h_a2[i]      = empty;
        h_a3[i]      = empty;

        h_old_a[i]   = empty;
        h_old_a1[i]  = empty;

        h_t[i]       = 0.0;
        h_dt[i]      = 0.0;

        h_move[i]    = 0;
    }
}

void NbodySystem::alloc_arrays_host()
{
    int f4_size = n * sizeof(float4);
    int d4_size = n * sizeof(double4);
    int d1_size = n * sizeof(double);
    int f1_size = n * sizeof(float);
    int i1_size = n * sizeof(int);

    h_r      = new double4[d4_size];
    h_v      = new double4[d4_size];
    h_f      = new Forces[sizeof(Forces) * n];
    h_a2     = new float4[f4_size];
    h_a3     = new float4[f4_size];
    h_old_a  = new float4[f4_size];
    h_old_a1 = new float4[f4_size];
    h_p      = new Predictor[sizeof(Predictor) * n];
    h_ekin   = new double[d1_size];
    h_epot   = new double[d1_size];
    h_t      = new double[d1_size];
    h_dt     = new double[d1_size];
    h_m      = new float[f1_size];
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
    delete h_m;
    delete h_move;

}

void NbodySystem::integration(Hermite4 h)
{

    double ATIME = 1.0e+10; // Actual integration time
    double ITIME = 0.0;     // Integration time
    int nact     = 0;       // Active particles
    int nsteps   = 0;       // Amount of steps per particles on the system
    static long long interactions = 0;

    //int max_threads = omp_get_max_threads();
    //omp_set_num_threads( max_threads - 1);

    h.init_acc_jrk(h_r, h_v, h_f);     // Initial calculation of a and a1
    h.init_dt(ATIME, h_f, h_dt, h_t);  // Initial calculation of time-steps using simple equation

    energy_ini = get_energy();   // Initial calculation of the energy of the system

    //print_log(ITIME, iterations, interactions, nsteps, energy_ini, energy_ini);

    while (ITIME < integration_time)
    {
        ITIME = ATIME;                         // New integration time
        nact = h.find_particles_to_move(ITIME, h_dt, h_t, h_move);  // Find particles to move (nact)

        //save_old(nact);                        // Save old information

        //h.predicted_pos_vel(ITIME);              // Predict all the particles
        //h.update_acc_jrk(nact);                  // Update a and a1 of nact particles
        //h.correction_pos_vel(ITIME, nact);       // Correct r and v of nact particles


        // Update the amount of interactions counter
        interactions += nact * n;

        // Find the next integration time
        //next_itime(&ATIME);

        // Print log every integer ITIME
        //if(std::ceil(ITIME) == ITIME)
        if(nact == n)          // Print log in every integer ITIME
        {
           //print_log(ITIME, iterations, interacions, nsteps, energy_ini, get_energy());
        }

        // Update nsteps with nact
        nsteps += nact;

        // Increase iteration counter
        iterations++;
    }
}

double NbodySystem::get_energy()
{
    total_epot = 0.0;
    total_ekin = 0.0;

    //#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        double epot_tmp = 0.0;
        for (int j = i+1; j < n; j++)
        {
            double rx = h_r[j].x - h_r[i].x;
            double ry = h_r[j].y - h_r[i].y;
            double rz = h_r[j].z - h_r[i].z;
            double r2 = rx*rx + ry*ry + rz*rz;

            epot_tmp -= (h_m[i] * h_m[j]) / sqrt(r2);
        }

        double vx = h_v[i].x * h_v[i].x;
        double vy = h_v[i].y * h_v[i].y;
        double vz = h_v[i].z * h_v[i].z;
        double v2 = vx + vy + vz;

        double ekin_tmp = 0.5 * h_m[i] * v2;

        //#pragma omp atomic
        total_ekin += ekin_tmp;
        //#pragma omp atomic
        total_epot += epot_tmp;
    }
    return total_epot + total_ekin;
}
