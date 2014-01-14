#include "NbodySystem.hpp"

/** Constructor */
NbodySystem::NbodySystem(OptionsParser op)
{
    input_filename   = op.input_filename;
    output_filename  = op.output_filename;
    integration_time = op.integration_time;
    e2               = op.softening * op.softening;
    eta              = op.eta;
    ops              = op.ops;

    iterations       = 0;
    en.ini = 0.0;
    en.end = 0.0;
    en.tmp = 0.0;
    gtime.gflops = 0.0;
    if (! ops.print_screen)
    {
        out_file.open(output_filename.c_str(), std::ios::out);
    }
}

/** Destructor */
NbodySystem::~NbodySystem()
{
    if (! ops.print_screen)
    {
        out_file.close();
    }
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
                tmp.m    = (float)strtod(tokens[0].c_str(), NULL);
                tmp.r[0] = strtod(tokens[1].c_str(), NULL);
                tmp.r[1] = strtod(tokens[2].c_str(), NULL);
                tmp.r[2] = strtod(tokens[3].c_str(), NULL);
                tmp.v[0] = strtod(tokens[4].c_str(), NULL);
                tmp.v[1] = strtod(tokens[5].c_str(), NULL);
                tmp.v[2] = strtod(tokens[6].c_str(), NULL);
            }
            reader.push_back(tmp);
            total_mass += tmp.m;
            getline(file, line);
        }
        file.close();
    }
    n = (int)reader.size();
}

void NbodySystem::alloc_base_attributes(int rank)
{

    #ifdef MPI
    // Sending the amount of particles
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Resizing input data vector, to be able to Bcast it
    if (rank > 0)
        reader.resize(n);

    MPI_Bcast(&reader[0], sizeof(file_data) * n, MPI_BYTE, 0, MPI_COMM_WORLD);
    #endif

    int d4_size = n * sizeof(double4);
    h_r      = new double4[d4_size];
    h_v      = new double4[d4_size];

}

void NbodySystem::free_base_attributes()
{
    int d4_size = n * sizeof(double4);

    h_r      = new double4[d4_size];
    h_v      = new double4[d4_size];

}
void NbodySystem::copy_input_data()
{
    for (int i = 0; i < (int)reader.size(); i++)
    {
        h_r[i].w    = reader[i].m;
        h_r[i].x    = reader[i].r[0];
        h_r[i].y    = reader[i].r[1];
        h_r[i].z    = reader[i].r[2];
        h_v[i].x    = reader[i].v[0];
        h_v[i].y    = reader[i].v[1];
        h_v[i].z    = reader[i].v[2];
    }
}


