#include "NbodySystem.hpp"

NbodySystem::NbodySystem(OptionsParser op)
{
    input_filename   = op.input_filename;
    output_filename  = op.output_filename;
    integration_time = op.integration_time;
    e2               = op.softening * op.softening;
    eta              = op.eta;
    print_log        = op.print_log;
    iterations       = 0;
    en.ini = 0.0;
    en.end = 0.0;
    en.tmp = 0.0;
    gtime.gflops = 0.0;
    if (print_log)
    {
        std::cout << "abrimos" << std::endl;
        out_file.open(output_filename.c_str(), std::ios::out);
    }
}

NbodySystem::~NbodySystem()
{
    if (print_log)
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

void NbodySystem::copy_initial_data()
{

    double4 empty = {0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < (int)reader.size(); i++) {

        h_r[i].w    = reader[i].m;
        h_p[i].m    = reader[i].m;

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

        h_old[i].a[0]   = 0.0;
        h_old[i].a[1]   = 0.0;
        h_old[i].a[2]   = 0.0;

        h_old[i].a1[0]  = 0.0;
        h_old[i].a1[1]  = 0.0;
        h_old[i].a1[2]  = 0.0;

        h_a2[i]      = empty;
        h_a3[i]      = empty;

        h_t[i]       = 0.0;
        h_dt[i]      = 0.0;

        h_move[i]    = 0;
    }
}

double NbodySystem::get_energy()
{
    en.potential = 0.0;
    en.kinetic   = 0.0;

    #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        double epot_tmp = 0.0;
        for (int j = i+1; j < n; j++)
        {
            double rx = h_r[j].x - h_r[i].x;
            double ry = h_r[j].y - h_r[i].y;
            double rz = h_r[j].z - h_r[i].z;
            double r2 = rx*rx + ry*ry + rz*rz + e2;

            epot_tmp -= (h_r[i].w * h_r[j].w) / sqrt(r2);
        }

        double vx = h_v[i].x * h_v[i].x;
        double vy = h_v[i].y * h_v[i].y;
        double vz = h_v[i].z * h_v[i].z;
        double v2 = vx + vy + vz;

        double ekin_tmp = 0.5 * h_r[i].w * v2;

        #pragma omp atomic
        en.kinetic += ekin_tmp;
        #pragma omp atomic
        en.potential += epot_tmp;
    }
    return en.kinetic + en.potential;
}
