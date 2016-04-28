#include "NbodySystem.hpp"

/** Constructor */
NbodySystem::NbodySystem(OptionsParser op)
{
    input_filename   = op.input_filename;
    output_filename  = op.output_filename;
    resume_filename  = op.resume_filename;
    snapshot_filename  = op.snapshot_filename;
    integration_time = op.integration_time;
    snapshot_time    = op.snapshot_time;
    snapshot_number  = op.snapshot_number;
    e2               = op.softening * op.softening;
    eta              = op.eta;
    ops              = op.ops;
    iterations       = 0;
    gpus             = op.gpus;
    resume           = op.resume;

    en.ini = 0.0;
    en.end = 0.0;
    en.tmp = 0.0;
    gtime.gflops = 0.0;

    r_cl = 0.0;
    dt_cl = 0.0;
    m_g = 0.0;
    t_rlx = 0.0;
    t_cr = 0.0;
    r_virial = 0.0;
}

/** Destructor */
NbodySystem::~NbodySystem()
{
    /* Empty constructor */
}

void NbodySystem::read_input_file()
{
    file_data tmp;
    total_mass = 0;
    max_mass = 0.0;
    bool got_time = false;
    std::string line;

    std::ifstream file;
    if (resume)
    {
        file.open(snapshot_filename.c_str());
    }
    else
    {
        file.open(input_filename.c_str());
    }

    if (file.is_open())
    {
        getline(file, line);
        while (file.good())
        {

            if(file.eof())
                break;

            // Lines that start with a "#"

            int first_no_space = line.find_first_not_of(" ");
            if (line[first_no_space] == std::string("#")[0])
            {
                if (!got_time)
                {
                    // TO DO: look for integration time if restart option

                    // remove all the "#" and ":"
                    char chars[] = "#:";
                    for (unsigned int i = 0; i < strlen(chars); i++)
                    {
                        line.erase(std::remove(line.begin(),
                                               line.end(),
                                               chars[i]),
                                               line.end());
                    }
                    // Lowering the case of the line
                    std::transform(line.begin(), line.end(), line.begin(), ::tolower);

                    // Split string
                    std::vector<std::string> ctokens;
                    std::istringstream ciss(line);
                    copy(std::istream_iterator<std::string>(ciss),
                         std::istream_iterator<std::string>(),
                         std::back_inserter<std::vector<std::string> >(ctokens));

                    for (int t = 0; t < (int)ctokens.size(); t++)
                    {
                        if (ctokens[t] == "time" && (t+1 < (int)ctokens.size()))
                        {
                            snapshot_time = strtod(ctokens[t+1].c_str(), NULL);
                            got_time = true;
                            std::cout << "It's a snapshot" << std::endl;
                    //        getline(file, line);
                            break;
                        }
                    }

                    getline(file, line);
                    continue;
                }
                else // commented line
                {
                    //getline(file, line);
                }
            }
            else
            {
                std::vector<std::string> tokens;
                std::istringstream iss(line);
                copy(std::istream_iterator<std::string>(iss),
                     std::istream_iterator<std::string>(),
                     std::back_inserter<std::vector<std::string> >(tokens));

                // We need only 7 parameters id, m, rx, ry, rz, vx, vy, vz
                if (tokens.size() != 8)
                {
                    std::cerr << "gravidy: wrong line format in"
                              <<input_filename
                              << ": only eight numbers per line are allowed"
                              << std::endl;
                    break;
                }
                else
                {
                    tmp.id   = (int)std::atoi(tokens[0].c_str());
                    tmp.m    = (float)strtod(tokens[1].c_str(), NULL);
                    tmp.r[0] = strtod(tokens[2].c_str(), NULL);
                    tmp.r[1] = strtod(tokens[3].c_str(), NULL);
                    tmp.r[2] = strtod(tokens[4].c_str(), NULL);
                    tmp.v[0] = strtod(tokens[5].c_str(), NULL);
                    tmp.v[1] = strtod(tokens[6].c_str(), NULL);
                    tmp.v[2] = strtod(tokens[7].c_str(), NULL);
                    if (tmp.m > max_mass)
                        max_mass = tmp.m;

                }
                reader.push_back(tmp);
                total_mass += tmp.m;
            }
            getline(file, line);
        }
        file.close();
    }
    n = (int)reader.size();
}

void NbodySystem::alloc_base_attributes(int rank)
{

    #if defined(_MPI) || defined(MPIGPU)
    // Sending the amount of particles
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Resizing input data vector, to be able to Bcast it
    if (rank > 0)
        reader.resize(n);

    MPI_Bcast(&reader[0], sizeof(file_data) * n, MPI_BYTE, 0, MPI_COMM_WORLD);
    #endif

    h_id = new int[n];
    h_r = new double4[n];
    h_v = new double4[n];
}

void NbodySystem::free_base_attributes()
{
    delete h_id;
    delete h_r;
    delete h_v;
}

void NbodySystem::copy_input_data()
{
    for (int i = 0; i < (int)reader.size(); i++)
    {
        h_id[i]     = reader[i].id;
        h_r[i].w    = reader[i].m;
        h_r[i].x    = reader[i].r[0];
        h_r[i].y    = reader[i].r[1];
        h_r[i].z    = reader[i].r[2];
        h_v[i].x    = reader[i].v[0];
        h_v[i].y    = reader[i].v[1];
        h_v[i].z    = reader[i].v[2];
    }
}


