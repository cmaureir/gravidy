#include "file_utils.hpp"

/*
 * @fn read_input_file(std::string path)
 *
 * @param path Input file route
 *
 * @brief
 *  Function that read an input file
 *  with the following structure:
 *
 *  mass rx ry rz vx vy vz
 *
 */
void read_input_file(std::string path)
{
    std::ifstream file(path.c_str());
    particle part_tmp;
    std::istream_iterator<std::string> word_iter(file), word_iter_end;

    total_mass = 0;

    for( /* empty */ ; word_iter != word_iter_end; ++word_iter)
    {
        part_tmp.m   = strtod((*word_iter).c_str(), NULL); word_iter++;
        part_tmp.r.x = strtod((*word_iter).c_str(), NULL); word_iter++;
        part_tmp.r.y = strtod((*word_iter).c_str(), NULL); word_iter++;
        part_tmp.r.z = strtod((*word_iter).c_str(), NULL); word_iter++;
        part_tmp.v.x = strtod((*word_iter).c_str(), NULL); word_iter++;
        part_tmp.v.y = strtod((*word_iter).c_str(), NULL); word_iter++;
        part_tmp.v.z = strtod((*word_iter).c_str(), NULL);
        part.push_back(part_tmp);
        total_mass += part_tmp.m;
    }

    file.close();
    n = (int)part.size();
}


/*
 * @fn write_output_file(std::string path)
 *
 * @param path Output file route
 *
 * @brief
 *  Function that writes an output file
 *  with the following structure:
 *
 *  rx ry rz vx vy vz ax ay az
 *
 */
void write_output_file(std::string path)
{
    if(path == "")
    {
        path = std::string("output/");
        path.append(input_file.begin()+6, input_file.end());
        path.append("_");
        std::ostringstream strs;
        strs << int_time;
        path.append(strs.str());
        path.append("_");
        path.append(run);
        path.append(".out");
    }

    std::ofstream file(path.c_str());
    std::stringstream sline;
    //sline.precision(std::numeric_limits<long float>::digits10);

    sline << "Run mode: " << run;
    file << sline.str(); sline.str("");
    sline << "\nInit Acc and Jerk calculation time: "  << init_time;
    file << sline.str(); sline.str("");
    sline << "\nEnergy Init: "  << energy_ini;
    file << sline.str(); sline.str("");
    sline << "\nEnergy End:  "  << energy_end;
    file << sline.str(); sline.str("");
    sline << "\nIntegration time: "  << end_time - ini_time;
    file << sline.str(); sline.str("");
    sline << "\nIterations number: "  << iterations;
    file << sline.str(); sline.str("");
    file << "\n";

//    for (int i = 0; i < n; i++) {
//        sline  << setw(10)  << i;
//        file << sline.str(); sline.str("");
//        sline  << setw(22) << h_r[i].x;
//        file << sline.str(); sline.str("");
//        sline  << setw(22) << h_r[i].y;
//        file << sline.str(); sline.str("");
//        sline  << setw(22) << h_r[i].z;
//        file << sline.str(); sline.str("");
//        file << "\n";
//    }

    file.close();
}
