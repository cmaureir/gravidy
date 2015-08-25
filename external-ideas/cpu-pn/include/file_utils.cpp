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
