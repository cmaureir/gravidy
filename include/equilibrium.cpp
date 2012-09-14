#include "equilibrium.hpp"

point get_center_of_density()
{
    std::vector<double> p;
    std::vector<distance> d;
    distance dist;
    double rx, ry, rz;
    float radius, dsum;
    float aa, bb;

    for (int i = 0; i < n; i++)
    {
        d.clear();
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                rx = (h_r[j].x - h_r[i].x);
                ry = (h_r[j].y - h_r[i].y);
                rz = (h_r[j].z - h_r[i].z);
                dist.value = sqrt(rx*rx + ry*ry + rz*rz);
                dist.index = j;
                d.push_back(dist);
            }
        }
        std::sort(d.begin(), d.end());
        radius = d[J-1].value;
        aa = (J-1) * h_m[i];
        bb = (4.0 * PI * radius * radius * radius)/3.0;
        p.push_back( aa / bb );
    }

    point density_center;
    density_center.x = 0.0f;
    density_center.y = 0.0f;
    density_center.z = 0.0f;

    dsum = 0.0f;

    for (int i = 0; i < n; i++) {
        dsum += p[i];
        density_center.x += h_r[i].x * p[i];
        density_center.y += h_r[i].y * p[i];
        density_center.z += h_r[i].z * p[i];
    }

    density_center.x /= dsum;
    density_center.y /= dsum;
    density_center.z /= dsum;

    return density_center;
}

double get_halfmass_radius(double dx, double dy, double dz)
{
    float half_mass;
    double tmp, r_h;
    int i, j;

    distance d_tmp;
    std::vector<distance> distances;

    half_mass = 0;
    j = 0;

    for (i = 0; i < n; i++)
    {
        tmp = sqrt( (dx - h_r[i].x) * (dx - h_r[i].x) +
                    (dy - h_r[i].y) * (dy - h_r[i].y) +
                    (dz - h_r[i].z) * (dz - h_r[i].z) );
        d_tmp.index = i;
        d_tmp.value  = tmp;
        distances.push_back(d_tmp);
    }

    std::sort(distances.begin(), distances.end());

    for (i = 0; i < n; i++) {
        if(half_mass == total_mass/2)
        {
            j = i;
            break;
        }
        half_mass += h_m[distances[i].index];
    }

    r_h = distances[j-1].value;

    return r_h;
}


float get_crossing_time()
{
    float M = total_mass;

    float Rv = (-G * M * M) / (4 * (energy_ini));
    float Ut = sqrt( (Rv * Rv * Rv) / G * M);
    float t_cr = 2 * sqrt(2) * Ut;

    std::cout << "Virial radius: " << Rv << std::endl;
    std::cout << "Unity of time: " << Ut << std::endl;

    return t_cr;
}

float get_relaxation_time()
{
    float t_rh;
    float r_h, a, b;

    point center = get_center_of_density();
    r_h = get_halfmass_radius(center.x, center.y, center.z);
    a = sqrt( (n * r_h * r_h * r_h) / ( G * h_m[0]) );
    b = 1/log(0.11 * n);

    t_rh = 0.138 * a * b;
    return t_rh;
}
