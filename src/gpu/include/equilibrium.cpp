#include "equilibrium.hpp"


/*
 * @fn void get_system_info()
 *
 * @brief Print some general information of the N-body system
 *
 */
void get_system_info()
{
    Point p_c = get_center_of_density();
    t_rh = get_relaxation_time();
    t_cr = get_crossing_time();

    std::cout << "Halfmass radius " << get_halfmass_radius(p_c) << std::endl;
    std::cout << "Core radius " << get_core_radius(p_c) << std::endl;
    std::cout << "Relax. time " << t_rh << std::endl;
    std::cout << "Cross. time " << t_cr << std::endl;
}

/*
 * @fn get_center_of_density()
 *
 * @brief Function to calculate the center of density of the system
 *
 */
Point get_center_of_density()
{
    std::vector<double> p;
    std::vector<Distance> d;
    Distance dist;
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
        bb = (4.0 * M_PI * radius * radius * radius)/3.0;
        p.push_back( aa / bb );
    }

    Point density_center;
    density_center.x = 0.0f;
    density_center.y = 0.0f;
    density_center.z = 0.0f;

    dsum = 0.0f;

    for (int i = 0; i < n; i++)
    {
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

/*
 * @fn get_halfmass_radius()
 *
 * @brief Function to calculate the half mass radius using the center of density
 *
 */
double get_halfmass_radius(Point dc)
{
    float half_mass;
    double tmp, r_h;
    int i, j;

    Distance d_tmp;
    std::vector<Distance> distances(n);

    half_mass = 0;
    j = 0;

    for (i = 0; i < n; i++)
    {
        tmp = sqrt( (dc.x - h_r[i].x) * (dc.x - h_r[i].x) +
                    (dc.y - h_r[i].y) * (dc.y - h_r[i].y) +
                    (dc.z - h_r[i].z) * (dc.z - h_r[i].z) );
        d_tmp.index = i;
        d_tmp.value  = tmp;
        distances[i] = d_tmp;
    }

    std::sort(distances.begin(), distances.end());

    for (i = 0; i < n; i++)
    {
        if(half_mass >= total_mass/2.0)
        {
            j = i;
            break;
        }
        half_mass += h_m[distances[i].index];
    }

    r_h = distances[j-1].value;

    return r_h;
}

/*
 * @fn get_core_radius()
 *
 * @brief Function to calculate the core radius
 *
 */
float get_core_radius(Point pc)
{
    std::vector<Distance> d(n);
    double core_mass = 0.0;
    double radius = 0.0;
    int i;

    for (i = 0; i < n; i++)
    {
        double rx  = h_r[i].x - pc.x;
        double ry  = h_r[i].y - pc.y;
        double rz  = h_r[i].z - pc.z;
        double r   = sqrt(rx*rx + ry*ry + rz*rz);
        d[i].index = i;
        d[i].value = r;
    }

    std::sort(d.begin(), d.end());

    for (i = 0; i < n; i++)
    {
        if (core_mass > RADIUS_MASS_PORCENTAGE)
        {
            i -= 1;
            break;
        }
        core_mass += h_m[d[i].index];
    }
    radius = d[i].value;

    return radius;
}

/*
 * @fn get_crossing_time()
 *
 * @brief Function to calculate the crossing time of the system.
 *
 */
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

/*
 * @fn get_relaxation_time()
 *
 * @brief Function to calculate the half-mass relaxation time of the system.
 *
 */
float get_relaxation_time()
{
    float t_rh;
    float r_h, a, b;

    Point center = get_center_of_density();
    r_h = get_halfmass_radius(center);
    a = sqrt( (n * r_h * r_h * r_h) / ( G * (total_mass/n) ));
    b = 1/log(0.11 * n);
    t_rh = 0.138 * a * b;

    return t_rh;
}
