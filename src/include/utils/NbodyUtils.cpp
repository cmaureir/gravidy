#include "NbodyUtils.hpp"

NbodyUtils::NbodyUtils(int n, double4 *r, float m)
{
    double3 tmp = {0.0,0.0,0.0};
    this->n = n;
    this->r = r;
    this->total_mass = m;
    this->cod = tmp;
    this->energy = 0.0;
}

NbodyUtils::~NbodyUtils()
{

}

void NbodyUtils::set_energy(double e)
{
    energy = e;
}
double NbodyUtils::get_core_radius(double3 pc)
{

    std::vector<Distance> d(n);
    double core_mass = 0.0;
    double radius = 0.0;
    int i;

    for (i = 0; i < n; i++)
    {
        double rx  = r[i].x - pc.x;
        double ry  = r[i].y - pc.y;
        double rz  = r[i].z - pc.z;
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
        core_mass += r[d[i].index].w;
    }
    radius = d[i].value;

    return radius;
}

double NbodyUtils::get_half_mass_relaxation_time()
{
    float t_rh;
    float r_h, a, b;

    cod = get_center_of_density();
    r_h = get_halfmass_radius(cod);
    a = sqrt( (n * r_h * r_h * r_h) / ( G * (total_mass/n) ));
    b = 1/log(0.11 * n);
    t_rh = 0.138 * a * b;
    return t_rh;
}

double NbodyUtils::get_crossing_time()
{
    float M = total_mass;
    float Rv = (-G * M * M) / (4 * (energy));
    float Ut = sqrt( (Rv * Rv * Rv) / G * M);
    float t_cr = 2 * sqrt(2) * Ut;

    return t_cr;
}

double3 NbodyUtils::get_center_of_density()
{
    std::vector<Distance> d(n);
    std::vector<double> p(n);
    Distance empty;
    float radius, dsum;
    float aa, bb;

    empty.index = 0;
    empty.value = 0.0;

    for (int i = 0; i < n; i++)
    {
        std::fill(d.begin(), d.end(), empty);
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                double rx = (r[j].x - r[i].x);
                double ry = (r[j].y - r[i].y);
                double rz = (r[j].z - r[i].z);
                d[j].value = sqrt(rx*rx + ry*ry + rz*rz);
                d[j].index = j;
            }
        }

        std::sort(d.begin(), d.end());
        radius = d[J].value;
        aa = (J-1) * r[i].w;
        bb = (4.0 * M_PI * radius * radius * radius)/3.0;
        p[i] = aa / bb;
    }

    double3 density_center = {0.0, 0.0, 0.0};
    dsum = 0.0;

    for (int i = 0; i < n; i++)
    {
        dsum += p[i];
        density_center.x += r[i].x * p[i];
        density_center.y += r[i].y * p[i];
        density_center.z += r[i].z * p[i];
    }

    density_center.x /= dsum;
    density_center.y /= dsum;
    density_center.z /= dsum;

    return density_center;
}

double NbodyUtils::get_halfmass_radius(double3 dc)
{
    float half_mass;
    double tmp, r_h;
    int i, j;

    std::vector<Distance> distances(n);

    half_mass = 0;
    j = 0;

    for (i = 0; i < n; i++)
    {
        tmp = sqrt( (dc.x - r[i].x) * (dc.x - r[i].x) +
                    (dc.y - r[i].y) * (dc.y - r[i].y) +
                    (dc.z - r[i].z) * (dc.z - r[i].z) );
        distances[i].index = i;
        distances[i].value = tmp;
    }

    std::sort(distances.begin(), distances.end());

    for (i = 0; i < n; i++)
    {
        if(half_mass >= total_mass/2.0)
        {
            j = i;
            break;
        }
        half_mass += r[distances[i].index].w;
    }

    r_h = distances[j-1].value;

    return r_h;
}
