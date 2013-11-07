#include "extra_utils.hpp"

double get_magnitude(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}

double get_timestep_normal(int i, double4 *a2, double4 *a3, double *dt, Forces *f, double eta)
{
    // Calculating a_{1,i}^{(2)} = a_{0,i}^{(2)} + dt * a_{0,i}^{(3)}
    double ax1_2 = a2[i].x + dt[i] * a3[i].x;
    double ay1_2 = a2[i].y + dt[i] * a3[i].y;
    double az1_2 = a2[i].z + dt[i] * a3[i].z;

    // |a_{1,i}|
    double abs_a1 = get_magnitude(f[i].a[0], f[i].a[1], f[i].a[2]);
    // |j_{1,i}|
    double abs_j1 = get_magnitude(f[i].a1[0], f[i].a1[1], f[i].a1[2]);
    // |j_{1,i}|^{2}
    double abs_j12  = abs_j1 * abs_j1;
    // a_{1,i}^{(3)} = a_{0,i}^{(3)} because the 3rd-order interpolation
    double abs_a1_3 = get_magnitude(a3[i].x, a3[i].y, a3[i].z);
    // |a_{1,i}^{(2)}|
    double abs_a1_2 = get_magnitude(ax1_2, ay1_2, az1_2);
    // |a_{1,i}^{(2)}|^{2}
    double abs_a1_22  = abs_a1_2 * abs_a1_2;

    double normal_dt = sqrt(eta * ((abs_a1 * abs_a1_2 + abs_j12) / (abs_j1 * abs_a1_3 + abs_a1_22)));

    return normal_dt;
}

double normalize_dt(double new_dt, double old_dt, double t, int i)
{
    if (new_dt <= old_dt/8)
    {
        new_dt = D_TIME_MIN;
    }
    else if ( old_dt/8 < new_dt && new_dt <= old_dt/4)
    {
        new_dt = old_dt / 8;
    }
    else if ( old_dt/4 < new_dt && new_dt <= old_dt/2)
    {
        new_dt = old_dt / 4;
    }
    else if ( old_dt/2 < new_dt && new_dt <= old_dt)
    {
        new_dt = old_dt / 2;
    }
    else if ( old_dt < new_dt && new_dt <= old_dt * 2)
    {
        new_dt = old_dt;
    }
    else if (2 * old_dt < new_dt)
    {
        double val = t/(2 * old_dt);
        //float val = t/(2 * old_dt);
        if(std::ceil(val) == val)
        {
            new_dt = 2.0 * old_dt;
        }
        else
        {
            new_dt = old_dt;
        }
    }
    else
    {
        //std::cerr << "this will never happen...I promise" << std::endl;
        new_dt = old_dt;
    }

    //if (new_dt <= D_TIME_MIN)
    if (new_dt < D_TIME_MIN)
    {
        new_dt = D_TIME_MIN;
    }
    //else if (new_dt >= D_TIME_MAX)
    else if (new_dt > D_TIME_MAX)
    {
        new_dt = D_TIME_MAX;
    }

    return new_dt;
}

