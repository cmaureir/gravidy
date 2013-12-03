#include "Hermite4.hpp"
#include <iostream>

Hermite4::Hermite4(int n, double e2, float eta)
{
    this->n   = n;
    this->e2  = e2;
    this->eta = eta;
}

Hermite4::~Hermite4()
{
    // Nothing
}

int Hermite4::find_particles_to_move(int *move, double ITIME, double *dt, double *t)
{

    int j = 0;
    for (int i = 0; i < n; i++)
    {
        move[i] = -1;
        double tmp_time = t[i] + dt[i];
        if(std::fabs(ITIME - tmp_time) < 2*std::numeric_limits<double>::epsilon())
        {
            move[j] = i;
            j++;
        }
    }
    return j;
}

void Hermite4::next_integration_time(double &ATIME, double *dt, double *t)
{
    // Big number to find the minimum
    //ATIME = 1.0e10;
    ATIME = t[0] + dt[0];
    for (int i = 1; i < n; i++)
    {
        double time = t[i] + dt[i];
        if(time < ATIME)
        {
            ATIME = time;
        }
    }
}

void Hermite4::init_dt(double &ATIME, double *dt, double *t, Forces *f)
{
    // Aarseth initial timestep
    // dt_{i} = ETA_S * sqrt( (|a|) / (|j|) )
    double tmp_dt;
    for (int i = 0; i < n; i++)
    {
        double a2 = get_magnitude(f[i].a[0],  f[i].a[1],  f[i].a[2]);
        double j2 = get_magnitude(f[i].a1[0], f[i].a1[1], f[i].a1[2]);
        tmp_dt = ETA_S * (a2/j2);

        // Adjusting to block timesteps
        // to the nearest-lower power of two
        int exp = (int)(std::ceil(log(tmp_dt)/log(2.0))-1);
        tmp_dt = pow(2,exp);
        //2 << (exp - 1);

        if (tmp_dt < D_TIME_MIN)
            tmp_dt = D_TIME_MIN;
        else if (tmp_dt > D_TIME_MAX)
            tmp_dt = D_TIME_MAX;

        dt[i] = tmp_dt;
        t[i] = 0.0;

        // Obtaining the first integration time
        if(tmp_dt < ATIME)
            ATIME = tmp_dt;
    }
}

void Hermite4::save_old_acc_jrk(int nact, int *move, Forces *old, Forces *f)
{
    for (int k = 0; k < nact; k++)
    {
        int i = move[k];
        old[i].a[0]  = f[i].a[0];
        old[i].a[1]  = f[i].a[1];
        old[i].a[2]  = f[i].a[2];

        old[i].a1[0] = f[i].a1[0];
        old[i].a1[1] = f[i].a1[1];
        old[i].a1[2] = f[i].a1[2];
    }

}

