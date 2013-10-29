#include "Hermite4.hpp"

Hermite4::Hermite4(int n, float e2)
{
    this->n = n;
    this->e2 = e2;
}

Hermite4::~Hermite4()
{
}

void Hermite4::init_acc_jrk(double4* r, double4* v, Forces* f)
{
    int i,j;
//    #pragma omp parallel for private(j) schedule(dynamic, 24)
    for (i = INIT_PARTICLE; i < n; i++)
    {
        for (j = INIT_PARTICLE; j < n; j++)
        {
            if(i == j) continue;
            force_calculation(i, j, r, v, f);
        }

    }
}

void Hermite4::force_calculation(int i, int j, double4 *r, double4 *v, Forces *f)
{
    double rx = r[j].x - r[i].x;
    double ry = r[j].y - r[i].y;
    double rz = r[j].z - r[i].z;

    double vx = v[j].x - v[i].x;
    double vy = v[j].y - v[i].y;
    double vz = v[j].z - v[i].z;

    double r2     = rx*rx + ry*ry + rz*rz + e2;
    double rinv   = 1.0/sqrt(r2);
    double r2inv  = rinv  * rinv;
    double r3inv  = r2inv * rinv;
    double r5inv  = r2inv * r3inv;
    double mr3inv = r3inv * r[j].w;
    double mr5inv = r5inv * r[j].w;

    double rv = rx*vx + ry*vy + rz*vz;

    f[i].a[0] += (rx * mr3inv);
    f[i].a[1] += (ry * mr3inv);
    f[i].a[2] += (rz * mr3inv);

    f[i].a1[0] += (vx * mr3inv - (3 * rv ) * rx * mr5inv);
    f[i].a1[1] += (vy * mr3inv - (3 * rv ) * ry * mr5inv);
    f[i].a1[2] += (vz * mr3inv - (3 * rv ) * rz * mr5inv);
}


void Hermite4::init_dt(double &ATIME, Forces *f, double *dt, double *t)
{
    // Aarseth initial timestep
    // dt_{i} = ETA_S * sqrt( (|a|) / (|j|) )
    double tmp_dt;
    for (int i = INIT_PARTICLE; i < n; i++)
    {
        double a2 = get_magnitude(f[i].a[0],  f[i].a[1],  f[i].a[2]);
        double j2 = get_magnitude(f[i].a1[0], f[i].a1[1], f[i].a1[2]);
        tmp_dt = ETA_S * (a2/j2);

        // Adjusting to block timesteps
        // to the nearest-lower power of two
        int exp = (int)(std::ceil(log(tmp_dt)/log(2.0))-1);
        tmp_dt = pow(2,exp);

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

int Hermite4::find_particles_to_move(double ITIME, double *dt, double *t, int *move)
{

    int j = 0;
    for (int i = 0; i < n; i++)
    {
        move[i] = -1;
        if (t[i] + dt[i] == ITIME)
        {
            move[j] = i;
            j++;
        }
    }
    return j;
}

void Hermite4::save_old_acc_jrk(Forces *f, int *move, int nact, float4 *old_a, float4 *old_a1)
{
    for (int k = INIT_PARTICLE; k < nact; k++)
    {
        int i = move[k];
        old_a[i].x = f[i].a[0];
        old_a[i].y = f[i].a[1];
        old_a[i].z = f[i].a[2];

        old_a1[i].x = f[i].a1[0];
        old_a1[i].y = f[i].a1[1];
        old_a1[i].z = f[i].a1[2];
    }

}

void Hermite4::predicted_pos_vel(double ITIME, Predictor *p, double4 *r, double4 *v, Forces *f, double *t)
{
    for (int i = INIT_PARTICLE; i < n; i++)
    {
        double dt  = ITIME - t[i];
        double dt2 = (dt  * dt);
        double dt3 = (dt2 * dt);

        p[i].r[0] = (dt3/6 * f[i].a1[0]) + (dt2/2 * f[i].a[0]) + (dt * v[i].x) + r[i].x;
        p[i].r[1] = (dt3/6 * f[i].a1[1]) + (dt2/2 * f[i].a[1]) + (dt * v[i].y) + r[i].y;
        p[i].r[2] = (dt3/6 * f[i].a1[2]) + (dt2/2 * f[i].a[2]) + (dt * v[i].z) + r[i].z;

        p[i].v[0] = (dt2/2 * f[i].a1[0]) + (dt * f[i].a[0]) + v[i].x;
        p[i].v[1] = (dt2/2 * f[i].a1[1]) + (dt * f[i].a[1]) + v[i].y;
        p[i].v[2] = (dt2/2 * f[i].a1[2]) + (dt * f[i].a[2]) + v[i].z;

    }

}

void Hermite4::update_acc_jrk(int nact, int *move, double4* r, double4* v, Forces* f)
{
    int i, j, k;
    //#pragma omp parallel for private(i,j)
    for (k = 0; k < nact; k++)
    {
        i = move[k];
        f[i].a[0]  = 0.0;
        f[i].a[1]  = 0.0;
        f[i].a[2]  = 0.0;
        f[i].a1[0] = 0.0;
        f[i].a1[1] = 0.0;
        f[i].a1[2] = 0.0;

        //#pragma omp parallel for
        for (j = INIT_PARTICLE; j < n; j++)
        {
            if(i == j) continue;
            force_calculation(i, j, r, v, f);
        }
    }

}

void Hermite4::correction_pos_vel(double ITIME, int nact, int *move, double4 *r, double4 *v, Forces *f, double *t, double *dt, Predictor *p, float4 *old_a, float4 *old_a1, float4 *a3, float4 *a2, float eta)
{
    for (int k = 0; k < nact; k++)
    {
        int i = move[k];

        double dt1 = dt[i];
        double dt2 = dt1 * dt1;
        double dt3 = dt2 * dt1;
        double dt4 = dt2 * dt2;
        double dt5 = dt4 * dt1;

        // Acceleration 2nd derivate
        a2[i].x = (-6 * (old_a[i].x - f[i].a[0] ) - dt1 * (4 * old_a1[i].x + 2 * f[i].a1[0]) ) / dt2;
        a2[i].y = (-6 * (old_a[i].y - f[i].a[1] ) - dt1 * (4 * old_a1[i].y + 2 * f[i].a1[1]) ) / dt2;
        a2[i].z = (-6 * (old_a[i].z - f[i].a[2] ) - dt1 * (4 * old_a1[i].z + 2 * f[i].a1[2]) ) / dt2;

        // Acceleration 3rd derivate
        a3[i].x = (12 * (old_a[i].x - f[i].a[0] ) + 6 * dt1 * (old_a1[i].x + f[i].a1[0]) ) / dt3;
        a3[i].y = (12 * (old_a[i].y - f[i].a[1] ) + 6 * dt1 * (old_a1[i].y + f[i].a1[1]) ) / dt3;
        a3[i].z = (12 * (old_a[i].z - f[i].a[2] ) + 6 * dt1 * (old_a1[i].z + f[i].a1[2]) ) / dt3;

        // Correcting position
        r[i].x = p[i].r[0] + (dt4/24)*a2[i].x + (dt5/120)*a3[i].x;
        r[i].y = p[i].r[1] + (dt4/24)*a2[i].y + (dt5/120)*a3[i].y;
        r[i].z = p[i].r[2] + (dt4/24)*a2[i].z + (dt5/120)*a3[i].z;

        // Correcting velocity
        v[i].x = p[i].v[0] + (dt3/6)*a2[i].x + (dt4/24)*a3[i].x;
        v[i].y = p[i].v[1] + (dt3/6)*a2[i].y + (dt4/24)*a3[i].y;
        v[i].z = p[i].v[2] + (dt3/6)*a2[i].z + (dt4/24)*a3[i].z;

        t[i] = ITIME;
        double normal_dt  = get_timestep_normal(i, a2, a3, dt, f, eta);
        normal_dt = normalize_dt(normal_dt, dt[i], t[i], i);
        dt[i] = normal_dt;

    }

}
