#include "Hermite4Kepler.hpp"

void Hermite4Kepler::force_calculation(int i, int j, Predictor *p, Forces *f)
{
    double rx = p[j].r[0] - p[i].r[0];
    double ry = p[j].r[1] - p[i].r[1];
    double rz = p[j].r[2] - p[i].r[2];

    double vx = p[j].v[0] - p[i].v[0];
    double vy = p[j].v[1] - p[i].v[1];
    double vz = p[j].v[2] - p[i].v[2];

    double r2     = rx*rx + ry*ry + rz*rz + e2;
    double rinv   = 1.0/sqrt(r2);
    double r2inv  = rinv  * rinv;
    double r3inv  = r2inv * rinv;
    double r5inv  = r2inv * r3inv;
    double mr3inv = r3inv * p[j].m;
    double mr5inv = r5inv * p[j].m;

    double rv = rx*vx + ry*vy + rz*vz;

    f[i].a[0] += (rx * mr3inv);
    f[i].a[1] += (ry * mr3inv);
    f[i].a[2] += (rz * mr3inv);

    f[i].a1[0] += (vx * mr3inv - (3 * rv ) * rx * mr5inv);
    f[i].a1[1] += (vy * mr3inv - (3 * rv ) * ry * mr5inv);
    f[i].a1[2] += (vz * mr3inv - (3 * rv ) * rz * mr5inv);
}

void Hermite4Kepler::init_acc_jrk(Predictor *p, Forces *f)
{
    int i,j;
    #pragma omp parallel for private(j) schedule(dynamic, 24)
    for (i = INIT_PARTICLE; i < n; i++)
    {
        for (j = INIT_PARTICLE; j < n; j++)
        {
            if(i == j) continue;
            force_calculation(i, j, p, f);
        }
    }
}

void Hermite4Kepler::update_acc_jrk(int nact, int *move, Predictor *p, Forces *f, Gtime &gtime)
{
    gtime.update_ini = omp_get_wtime();
    int i, j, k;
    #pragma omp parallel for private(i,j)
    for (k = 0; k < nact; k++)
    {
        i = move[k];
        f[i].a[0]  = 0.0;
        f[i].a[1]  = 0.0;
        f[i].a[2]  = 0.0;
        f[i].a1[0] = 0.0;
        f[i].a1[1] = 0.0;
        f[i].a1[2] = 0.0;

        #pragma omp parallel for
        for (j = INIT_PARTICLE; j < n; j++)
        {
            if(i == j) continue;
            force_calculation(i, j, p, f);
        }
    }
    gtime.update_end += omp_get_wtime() - gtime.update_ini;

}
