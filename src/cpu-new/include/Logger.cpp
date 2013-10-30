#include "Logger.hpp"

void Logger::print_all(double ITIME, int n, double4 *r, double4 *v, Forces *f, double *dt)
{
    for (int i = 0; i < n; i++) {
    printf("% .10f % 6d % .10f % .10f % .10f % .10f % .10f % .10f % .10f % .10f % .10f % .10f % .10f % .10f % .10f\t% .10f\n",
            ITIME, i, r[i].w,
            r[i].x, r[i].y, r[i].z,
            v[i].x, v[i].y, v[i].z,
            f[i].a[0],  f[i].a[1],  f[i].a[2],
            f[i].a1[0], f[i].a1[1], f[i].a1[2],
            dt[i]);
    }
}
