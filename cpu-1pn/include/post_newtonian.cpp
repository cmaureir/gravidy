#include "post_newtonian.hpp"

void first_pn_calculation(int i)
{
    int c2 = C*C; // CL2
    double general1[3]  = {0.0, 0.0, 0.0};
    double general1d[3] = {0.0, 0.0, 0.0};

    double rx = h_r[i].x - h_r[0].x;
    double ry = h_r[i].y - h_r[0].y;
    double rz = h_r[i].z - h_r[0].z;
    double vx = h_v[i].x - h_v[0].x;
    double vy = h_v[i].y - h_v[0].y;
    double vz = h_v[i].z - h_v[0].z;

    double r2 = rx*rx + ry*ry + rz*rz; //RIJ2
    double rinv = 1/sqrt(r2); //RIJ but ^{-1}
    double r2inv = rinv  * rinv;
    double r3inv = r2inv  * rinv;
    double r4inv = r2inv  * r2inv;

    double nx = rx * rinv; // NVEC
    double ny = ry * rinv;
    double nz = rz * rinv;

    double at1x = -h_m[0] * nx * r2inv; // AT1
    double at1y = -h_m[0] * ny * r2inv;
    double at1z = -h_m[0] * nz * r2inv;

    double at2x = h_m[i] * nx * r2inv; // AT2
    double at2y = h_m[i] * ny * r2inv;
    double at2z = h_m[i] * nz * r2inv;

    double rv = rx*vx + ry*vy + rz*vz;
    double rp = rv * rinv; // RP

    double v1_2 = h_v[i].x * h_v[i].x + h_v[i].y * h_v[i].y + h_v[i].z * h_v[i].z; // V1_2
    double v2_2 = h_v[0].x * h_v[0].x + h_v[0].y * h_v[0].y + h_v[0].z * h_v[0].z; // V2_2

    double v1v2 = h_v[0].x * h_v[i].x + h_v[0].y * h_v[i].y + h_v[0].z * h_v[i].z; // V1V2

    // V1A
    double v1a = 2.0 *  (h_v[i].x * at1x +
                         h_v[i].y * at1y +
                         h_v[i].z * at1z);
    // V2A
    double v2a = 2.0 *  (h_v[0].x * at2x +
                          h_v[0].y * at2y +
                          h_v[0].z * at2z);
    // NV1
    double nv1 = nx * h_v[i].x + ny * h_v[i].y + nz * h_v[i].z;
    // NV2
    double nv2 = nx * h_v[0].x + ny * h_v[0].y + nz * h_v[0].z;
    double nv2_2 = nv2*nv2;



    double ndotx = (vx - nx * rp) * rinv;
    double ndoty = (vy - ny * rp) * rinv;
    double ndotz = (vz - nz * rp) * rinv;

    double nv1dot = ndotx * h_v[i].x + ndoty * h_v[i].y + ndotz * h_v[i].z +\
                    nx * at1x + ny * at1y + nz * at1z;

    double nv2dot = ndotx * h_v[0].x + ndoty * h_v[0].y + ndotz * h_v[0].z +\
                    nx * at2x + ny * at2y + nz * at2z;

    // G is ommited because is 1
    general1[0]= nx * ((5.0*h_m[i]*h_m[0] + 4.0*h_m[0]*h_m[0]) * r3inv +\
                 h_m[0] * r2inv * (1.5*nv2*nv2 - v1_2 + 4.0*v1v2 -\
                 2.0*v2_2)) + \
                 vx * h_m[0] * r2inv  * (4.0*nv1-3*nv2);

    general1[1]= ny * ((5.0*h_m[i]*h_m[0] + 4.0*h_m[0]*h_m[0]) * r3inv +\
                 h_m[0] * r2inv * (1.5*nv2*nv2 - v1_2 + 4.0*v1v2 -\
                 2.0*v2_2)) + \
                 vy * h_m[0] * r2inv  * (4.0*nv1-3*nv2);

    general1[2]= nz * ((5.0*h_m[i]*h_m[0] + 4.0*h_m[0]*h_m[0]) * r3inv +\
                 h_m[0] * r2inv * (1.5*nv2*nv2 - v1_2 + 4.0*v1v2 -\
                 2.0*v2_2)) + \
                 vz * h_m[0] * r2inv  * (4.0*nv1-3*nv2);


    general1d[0] = nx * (-3.0 * r4inv * rp * (5.0 *h_m[i]*h_m[0] + 4.0 * \
                   h_m[0]*h_m[0]) - 2.0 * h_m[0] * r3inv * rp * (1.5*nv2_2 - v1_2 + 4.0*v1v2-2.0*v2_2)+\
        h_m[0] * r2inv * (3.0*nv2 * nv2dot - v1a + 4.0 * ((at1x*h_v[0].x + at1y*h_v[0].y + at1z*h_v[0].z)+\
        (at2x*h_v[i].x + at2y*h_v[i].y + at2z*h_v[i].z)) - \
        2.0*v2a)) + ndotx * ((5.0*h_m[i] * h_m[0] + 4.0 * h_m[0] * h_m[0]) * r3inv + \
        h_m[0] * r2inv * (1.5*nv2_2 - v1_2 + 4.0*v1v2 - 2.0*v2_2)) \
        + vx * (-2.0*h_m[0] * rp * r3inv * (4.0 * nv1 - 3.0 *nv2)+ h_m[0] * r2inv * (4.0 * nv1dot \
        -3.0 * nv2dot)) + ( at1x - at2x ) * h_m[0] * r2inv * (4.0 * nv1 - 3*nv2);

    general1d[1] = ny * (-3.0 * r4inv * rp * (5.0 *h_m[i]*h_m[0] + 4.0 * \
                   h_m[0]*h_m[0]) - 2.0 * h_m[0] * r3inv * rp * (1.5*nv2_2 - v1_2 + 4.0*v1v2-2.0*v2_2)+\
        h_m[0] * r2inv * (3.0*nv2 * nv2dot - v1a + 4.0 * ((at1x*h_v[0].x + at1y*h_v[0].y + at1z*h_v[0].z)+\
        (at2x*h_v[i].x + at2y*h_v[i].y + at2z*h_v[i].z)) - \
        2.0*v2a)) + ndoty * ((5.0*h_m[i] * h_m[0] + 4.0 * h_m[0] * h_m[0]) * r3inv + \
        h_m[0] * r2inv * (1.5*nv2_2 - v1_2 + 4.0*v1v2 - 2.0*v2_2)) \
        + vy * (-2.0*h_m[0] * rp * r3inv * (4.0 * nv1 - 3.0 *nv2)+ h_m[0] * r2inv * (4.0 * nv1dot \
        -3.0 * nv2dot)) + ( at1y - at2y ) * h_m[0] * r2inv * (4.0 * nv1 - 3*nv2);

    general1d[2] = nz * (-3.0 * r4inv * rp * (5.0 *h_m[i]*h_m[0] + 4.0 * \
                   h_m[0]*h_m[0]) - 2.0 * h_m[0] * r3inv * rp * (1.5*nv2_2 - v1_2 + 4.0*v1v2-2.0*v2_2)+\
        h_m[0] * r2inv * (3.0*nv2 * nv2dot - v1a + 4.0 * ((at1x*h_v[0].x + at1y*h_v[0].y + at1z*h_v[0].z)+\
        (at2x*h_v[i].x + at2y*h_v[i].y + at2z*h_v[i].z)) - \
        2.0*v2a)) + ndotz * ((5.0*h_m[i] * h_m[0] + 4.0 * h_m[0] * h_m[0]) * r3inv + \
        h_m[0] * r2inv * (1.5*nv2_2 - v1_2 + 4.0*v1v2 - 2.0*v2_2)) \
        + vz * (-2.0*h_m[0] * rp * r3inv * (4.0 * nv1 - 3.0 *nv2)+ h_m[0] * r2inv * (4.0 * nv1dot \
        -3.0 * nv2dot)) + ( at1z - at2z ) * h_m[0] * r2inv * (4.0 * nv1 - 3*nv2);

    for (int k=0;k<3; k++)
    {
        h_f[i].a[k]  += general1[k]/c2;
        h_f[i].a1[k] += general1d[k]/c2;
    }
    //printf("1PNa:  %.8f %.8f %.8f\n", general1[0], general1[1], general1[2]);
    //printf("a_{%d}:  %.8f %.8f %.8f\n", i, h_f[i].a[0], h_f[i].a[1], h_f[i].a[2]);
    //printf("---:  %.8f %.8f %.8f\n", general1[0]/h_f[i].a[0], general1[1]/h_f[i].a[1], general1[2]/h_f[i].a[2]);

    //printf("1PNb:  %.8f %.8f %.8f\n", general1d[0], general1d[1], general1d[2]);
    //printf("j_{%d}:  %.8f %.8f %.8f\n", i, h_f[i].a1[0], h_f[i].a1[1], h_f[i].a1[2]);
    //printf("---:  %.8f %.8f %.8f\n", general1d[0]/h_f[i].a1[0], general1d[1]/h_f[i].a1[1], general1d[2]/h_f[i].a1[2]);
    //getchar();
}


void update_acc_jrk_1pn(int total)
{
    int i, j;
    //#pragma omp parallel for private(i,j)
    for (i = 0; i < total; i++)
    {
        j = h_move[i];
        first_pn_calculation(j);
    }
}
