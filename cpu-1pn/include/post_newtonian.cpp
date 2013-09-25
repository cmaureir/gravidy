#include "post_newtonian.hpp"

void first_pn_calculation(int i)
{
    int c2 = C*C; // CL2
    double general1[3]  = {0.0, 0.0, 0.0};
    double general1d[3] = {0.0, 0.0, 0.0};

    double rx = h_r[0].x - h_r[i].x;
    double ry = h_r[0].y - h_r[i].y;
    double rz = h_r[0].z - h_r[i].z;

    double vx = h_v[0].x - h_v[i].x;
    double vy = h_v[0].y - h_v[i].y;
    double vz = h_v[0].z - h_v[i].z;

    double r2 = rx*rx + ry*ry + rz*rz; //RIJ2
    double rinv = 1/sqrt(r2); //RIJ but ^{-1}
    double r2inv = rinv  * rinv;
    double r3inv = r2inv  * rinv;
    double r4inv = r2inv  * r2inv;

    double nx = rx * rinv; // NVEC
    double ny = ry * rinv;
    double nz = rz * rinv;


    double mnrbhx = -h_m[0] * nx * r2inv; // AT1
    double mnrbhy = -h_m[0] * ny * r2inv;
    double mnrbhz = -h_m[0] * nz * r2inv;

    double mnrpx = -h_m[i] * nx * r2inv; // AT2
    double mnrpy = -h_m[i] * ny * r2inv;
    double mnrpz = -h_m[i] * nz * r2inv;

    double rv = rx*vx + ry*vy + rz*vz;
    double rvr = rv * rinv; // RP

    double v12_2 = vx*vx + vy*vy + vz*vz; //V1_V22
    double vp2 = h_v[i].x * h_v[i].x + h_v[i].y * h_v[i].y + h_v[i].z * h_v[i].z; // V1_2
    double vbh2 = h_v[0].x * h_v[0].x + h_v[0].y * h_v[0].y + h_v[0].z * h_v[0].z; // V2_2

    double v1v2 = h_v[0].x * h_v[i].x + h_v[0].y * h_v[i].y + h_v[0].z * h_v[i].z; // V1V2

    // V12A
    double vmnr = 2.0 * (vx * (mnrbhx - mnrpx) +
                          vy * (mnrbhy - mnrpy) +
                          vz * (mnrbhz - mnrpz));
    // V1A
    double vpa = 2.0 * (h_v[i].x * mnrpx +
                         h_v[i].y * mnrpy +
                         h_v[i].z * mnrpz);
    // V2A
    double vbha = 2.0 * (h_v[0].x * mnrbhx +
                          h_v[0].y * mnrbhy +
                          h_v[0].z * mnrbhz);
    // NV1
    double nv1 = nx * h_v[i].x + ny * h_v[i].y + nz * h_v[i].z;
    // NV2
    double nv2 = nx * h_v[0].x + ny * h_v[0].y + nz * h_v[0].z;
    double nv2_2 = nv2*nv2;
    // NV
    double nv  = nx * vx + ny * vy + nz * vz;

    double ndotx = vx - nx * rvr * rinv;
    double ndoty = vy - ny * rvr * rinv;
    double ndotz = vz - nz * rvr * rinv;

    double nv1dot = ndotx * h_v[i].x + ndoty * h_v[i].y + ndotz * h_v[i].z +\
                    nx * mnrbhx + ny * mnrbhy + nz * mnrbhz;

    double nv2dot = ndotx * h_v[0].x + ndoty * h_v[0].y + ndotz * h_v[0].z +\
                    nx * mnrpx + ny * mnrpy + nz * mnrpz;

    double nvdot = ndotx * vx + ndoty * vy + ndotz * vz + \
                   nx * (mnrbhx - mnrpx) + ny * (mnrbhy - mnrpy) + nz * (mnrbhz - mnrpz);

    // G is ommited because is 1
    general1[0]= h_m[0] * r2inv * (nx * (-vp2 - 2*vbh2 + 4*v1v2 + 1.5*nv2_2 +\
                 5*(h_m[i] * rinv) + 4*(h_m[0] * rinv) ) + vx*(4*nv1 - 3 * nv2));

    general1[1]= h_m[0] * r2inv * (ny * (-vp2 - 2*vbh2 + 4*v1v2 + 1.5*nv2_2 +\
                 5*(h_m[i] * rinv) + 4*(h_m[0] * rinv) ) + vy*(4*nv1 - 3 * nv2));

    general1[2]= h_m[0] * r2inv * (nz * (-vp2 - 2*vbh2 + 4*v1v2 + 1.5*nv2_2 +\
                 5*(h_m[i] * rinv) + 4*(h_m[0] * rinv) ) + vz*(4*nv1 - 3 * nv2));


    general1d[0] = nx * (-3.0 * r4inv * rvr * (5.0 *h_m[i]*h_m[0] + 4.0 * \
                   h_m[0]*h_m[0]) - 2.0 * h_m[0] * r3inv * rvr * (1.5*nv2_2 - vp2 + 4.0*v1v2-2.0*vbh2)+\
        h_m[0] * r2inv * (3.0*nv2 * nv2dot - vpa + 4.0 * (mnrbhx*h_v[0].x + mnrbhy*h_v[0].y + mnrbhz*h_v[0].z)+\
        (mnrpx*h_v[i].x + mnrpy*h_v[i].y + mnrpz*h_v[i].z)) - \
        2.0*vbha) + ndotx * ((5.0*h_m[i] * h_m[0] + 4.0 * h_m[0] * h_m[0]) * r3inv + \
        h_m[0] * r2inv * (1.5*nv2_2 - vp2 + 4.0*v1v2 - 2.0*vbh2 + 4.0 * v1v2 - 2*vbh2)) \
        + vx * (-2.0*h_m[0] * rvr * r3inv * (4.0 * nv1 - 3.0 *nv2)+ h_m[0] * r2inv * (4.0 * nv1dot \
        -3.0 * nv2dot)) + ( mnrbhx - mnrpx ) * h_m[0] * r2inv * (4.0 * nv1 - 3*nv2);

    general1d[1] = ny * (-3.0 * r4inv * rvr * (5.0 *h_m[i]*h_m[0] + 4.0 * \
                   h_m[0]*h_m[0]) - 2.0 * h_m[0] * r3inv * rvr * (1.5*nv2_2 - vp2 + 4.0*v1v2-2.0*vbh2)+\
        h_m[0] * r2inv * (3.0*nv2 * nv2dot - vpa + 4.0 * (mnrbhy*h_v[0].y + mnrbhy*h_v[0].y + mnrbhz*h_v[0].z)+\
        (mnrpy*h_v[i].y + mnrpy*h_v[i].y + mnrpz*h_v[i].z)) - \
        2.0*vbha) + ndoty * ((5.0*h_m[i] * h_m[0] + 4.0 * h_m[0] * h_m[0]) * r3inv + \
        h_m[0] * r2inv * (1.5*nv2_2 - vp2 + 4.0*v1v2 - 2.0*vbh2 + 4.0 * v1v2 - 2*vbh2)) \
        + vy * (-2.0*h_m[0] * rvr * r3inv * (4.0 * nv1 - 3.0 *nv2)+ h_m[0] * r2inv * (4.0 * nv1dot \
        -3.0 * nv2dot)) + ( mnrbhy - mnrpy ) * h_m[0] * r2inv * (4.0 * nv1 - 3*nv2);

    general1d[2] = nz * (-3.0 * r4inv * rvr * (5.0 *h_m[i]*h_m[0] + 4.0 * \
                   h_m[0]*h_m[0]) - 2.0 * h_m[0] * r3inv * rvr * (1.5*nv2_2 - vp2 + 4.0*v1v2-2.0*vbh2)+\
        h_m[0] * r2inv * (3.0*nv2 * nv2dot - vpa + 4.0 * (mnrbhz*h_v[0].z + mnrbhy*h_v[0].y + mnrbhz*h_v[0].z)+\
        (mnrpz*h_v[i].z + mnrpy*h_v[i].y + mnrpz*h_v[i].z)) - \
        2.0*vbha) + ndotz * ((5.0*h_m[i] * h_m[0] + 4.0 * h_m[0] * h_m[0]) * r3inv + \
        h_m[0] * r2inv * (1.5*nv2_2 - vp2 + 4.0*v1v2 - 2.0*vbh2 + 4.0 * v1v2 - 2*vbh2)) \
        + vz * (-2.0*h_m[0] * rvr * r3inv * (4.0 * nv1 - 3.0 *nv2)+ h_m[0] * r2inv * (4.0 * nv1dot \
        -3.0 * nv2dot)) + ( mnrbhz - mnrpz ) * h_m[0] * r2inv * (4.0 * nv1 - 3*nv2);

    for (int k=0;k<3; k++)
    {
        h_f[i].a[k]  += general1[k]/c2;
        h_f[i].a1[k] += general1d[k]/c2;
    }
    //printf("1PNa:  %.8f %.8f %.8f\n", general1[0], general1[1], general1[2]);
    //printf("a_{%d}:  %.8f %.8f %.8f\n", i, h_f[i].a[0], h_f[i].a[1], h_f[i].a[2]);

    //printf("1PNb:  %.8f %.8f %.8f\n", general1d[0], general1d[1], general1d[2]);
    //printf("j_{%d}:  %.8f %.8f %.8f\n", i, h_f[i].a1[0], h_f[i].a1[1], h_f[i].a1[2]);
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
