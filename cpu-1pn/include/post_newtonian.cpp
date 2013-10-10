#include "post_newtonian.hpp"

void first_pn_calculation(int i)
{

    //int c2 = C*C; // CL2
    int c2 = 5*5;
    int c4 = c2*c2;
    int c5 = c4*5;
    double pn1[3]  = {0.0, 0.0, 0.0};
    double pn1d[3] = {0.0, 0.0, 0.0};
    double pn2[3]  = {0.0, 0.0, 0.0};
    double pn2d[3] = {0.0, 0.0, 0.0};
    double pn25[3]  = {0.0, 0.0, 0.0};
    double pn25d[3] = {0.0, 0.0, 0.0};


    //double r[0] = h_r[i].x - h_r[0].x;
    //double r[1] = h_r[i].y - h_r[0].y;
    //double r[2] = h_r[i].z - h_r[0].z;
    double r[3];
    r[0] = 1.0;
    r[1] = 0.5;
    r[2] = 0.1;
    h_v[i].x = 1.0;
    h_v[i].y = 0.2;
    h_v[i].z = 1.5;
    h_v[0].x = 0.3;
    h_v[0].y = 1.0;
    h_v[0].z = 0.5;
    h_m[i] = 1.0;
    h_m[0] = 3.0;

    double v[3];
    v[0] = h_v[i].x - h_v[0].x;
    v[1] = h_v[i].y - h_v[0].y;
    v[2] = h_v[i].z - h_v[0].z;

    double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2]; //RIJ2
    double rinv = 1/sqrt(r2); //RIJ but ^{-1}
    double r2inv = rinv  * rinv;
    double r3inv = r2inv  * rinv;
    double r4inv = r2inv  * r2inv;
    double r5inv = r3inv  * r2inv;

    double v1_v22 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];

    double nvec[3];
    nvec[0] = r[0] * rinv; // NVEC
    nvec[1] = r[1] * rinv;
    nvec[2] = r[2] * rinv;

    double at1[3];
    at1[0] = -h_m[0] * nvec[0] * r2inv; // AT1
    at1[1] = -h_m[0] * nvec[1] * r2inv;
    at1[2] = -h_m[0] * nvec[2] * r2inv;

    double at2[3];
    at2[0] = h_m[i] * nvec[0] * r2inv; // AT2
    at2[1] = h_m[i] * nvec[1] * r2inv;
    at2[2] = h_m[i] * nvec[2] * r2inv;

    double rv = r[0]*v[0] + r[1]*v[1] + r[2]*v[2];
    double rp = rv * rinv; // RP

    // V!_2
    double v1_2 = h_v[i].x * h_v[i].x + h_v[i].y * h_v[i].y + h_v[i].z * h_v[i].z;
    // V2_2
    double v2_2 = h_v[0].x * h_v[0].x + h_v[0].y * h_v[0].y + h_v[0].z * h_v[0].z;
    // V1V2
    double v1v2 = h_v[0].x * h_v[i].x + h_v[0].y * h_v[i].y + h_v[0].z * h_v[i].z;
    // V1A
    double v1a = 2.0 *  (h_v[i].x * at1[0] + h_v[i].y * at1[1] + h_v[i].z * at1[2]);
    // V2A
    double v2a = 2.0 *  (h_v[0].x * at2[0] + h_v[0].y * at2[1] + h_v[0].z * at2[2]);
    // V1V2DOT
    double v1v2dot = (at1[0] * h_v[0].x + at1[1] * h_v[0].y + at1[2] * h_v[0].z) +
                     (at2[0] * h_v[i].x + at2[1] * h_v[i].y + at2[2] * h_v[i].z);

    double v12a = 2. * (v[0] * (at1[0] - at2[0]) +
                        v[1] * (at1[1] - at2[1]) +
                        v[2] * (at1[2] - at2[2]));
    // NV1
    double nv1 = nvec[0] * h_v[i].x + nvec[1] * h_v[i].y + nvec[2] * h_v[i].z;
    double nv1_2 = nv1 * nv1;
    // NV2
    double nv2 = nvec[0] * h_v[0].x + nvec[1] * h_v[0].y + nvec[2] * h_v[0].z;
    double nv2_2 = nv2 * nv2;
    double nv2_3 = nv2_2 * nv2;
    double nv2_4 = nv2_2 * nv2_2;

    double nvdot = (nvec[0] * v[0] + nvec[1] * v[1] + nvec[2] * v[2]);


    double ndot[3];
    ndot[0] = (v[0] - nvec[0] * rp) * rinv;
    ndot[1] = (v[1] - nvec[1] * rp) * rinv;
    ndot[2] = (v[2] - nvec[2] * rp) * rinv;

    double nv1dot = ndot[0] * h_v[i].x + ndot[1] * h_v[i].y + ndot[2] * h_v[i].z + \
                    nvec[0] * at1[0] + nvec[1] * at1[1] + nvec[2] * at1[2];

    double nv2dot = ndot[0] * h_v[0].x + ndot[1] * h_v[0].y + ndot[2] * h_v[0].z + \
                    nvec[0] * at2[0] + nvec[1] * at2[1] + nvec[2] * at2[2];

    double at1v = at1[0]*h_v[0].x + at1[1]*h_v[0].y + at1[2]*h_v[0].z;
    double at2v = at2[0]*h_v[i].x + at2[1]*h_v[i].y + at2[2]*h_v[i].z;

    double mbh2 = h_m[0] * h_m[0];
    double mp2  = h_m[i] * h_m[i];

    double mpbh = h_m[i] * h_m[0];
    double mpbh2 = h_m[i] * mbh2;
    double mp2bh = mp2 * h_m[0];

    double mbh3 = h_m[0] * mbh2;

    // G is ommited because is 1
    for (int k = 0; k < 3; k++)
    {
        pn1[k] = nvec[k] * ((5.0 * mpbh + 4.0 * mbh2) * r3inv + h_m[0] * r2inv \
                 * (1.5 * nv2 * nv2 - v1_2 + 4.0 * v1v2 - 2.0 * v2_2)) + v[k]  \
                 * h_m[0] * r2inv  * (4.0 * nv1 - 3 * nv2);

        pn1d[k] = nvec[k] * (-3.0 * r4inv * rp * (5.0 * mpbh + 4.0 * mbh2) -    \
                  2.0 * h_m[0] * r3inv * rp * (1.5 * nv2_2 - v1_2 + 4.0 * v1v2  \
                  - 2.0 * v2_2) + h_m[0] * r2inv * (3.0 * nv2 * nv2dot - v1a    \
                  + 4.0 * (at1v + at2v) - 2.0 * v2a)) + ndot[k] * ((5.0 * mpbh  \
                  + 4.0 * mbh2) * r3inv + h_m[0] * r2inv * (1.5 * nv2_2 - v1_2  \
                  + 4.0 * v1v2 - 2.0 * v2_2)) + v[k] * (-2.0 * h_m[0] * rp *    \
                  r3inv * (4.0 * nv1 - 3.0 * nv2) + h_m[0] * r2inv * (4.0 *     \
                  nv1dot - 3.0 * nv2dot)) + (at1[k] - at2[k] ) * h_m[0] * r2inv \
                  * (4.0 * nv1 - 3*nv2);

        pn2[k] = nvec[k] * ((-57./4.* mp2bh - 69./2. * mpbh2 - 9.* mbh3) * r4inv + \
                 h_m[0] * r2inv *(-15./8. * nv2_4 + 1.5 * nv2_2 * v1_2 - 6. * nv2_2\
                 * v1v2 - 2. * v1v2 * v1v2 + 4.5 * nv2_2 * v2_2 + 4. * v1v2 * v2_2 \
                 - 2. * v2_2 * v2_2) + mpbh * r3inv * (19.5 * nv1_2 - 39. * nv1 *  \
                 nv2 + 8.5 * nv2 * nv2 - 15./4. * v1_2 - 2.5 * v1v2 + 5./4. * v2_2)\
                 + h_m[0] * h_m[0] * r3inv * (2. * nv1_2 - 4. * nv1 * nv2 - 6. *   \
                 nv2_2 - 8. * v1v2 + 4. * v2_2)) + v[k] * (mbh2 * r3inv * (-2. *   \
                 nv1 - 2. * nv2) + mpbh * r3inv *(-63./4. * nv1 + 55./4. * nv2) +  \
                 h_m[0] * r2inv * (-6. * nv1 * nv2_2 + 4.5 * nv2_3 + nv2 * v1_2 -  \
                 4.* nv1 * v1v2 + 4. * nv2 * v1v2 + 4. * nv1 * v2_2 - 5. * nv2 *   \
                 v2_2));

        pn2d[k] = ndot[k] * ((-57./4. * mp2bh - 69./2. * mpbh2 - 9. * mbh3) * r4inv\
                  + h_m[0] * r2inv * (-15./8. * nv2_4 + 1.5 * nv2_2 * v1_2 - 6.*   \
                  nv2_2 * v1v2 - 2. * v1v2 * v1v2 + 4.5 * nv2_2 * v2_2 + 4. * v1v2 \
                  * v2_2 - 2. * v2_2 * v2_2) + mpbh * r3inv * (19.5 * nv1_2 - 39. *\
                  nv1 * nv2 + 8.5 * nv2 * nv2 - 15./4. * v1_2 - 2.5 * v1v2 + 5./4. \
                  * v2_2) + mbh2 * r3inv *(2. * nv1_2 - 4. * nv1 *nv2 - 6. * nv2_2 \
                  - 8. * v1v2 + 4. * v2_2)) + nvec[k] * (4. * rp * r5inv * (57./4. \
                  * mp2bh + 69./2. * mpbh2 + 9. * mbh3) - 2.* h_m[0] * r3inv *rp * \
                  (-15./8. * nv2_4 + 1.5 * nv2_2 * v1_2 - 6. * nv2_2 * v1v2 - 2. * \
                  v1v2 * v1v2 + 4.5 * nv2_2 * v2_2 + 4. * v1v2 * v2_2 - 2. * v2_2 *\
                  v2_2) + h_m[0] * r2inv * (-7.5 * nv2_3 * nv2dot + 1.5 * (2. *    \
                  nv2 * nv2dot * v1_2 + nv2_2 * v1a) - 6. * (2. * nv2 * nv2dot *   \
                  v1v2 + nv2_2 * v1v2dot - 4. * v1v2dot * v1v2 + 4.5 * (2. * nv2 * \
                  nv2dot * v2_2 + nv2_2 * v2a) +4. * (v1v2dot * v2_2 + v1v2 * v2a) \
                  - 4. * v2a)) -3. * mpbh * r4inv * rp * (19.5 * nv1 * nv1 - 39. * \
                  nv1 * nv2 + 8.5 * nv2 * nv2 - 15./4. * v1_2 - 2.5 * v1v2 + 5./4. \
                  * v2_2) + mpbh * r3inv * (39. * nv1 * nv1dot - 39. * (nv1dot *   \
                  nv2 + nv1 * nv2dot) + 17. * nv2 * nv2dot - 15./4. * v1a - 2.5 *  \
                  v1v2dot + 1.25 * v2a) - 3. * mbh2 * rp * r4inv * (2. * nv1_2 -   \
                  4. * nv1 * nv2 - 6. * nv2_2 - 8. * v1v2 + 4. * v2_2) + mbh2 *    \
                  r3inv * (4. * nv1 * nv1dot -4. * (nv1dot * nv2 + nv1 * nv2dot) - \
                  12. * nv2 * nv2dot - 8. * v1v2dot + 4. * v2a)) + (at1[k] -       \
                  at2[k]) * (mbh2 * r3inv * (-2. * nv1 - 2. * nv2) + mpbh * r3inv *\
                  (-63./4. * nv1 + 55./4. * nv2) + h_m[0] * r2inv * (-6. * nv1 *   \
                  nv2_2 + 4.5 * nv2_3 + nv2 * v1_2 - 4. * nv1 * v1v2 + 4. * nv2 *  \
                  v1v2 + 4. * nv1 * v2_2 - 5. * nv2 * v2_2)) + v[k] *(-3. * rp *   \
                  r4inv * ((-2. * nv1 - 2. * nv2) * mbh2 + mpbh * (-63./4. * nv1 + \
                  55./4. * nv2)) + 1. * r3inv * (mbh2 * (-2. * nv1dot - 2. * nv2dot\
                  ) + mpbh * (-63./4. * nv1dot + 55./4. * nv2dot)) - 2. * h_m[0] * \
                  r3inv * rp * (-6. * nv1 * nv2_2 + 4.5 * nv2_3 + nv2 * v1_2 - 4. *\
                  nv1 * v1v2 + 4. * nv2 * v1v2 + 4. * nv1 * v2_2 - 5. * nv2 * v2_2)\
                  + h_m[0] * r2inv * (-6. * (nv1dot * nv2_2 + 2. * nv1 * nv2 *     \
                  nv2dot) + 27./2. * nv2_2 * nv2dot + (nv2dot * v1_2 + nv2 * v1a) -\
                  4. * (nv1dot * v1v2 + nv1 * v1v2dot) + 4. * (nv2dot * v1v2 + nv2 \
                  * v1v2dot) + 4. * (nv1dot * v2_2 + nv1 * v2a) - 5. * (nv2dot *   \
                  v2_2 + nv2 * v2a)));

        pn25[k] = ((208. * mpbh2 * nvdot/15.-24./5. * mp2bh * nvdot) * r4inv + 12. \
                  * mpbh/5. * r3inv * nvdot * v1_v22) * nvec[k] + (8./5. * mp2bh * \
                  r4inv - 32./5. * mpbh2 * r4inv - 4./5. * mpbh * r3inv * v1_v22) *\
                  v[k];

        pn25d[k] = ndot[k] * ((208./15. * mpbh2 - 24./5. * mp2bh) * nvdot * r4inv +\
                   12./5. * mpbh * r3inv * nvdot * v1_v22) + nvec[k] * (-4. * r5inv\
                   * rp * nvdot * (208./15. * mpbh2 - 24./5. * mp2bh) + nvdot *    \
                   r4inv * (208./15. * mpbh2 - 24./5. * mp2bh) - 36./5. * mpbh *   \
                   r4inv * rp * nvdot * v1_v22 + 12./5. * mpbh * r3inv * (nvdot *  \
                   v1_v22 + nvdot * v12a)) + (at1[k] - at2[k]) * (mpbh * r4inv *   \
                   (8./5. * h_m[i] - 32./5. * h_m[0]) - 4./5. * mpbh * r3inv *     \
                   v1_v22) + v[k] * (-4. * mpbh * r5inv * rp * (8./5. * h_m[i] -   \
                   32./5. * h_m[0]) + 12./5. * mpbh * r4inv * rp * v1_v22 - 4./5. *\
                   mpbh * r3inv * v12a);
    }

    for (int k=0;k<3; k++)
    {
        h_f[i].a[k]  += pn1[k]/c2;
        h_f[i].a1[k] += pn1d[k]/c2;
    }
    printf("1PNa:  %.10f %.10f %.10f\n", pn1[0]/c2, pn1[1]/c2, pn1[2]/c2);
    printf("1PNb:  %.10f %.10f %.10f\n", pn1d[0]/c2, pn1d[1]/c2, pn1d[2]/c2);

    printf("2PNa:  %.10f %.10f %.10f\n", pn2[0]/c4, pn2[1]/c4, pn2[2]/c4);
    printf("2PNb:  %.10f %.10f %.10f\n", pn2d[0]/c4, pn2d[1]/c4, pn2d[2]/c4);

    printf("25PNa:  %.10f %.10f %.10f\n", pn25[0]/c5, pn25[1]/c5, pn25[2]/c5);
    printf("25PNb:  %.10f %.10f %.10f\n", pn25d[0]/c5, pn25d[1]/c5, pn25d[2]/c5);
    CHECK!
    getchar();
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
