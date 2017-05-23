/*
 * Copyright (c) 2016
 *
 * Cristián Maureira-Fredes <cmaureirafredes@gmail.com>
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. The name of the author may not be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#include "Hermite4CPU.hpp"


/** Constructor that uses its parent constructor */
Hermite4CPU::Hermite4CPU(NbodySystem *ns, Logger *logger, NbodyUtils *nu)
                         : Hermite4(ns, logger, nu)
{
    /* Empty */
}


/** Destructor */
Hermite4CPU::~Hermite4CPU()
{
    /* Empty */
}

/** Method that calculate the gravitational interaction between two particles
 */
void Hermite4CPU::force_calculation(const Predictor &pi, const Predictor &pj, Forces &fi)
{
    double rx = pj.r[0] - pi.r[0];
    double ry = pj.r[1] - pi.r[1];
    double rz = pj.r[2] - pi.r[2];

    double vx = pj.v[0] - pi.v[0];
    double vy = pj.v[1] - pi.v[1];
    double vz = pj.v[2] - pi.v[2];

    double r2     = rx*rx + ry*ry + rz*rz + ns->e2;
    double rinv   = 1.0/sqrt(r2);
    double r2inv  = rinv  * rinv;
    double r3inv  = r2inv * rinv;
    double r5inv  = r2inv * r3inv;
    double mr3inv = r3inv * pj.m;
    double mr5inv = r5inv * pj.m;

    double rv = rx*vx + ry*vy + rz*vz;

    fi.a[0] += (rx * mr3inv);
    fi.a[1] += (ry * mr3inv);
    fi.a[2] += (rz * mr3inv);

    fi.a1[0] += (vx * mr3inv - (3 * rv ) * rx * mr5inv);
    fi.a1[1] += (vy * mr3inv - (3 * rv ) * ry * mr5inv);
    fi.a1[2] += (vz * mr3inv - (3 * rv ) * rz * mr5inv);
}

#ifdef PN
// Adding PN1, PN2 and PN2 terms to the acceleration
// including its derivatives
void Hermite4CPU::force_calculation_pn(const Predictor &pi, const Predictor &pj, Forces &fi, Forces &fj, int id)
{
    double rx = pi.r[0] - pj.r[0];
    double ry = pi.r[1] - pj.r[1];
    double rz = pi.r[2] - pj.r[2];

    double vx = pi.v[0] - pj.v[0];
    double vy = pi.v[1] - pj.v[1];
    double vz = pi.v[2] - pj.v[2];
    double vv[3] = {vx, vy, vz};
    double zvv[3] = {-vx, -vy, -vz};

    double r2    = rx * rx + ry * ry + rz * rz;// + ns->e2;
    double r     = sqrt(r2);
    double rinv  = 1.0/r;
    double r2inv = rinv  * rinv;
    double r3inv = r2inv * rinv;
    double r4inv = r2inv * r2inv;
    double r5inv = r2inv * r3inv;

    double m1 = pi.m;
    double m2 = pj.m;
    double nv[3] = {rx * rinv, ry * rinv, rz * rinv};
    double G2 = G * G;
    double G3 = G2 * G;

    double m1_2 = m1 * m1;
    double m2_2 = m2 * m2;
    double m1m2 = m1 * m2;

    // Relative velocity
    double v1_2 = pi.v[0] * pi.v[0] + pi.v[1] * pi.v[1] + pi.v[2] * pi.v[2];
    double v2_2 = pj.v[0] * pj.v[0] + pj.v[1] * pj.v[1] + pj.v[2] * pj.v[2];
    double v1v2 = pi.v[0] * pj.v[0] + pi.v[1] * pj.v[1] + pi.v[2] * pj.v[2];
    double v2v1 = pj.v[0] * pi.v[0] + pj.v[1] * pi.v[1] + pj.v[2] * pi.v[2];
    double v2v1_2 = v2v1 * v2v1;

    double nv1 = nv[0] * pi.v[0] + nv[1] * pi.v[1] + nv[2] * pi.v[2];
    double nv2 = nv[0] * pj.v[0] + nv[1] * pj.v[1] + nv[2] * pj.v[2];
    double na1 = nv[0] * fi.a[0] + nv[1] * fi.a[1] + nv[2] * fi.a[2];
    double na2 = nv[0] * fj.a[0] + nv[1] * fj.a[1] + nv[2] * fj.a[2];

    double v1a1 = pi.v[0] * fi.a[0] + pi.v[1] * fi.a[1] + pi.v[2] * fi.a[2];
    double v1a2 = pi.v[0] * fj.a[0] + pi.v[1] * fj.a[1] + pi.v[2] * fj.a[2];
    double v2a1 = pj.v[0] * fi.a[0] + pj.v[1] * fi.a[1] + pj.v[2] * fi.a[2];
    double v2a2 = pj.v[0] * fj.a[0] + pj.v[1] * fj.a[1] + pj.v[2] * fj.a[2];
    double a1v2 = fi.a[0] * pj.v[0] + fi.a[1] * pj.v[1] + fi.a[2] * pj.v[2];
    double a2v2 = fj.a[0] * pj.v[0] + fj.a[1] * pj.v[1] + fj.a[2] * pj.v[2];
    double a2v1 = fj.a[0] * pi.v[0] + fj.a[1] * pi.v[1] + fj.a[2] * pi.v[2];

    double aa[3] = {fi.a[0] - fj.a[0], fi.a[1] - fj.a[1], fi.a[2] - fj.a[2]};
    double zaa[3] = {-aa[0], -aa[1], -aa[2]};

    double v2_4 = v2_2 * v2_2;      // V2^4
    double v1v2_2 = v1v2 * v1v2;    // (V1 * V2)^2
    double nv1_2 = nv1 * nv1;       // (N * V1)^2
    double nv2_2 = nv2 * nv2;       // (N * V2)^2
    double nv2_3 = nv2_2 * nv2;     // (N * V2)^3
    double nv2_4 = nv2_2 * nv2_2;   // (N * V2)^4

    double pn1[3] = {0};
    double pn2[3] = {0};
    double pn25[3] = {0};
    double dpn1[3] = {0};
    double dpn2[3] = {0};
    double dpn25[3] = {0};

    double c = SPEED_OF_LIGHT;
    double c1inv = 1.0 / c;
    double c2inv = c1inv * c1inv;
    double c4inv = c2inv * c2inv;
    double c5inv = c4inv * c1inv;

    for (int k = 0; k < 3; k++)
    {
    // Precession (PN1)
    pn1[k] = G * m2 * r2inv * (nv[k] * (-v1_2 - 2 * v2_2 + 4 * v1v2 + 1.5 * (nv2_2) +
             5 * (G * m1 * rinv) + 4 * (G * m2 * rinv)) + vv[k] * (4*nv1 - 3*nv2 ));
    pn1[k] *= c2inv;

    // Precession Correction (PN2)
    pn2[k] = G * m2 * r2inv * (nv[k] * (-2 * v2_4 + 4 * v2_2 * v1v2 - 2 * v1v2_2 +
             1.5 * v1_2 * nv2_2 + 4.5 * v2_2 * nv2_2 - 6 * v1v2 * nv2_2 - 1.875 *
             nv2_4 + (G * m1 * rinv) * (-3.75 * v1_2 + 1.25 * v2_2 - 2.5 *
             v1v2 + 19.5 * nv1_2 - 39 * nv1 * nv2 + 8.5 * nv2_2) + G * m2 *
             rinv * (4 * v2_2 - 8 * v1v2 + 2 * nv1_2 - 4 * nv1 * nv2 - 6 * nv2_2)) +
             vv[k] * (v1_2 * nv2 + 4 * v2_2 * nv1 - 5 * v2_2 * nv2 - 4 * v1v2 * nv1 +
             4 * v1v2 * nv2 - 6 * nv1 * nv2_2 + 4.5 * nv2_3 + (G * m1 * rinv) *
             (-15.75 * nv1 + 13.75 * nv2) + G * m2 * rinv * (-2 * nv1 - 2 *
             nv2))) + G3 * m2 * r4inv * nv[k] * (-14.25 * m1_2 - 9 * m2_2 -
             34.5 * m1 * m2);
    pn2[k] *= c4inv;

    // GW Decay (PN2.5)
    pn25[k] = 0.8 * G2 * m1 * m2 * r3inv * (vv[k] * (-(vv[k] * vv[k]) + 2 * G *
              m1 * rinv - 8 * G * m2 * rinv) + nv[k] * (nv1 - nv2) * (3 * vv[k] *
              vv[k] - 6 * G * m1 * rinv + 52.0/3.0 * G * m2 * rinv));
    pn25[k] *= c5inv;

    // Derivative Precession (dPN1)
    dpn1[k] = G * m2 * (-(vv[k] * v1_2 * r3inv + 2 * nv[k] * v1a1 * r2inv + 3 *
              nv[k] * v1_2 * (nv2 - nv1) * r3inv) -2 * (vv[k] * v2_2 * r3inv + 2 *
              nv[k] * v2a2 * r2inv + 3 * nv[k] * v2_2 * (nv2 - nv1) * r3inv) + 4 *
              (vv[k] * v1v2 * r3inv + nv[k] * (a1v2 + a2v1) * r2inv + 3 * nv[k] *
              v1v2 * (nv2 - nv1) * r3inv) + 1.5 * (vv[k] * nv2_2 * r3inv + 2 *
              nv[k] * nv2 * (r * na2 + v1v2 - v2_2) * r3inv) + 5 * nv[k] * (nv2_2 *
              (nv2 - nv1)) * r3inv + G * (vv[k] * r4inv + 4 * nv[k] * (nv2 - nv1) *
              r4inv) * (5 * m1 + 4 * m2) + 4 * nv1 * r2inv * aa[k] + 3 * nv2 *
              r2inv * zaa[k] + 4 * (v1_2 - v1v2 + r * na1 + 3 *(nv2 - nv1) * nv1) *
              r3inv * vv[k] + 3 * (v1v2 - v2_2 + r * na2 + 3 * (nv2 - nv1) * nv2) *
              r3inv * zvv[k]);
    dpn1[k] *= c2inv;

    // Derivative Precession Correction (dPN2)
    dpn2[k] = G * m2 * (-2 * (vv[k] * v2_4 * r3inv + nv[k] * (4 * v2_2 * a2v2) *
              r2inv + 3 * nv[k] * (v2_4 * (nv2 - nv1)) * r3inv) + 4 *(vv[k] * v2_2 *
              v1v2 * r3inv + 2 * nv[k] * (v2a2 * v1v2) * r2inv + nv[k] * (v2_2 *
              (a1v2 + v1a2)) * r2inv * 3 * nv[k] * (v2_2 * v1v2 * (nv2 - nv1)) *
              r3inv) - 2 * (vv[k] * v2v1_2 * r3inv + 2 * nv[k] * (v1v2 * (a1v2 +
              a2v1)) * r2inv + 3 * nv[k] * (v1v2_2 * (nv2 - nv1)) * r3inv) + 1.5 *
              (vv[k] * (v1_2 * nv2_2) * r3inv + 2 * nv[k] * (v1a1 * nv2_2) * r2inv +
              2 * nv[k] * (v1_2 * nv2) * r2inv * (na2 + (v1v2 - v2_2) * rinv) + 5 *
              nv[k] * (v1_2 * nv2_2) * r3inv * (nv2 - nv1)) + 4.5 * (vv[k] * (v2_2 *
              nv2_2) * r3inv + 2 * nv[k] * (v2a2 * nv2_2) * r2inv + 2 * nv[k] *
              (v2_2 * nv2) * r2inv * (na2 + (v1v2 - v2_2) * rinv) + 5 * nv[k] *
              (v2_2 * nv2_2) * r3inv * (nv2 - nv1)) -6 * (vv[k] * (v1v2 * nv2_2) *
              r3inv + nv[k] * ((a1v2 + v1a2) * nv2_2) * r2inv + 2 * nv[k] * (v1v2 *
              nv2) * r2inv * (na2 + (v1v2 - v2_2) * rinv) + 5 * nv[k] * (v1v2 *
              nv2_2 * (nv2 - nv1)) * r3inv) - 1.875 * (vv[k] * nv2_4 * r3inv +
              4 * nv[k] * nv2_3 * r2inv * (na2 + (v1v2 - v2_2) * rinv) + 7 * nv[k] *
              nv2_4 * r3inv * (nv2 - nv1)) + G * m1 * (-3.75 * (vv[k] * v1_2 *
              r4inv + 2 * nv[k] * v1a1 * r3inv + 4 * nv[k] * (v1_2 * (nv2 - nv1)) *
              r4inv) + 1.25 * (vv[k] * v2_2 * r4inv + 2 * nv[k] * v2a2 * r3inv +
              4 * nv[k] * v2_2 * (nv2 - nv1) *r4inv) - 2.5 * (vv[k] * v1v2 *
              r4inv + nv[k] * (a1v2 + v1a2) * r3inv + 4 * nv[k] * v1v2 * (nv2 -
              nv1) * r4inv) + 19.5 * (vv[k] * nv1_2 * r4inv + 2 * nv[k] * nv1 *
              r3inv * (na1 + (v1_2 - v1v2) * rinv) + 6 * nv[k] * (nv1_2 * (nv2 -
              nv1)) * r4inv) - 39.0 * (vv[k] * nv1 * nv2 * r4inv + nv[k] * r3inv *
              (nv1 * na2 + nv2 * na1 + (nv1 * (v1v2 - v2_2) * rinv) + (nv2 * (v1_2 -
              v1v2) * rinv)) + 6 * nv[k] * nv1 * nv2 * r4inv * (nv2 - nv1)) +
              8.5 * (vv[k] * nv2_2 * r4inv + 2 * nv[k] * nv2 * r3inv * (na2 +
              (v1v2 - v2_2) * rinv) + 6 * nv[k] * (nv2_2 * (nv2 - nv1)) * r4inv)) +
              G * m2 * (4.0 * (vv[k] * v2_2 * r4inv + 2 * nv[k] * v2a2 * r3inv + 4 *
              nv[k] * (v2_2 * (nv2 - nv1)) * r4inv) - 8.0 * (vv[k] * v1v2 * r4inv +
              nv[k] * (a1v2 + v1a2) * r3inv + 4 * nv[k] * (v1v2 * (nv2 - nv1)) *
              r4inv) + 2.0 * (vv[k] * nv1_2 * r4inv + 2 * nv[k] * nv1 * r3inv *
              (na1 + (v1_2 * v1v2) * rinv) + 6 * nv[k] * (nv1_2 * (nv2 - nv1)) *
              r4inv) - 4.0 * (vv[k] * nv1 * nv2 * r4inv + nv[k] * r3inv * (nv1 *
              na2 + nv2 * na1 + (nv1 * (v1v2 - v2_2)) * rinv + (nv2 * (v1_2 *
              v1v2)) * rinv) + 6 * nv[k] * nv1 * nv2 * r4inv * (nv2 - nv1)) - 6.0 *
              (vv[k] * nv2_2 * r4inv + 2 * nv[k] * nv2 * r3inv * (na2 + (v1v2 -
              v2_2) * rinv) + 6 * nv[k] * (nv2_2 * (nv2 - nv1)) * r4inv)) + aa[k] *
              v1_2 * nv2 * r2inv + vv[k] * (2 * v1a2 * nv2 * r2inv + v1_2 * r2inv *
              (na2 + (v1v2 - v2_2) * rinv + 3 * nv2 * rinv * (nv2 - nv1))) + 4.0 *
              aa[k] * v2_2 * nv1 * r2inv + 4 * vv[k] * (2 * v2a2 * nv1 * r2inv +
              v2_2 * r2inv * (na1 + (v1_2 - v1v2) * rinv + 3 * nv1 * rinv * (nv2 -
              nv1))) - 5.0 * aa[k] * v2_2 * nv2 * r2inv - 5 * vv[k] * (2 * v2a2 *
              nv2 * r2inv + v2_2 * r2inv * (na2 + (v1v2 - v2_2) * rinv + 3 * nv2 *
              rinv * (nv2 - nv1))) - 4.0 * aa[k] * v1v2 * nv1 * r2inv - 4 * vv[k] *
              ((a1v2 + v1a2) * nv1 * r2inv + v1v2 * r2inv * (na1 + (v1_2 - v1v2) *
              rinv + 3 * v1v2 * nv1 * rinv * (nv2 - nv1))) + 4.0 * aa[k] * v1v2 *
              nv2 * r2inv + 4 * vv[k] * ((a1v2 + v1a2) * nv2 * r2inv + v1v2 * r2inv *
              (na2 + (v1v2 - v2_2) * rinv + 3 * v1v2 * nv2 * rinv * (nv2 - nv1))) -
              6 * aa[k] * nv1 * nv2_2 * r2inv - 6 * vv[k] * (nv2_2 * r2inv * (na1 +
              (v1_2 - v1v2) * rinv) + 2 * nv1 * nv2 * r2inv * (na2 + (v1v2 - v2_2) *
              rinv) + 5 * nv1 * nv2_2 * r3inv * (nv2 - nv1)) + 4.5 * aa[k] * nv2_3 *
              r2inv + 4.5 * vv[k] * (3 * nv2_2 * r2inv * (na2 + (v1v2 - v2_2) *
              rinv) + 5 * r * nv2_3 * r3inv * (nv2 - nv1)) + G * ((aa[k] * nv1 *
              r3inv * vv[k] * (na1 * r3inv + (v1_2 - v1v2) * r4inv + 4 * nv1 *
              r4inv *(nv2 - nv1))) * (-15.75 * m1 - 2 * m2) + (aa[k] * nv2 *
              r3inv + vv[k] * (na2 * r3inv + (v1v2 - v2_2) * r4inv + 4 * nv2 *
              r4inv * (nv2 - nv1))) * (13.75 * m1 - 2 * m2))) + G3 * m2 *
              (-14.25 * m1_2 - 9 * m2_2 - 34.5 * m1m2) * (vv[k] * r5inv +
              5 * nv[k] * (nv2 - nv1) * r5inv);
    dpn2[k] *= c4inv;

    // Derivative GW Decay (dPN2.5)
    dpn25[k] = 0.8 * G2 * m1m2 * (-(aa[k] * vv[k] * vv[k]) * r3inv - 2 * vv[k] *
               r3inv * (v1a1 + v2a2 - v2a1 - v1a2) + 6 * vv[k] * r4inv * vv[k] *
               vv[k] * (nv1 - nv2) + G * (2 * m1 - 8 * m2) * (aa[k] * r4inv + 4 *
               vv[k] * r5inv * (nv2 - nv1)) + 3 * (nv[k] * vv[k] * vv[k] * r3inv *
               (na1 - na2 + (v1_2 + v2_2 - 2 * v1v2) * rinv) + 2 * nv[k] * (nv1 -
               nv2) * r3inv * (v1a2 + v2a2 - v2a1 - v1a2) - 5 * nv[k] * (vv[k] *
               vv[k] * (nv1 - nv2) * (nv1 - nv2)) * r4inv) + G * (52.0/3.0 * m2 -
               6 * m1) * (vv[k] * (nv1-nv2) * r5inv + nv[k] * r4inv * (na1 - na2 +
               (v1_2 - 2 * v1v2 + v2_2) * rinv) - 6 *nv[k] * (nv2 - nv1) * (nv2 -
               nv1) * r5inv));
    dpn25[k] *= c5inv;
    }

    fi.a[0] += pn1[0] + pn2[0] + pn25[0];
    fi.a[1] += pn1[1] + pn2[1] + pn25[1];
    fi.a[2] += pn1[2] + pn2[2] + pn25[2];

    fi.a1[0] += dpn1[0] + dpn2[0] + dpn25[0];
    fi.a1[1] += dpn1[1] + dpn2[1] + dpn25[1];
    fi.a1[2] += dpn1[2] + dpn2[2] + dpn25[2];
}
#endif

/** Method that initializes the acceleration and it first derivative
 */
void Hermite4CPU::init_acc_jrk()
{
    unsigned int j = 0;
    #pragma omp parallel for private(j)
    for (unsigned int i = 0; i < ns->n; i++)
    {
        for (j = 0; j < ns->n; j++)
        {
            if(i == j) continue;
            force_calculation(ns->h_p[i], ns->h_p[j], ns->h_f[i]);
            #ifdef PN
            force_calculation_pn(ns->h_p[i], ns->h_p[j], ns->h_f[i], ns->h_f[j], i);
            #endif
        }
    }
}

/** Method that call the force_calculation method for every \f$i-\f$ and \f$j\f$
 * particles interaction of the \f$N_{act}\f$ ones.
 */
void Hermite4CPU::update_acc_jrk(unsigned int nact)
{
    ns->gtime.update_ini = omp_get_wtime();
    unsigned int i = 0;
    unsigned int j = 0;
    #pragma omp parallel for private(i,j)
    for (unsigned int k = 0; k < nact; k++)
    {
        i = ns->h_move[k];
        ns->h_f[i].a[0]  = 0.0;
        ns->h_f[i].a[1]  = 0.0;
        ns->h_f[i].a[2]  = 0.0;
        ns->h_f[i].a1[0] = 0.0;
        ns->h_f[i].a1[1] = 0.0;
        ns->h_f[i].a1[2] = 0.0;

        #pragma omp parallel for
        for (j = 0; j < ns->n; j++)
        {
            if(i == j) continue;
            force_calculation(ns->h_p[i], ns->h_p[j], ns->h_f[i]);
            #ifdef PN
            force_calculation_pn(ns->h_p[i], ns->h_p[j], ns->h_f[i], ns->h_f[j], i);
            #endif
        }
    }
    ns->gtime.update_end += omp_get_wtime() - ns->gtime.update_ini;
}

/** Method that predict all the particles to the current integration time
 */
void Hermite4CPU::predicted_pos_vel(double ITIME)
{

    ns->gtime.prediction_ini = omp_get_wtime();
    for (unsigned int i = 0; i < ns->n; i++)
    {
        double dt  = ITIME - ns->h_t[i];
        double dt2 = (dt  * dt);
        double dt3 = (dt2 * dt);

        double dt2c = dt2/2.0;
        double dt3c = dt3/6.0;

        ns->h_p[i].r[0] = (dt3c * ns->h_f[i].a1[0]) + (dt2c * ns->h_f[i].a[0]) + (dt * ns->h_v[i].x) + ns->h_r[i].x;
        ns->h_p[i].r[1] = (dt3c * ns->h_f[i].a1[1]) + (dt2c * ns->h_f[i].a[1]) + (dt * ns->h_v[i].y) + ns->h_r[i].y;
        ns->h_p[i].r[2] = (dt3c * ns->h_f[i].a1[2]) + (dt2c * ns->h_f[i].a[2]) + (dt * ns->h_v[i].z) + ns->h_r[i].z;

        ns->h_p[i].v[0] = (dt2c * ns->h_f[i].a1[0]) + (dt * ns->h_f[i].a[0]) + ns->h_v[i].x;
        ns->h_p[i].v[1] = (dt2c * ns->h_f[i].a1[1]) + (dt * ns->h_f[i].a[1]) + ns->h_v[i].y;
        ns->h_p[i].v[2] = (dt2c * ns->h_f[i].a1[2]) + (dt * ns->h_f[i].a[2]) + ns->h_v[i].z;

        ns->h_p[i].m = ns->h_r[i].w;

    }
    ns->gtime.prediction_end += omp_get_wtime() - ns->gtime.prediction_ini;
}

/** Method that correct the positions and velocities of the particles at the
 * end of every integration step
 */
void Hermite4CPU::correction_pos_vel(double ITIME, unsigned int nact)
{
    ns->gtime.correction_ini = omp_get_wtime();
    for (unsigned int k = 0; k < nact; k++)
    {
        int i = ns->h_move[k];

        double dt1 = ns->h_dt[i];
        double dt2 = dt1 * dt1;
        double dt3 = dt2 * dt1;
        double dt4 = dt2 * dt2;
        double dt5 = dt4 * dt1;

        // Acceleration 2nd derivate
        ns->h_a2[i].x = (-6 * (ns->h_old[i].a[0] - ns->h_f[i].a[0] ) - dt1 * (4 * ns->h_old[i].a1[0] + 2 * ns->h_f[i].a1[0]) ) / dt2;
        ns->h_a2[i].y = (-6 * (ns->h_old[i].a[1] - ns->h_f[i].a[1] ) - dt1 * (4 * ns->h_old[i].a1[1] + 2 * ns->h_f[i].a1[1]) ) / dt2;
        ns->h_a2[i].z = (-6 * (ns->h_old[i].a[2] - ns->h_f[i].a[2] ) - dt1 * (4 * ns->h_old[i].a1[2] + 2 * ns->h_f[i].a1[2]) ) / dt2;

        // Acceleration 3rd derivate
        ns->h_a3[i].x = (12 * (ns->h_old[i].a[0] - ns->h_f[i].a[0] ) + 6 * dt1 * (ns->h_old[i].a1[0] + ns->h_f[i].a1[0]) ) / dt3;
        ns->h_a3[i].y = (12 * (ns->h_old[i].a[1] - ns->h_f[i].a[1] ) + 6 * dt1 * (ns->h_old[i].a1[1] + ns->h_f[i].a1[1]) ) / dt3;
        ns->h_a3[i].z = (12 * (ns->h_old[i].a[2] - ns->h_f[i].a[2] ) + 6 * dt1 * (ns->h_old[i].a1[2] + ns->h_f[i].a1[2]) ) / dt3;

        // Correcting position
        ns->h_r[i].x = ns->h_p[i].r[0] + (dt4/24)*ns->h_a2[i].x + (dt5/120)*ns->h_a3[i].x;
        ns->h_r[i].y = ns->h_p[i].r[1] + (dt4/24)*ns->h_a2[i].y + (dt5/120)*ns->h_a3[i].y;
        ns->h_r[i].z = ns->h_p[i].r[2] + (dt4/24)*ns->h_a2[i].z + (dt5/120)*ns->h_a3[i].z;

        // Correcting velocity
        ns->h_v[i].x = ns->h_p[i].v[0] + (dt3/6)*ns->h_a2[i].x + (dt4/24)*ns->h_a3[i].x;
        ns->h_v[i].y = ns->h_p[i].v[1] + (dt3/6)*ns->h_a2[i].y + (dt4/24)*ns->h_a3[i].y;
        ns->h_v[i].z = ns->h_p[i].v[2] + (dt3/6)*ns->h_a2[i].z + (dt4/24)*ns->h_a3[i].z;


        ns->h_t[i] = ITIME;
        double normal_dt  = nu->get_timestep_normal(i, ns->eta);
        normal_dt = nu->normalize_dt(normal_dt, ns->h_dt[i], ns->h_t[i], i);
        ns->h_dt[i] = normal_dt;

    }
    ns->gtime.correction_end += omp_get_wtime() - ns->gtime.correction_ini;
}
