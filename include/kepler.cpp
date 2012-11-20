#include "kepler.hpp"

static double
    _1_17_16 = 1./17./16.,
    _1_16_15 = 1./16./15.,
    _1_15_14 = 1./15./14.,
    _1_14_13 = 1./14./13.,
    _1_13_12 = 1./13./12.,
    _1_12_11 = 1./12./11.,
    _1_11_10 = 1./11./10.,
    _1_10_9 = 1./10./9.,
    _1_9_8 = 1./9./8.,
    _1_8_7 = 1./8./7.,
    _1_7_6 = 1./7./6.,
    _1_6_5 = 1./6./5.,
    _1_4_3 = 1./4./3.,
    _1_3_2 = 1./3./2.;

void kepler_prediction(double *rx, double *ry, double *rz,
                       double *vx, double *vy, double *vz,
                       double dt, int i)
{

    /*
     * Orbital parameters
     */
    double jx, jy, jz; // Angular momentum vector
    double ex, ey, ez; // Runge-Lenz vector
    double  a,  b,  w; // Semi-major axis, Semi-minor axis, Frequency of the orbit
    double ax, ay, az; // Semi-major axis vector
    double bx, by, bz; // Semi-minor axis vector
    double ecc;        // Eccentricity
    double e_anomaly, m_anomaly; // Eccentric anomaly, Mean anomaly
    double r_const, v_const;
    double r = sqrt((*rx) * (*rx) + (*ry) * (*ry) + (*rz) * (*rz));
    double cos_e, sin_e, e_tmp, ecc_const, cos_const;
    double m0 = h_m[0];

    #ifdef DEBUG_KEPLER
    printf("[M_bh] %.15f\n",m0);
    #endif

    #ifdef DEBUG_KEPLER
    printf("[relative] Old Position: (%.10f,%.10f,%.10f)\n", *rx, *ry, *rz);
    printf("[relative] Old Velocity: (%.10f,%.10f,%.10f)\n", *vx, *vy, *vz);
    #endif

    // Angular momentum
    // j = r x v
    jx = (*ry) * (*vz) - (*rz) * (*vy);
    jy = (*rz) * (*vx) - (*rx) * (*vz);
    jz = (*rx) * (*vy) - (*ry) * (*vx);
    #ifdef DEBUG_KEPLER
    printf("Angular momentum vector (j): %.10f %.10f %.10f\n", jx, jy, jz);
    #endif

    // Runge-Lenz-vector
    // e = { (v x j) / (G * m) }  - { r / |r| }
    ex = (*vy) * (jz) - (*vz) * (jy);
    ey = (*vz) * (jx) - (*vx) * (jz);
    ez = (*vx) * (jy) - (*vy) * (jx);

    // G, omitted due its value of 1
    ex = ex/m0 - *rx/r;
    ey = ey/m0 - *ry/r;
    ez = ez/m0 - *rz/r;
    #ifdef DEBUG_KEPLER
    printf("Runge-Lenz-vector (e): %.10f %.10f %.10f\n", ex, ey, ez);
    #endif

    // Eccentricity
    ecc = sqrt(ex*ex + ey*ey + ez*ez);
    #ifdef DEBUG_KEPLER
    printf("Eccentricity (ecc): %.10f\n", ecc);
    #endif

    // Semi-major axis
    // a = ( j * j ) / (G * m * | 1 - ecc^2 | )
    a = (jx*jx + jy*jy + jz*jz) / (G * m0 * fabs(1 - ecc*ecc));
    #ifdef DEBUG_KEPLER
    printf("Semi-major axis (a): %.10f\n", a);
    #endif

    // Frequency of the orbit
    // w = sqrt(( G * m )/ a^3 )
    w = sqrt((G * m0) / (a*a*a));
    #ifdef DEBUG_KEPLER
    printf("Frequency of the orbit (w): %.10f\n", w);
    #endif

    // Semi-major vector
    ax = a * ex/ecc;
    ay = a * ey/ecc;
    az = a * ez/ecc;
    #ifdef DEBUG_KEPLER
    printf("Semi-major vector: %.10f %.10f %.10f\n", ax, ay, az);
    #endif

    // Semi-minor vector
    bx = (jy * ez - jz * ey);
    by = (jz * ex - jx * ez);
    bz = (jx * ey - jy * ex);

    double b_vmag = sqrt(bx*bx + by*by + bz*bz);
    // Semi-minor axis
    b = a * sqrt(fabs(1 - ecc*ecc));
    //b = a * sqrt(fabs(1 - ecc*ecc))/b_vmag;
    #ifdef DEBUG_KEPLER
    printf("Semi-minor axis (b): %.10f\n", b);
    #endif

    bx *= b/b_vmag;
    by *= b/b_vmag;
    bz *= b/b_vmag;
    #ifdef DEBUG_KEPLER
    printf("Semi-minor vector:  %.10f %.10f %.10f\n", bx, by, bz);
    #endif
    // End of the variables calculation.


    if(ecc < 1)
    {
        /*
         * Elliptical Orbit case
         */

        // Calculate Eccentric Anomaly
        e_anomaly = (a - r) / (ecc * a);

        if(e_anomaly >= 1.0)
            e_anomaly = 0.0;
        else if(e_anomaly <= -1.0)
            e_anomaly = M_PI;
        else
            e_anomaly = acos(e_anomaly);

        // TO DO: Why?
        if (sqrt(*rx*bx + *ry*by + *rz*bz) < 0)
        {
            e_anomaly = 2*M_PI - e_anomaly;
        }
        // Calculate Mean Anomaly
        m_anomaly = (e_anomaly - ecc*sin(e_anomaly)) + dt * w;

        // Adjusting M anomaly to be < 2 * pi
        //while(m_anomaly >= 2*M_PI)
        //    m_anomaly -= 2 * M_PI;
        if(m_anomaly >= 2 * M_PI)
            m_anomaly = fmod(m_anomaly, 2 * M_PI);

        // Solving the Kepler equation for elliptical orbits
        e_anomaly = solve_kepler(m_anomaly, ecc);


        #ifdef DEBUG_KEPLER
        printf("[Elliptic orbit] Mean anomaly: %.10f\n", m_anomaly);
        printf("[Elliptic orbit] Ecce anomaly: %.10f\n", e_anomaly);
        #endif
        cos_e = cos(e_anomaly);
        sin_e = sin(e_anomaly);

        r_const = cos_e - ecc;
        v_const = w / (1.0 - ecc * cos_e);

        #ifdef DEBUG_KEPLER
        printf("[Elliptic orbit] r_const: %.10f\n", r_const);
        printf("[Elliptic orbit] v_const: %.10f\n", v_const);
        #endif

        // Based on Loeckmann and Baumgardt simulations criteria
        // This work better with 0.99 < e < 1 and |E| < 1e-3
        if(ecc > 0.99)
        {
            #ifdef DEBUG_KEPLER
            printf("[Elliptic orbit] Pericentre passage. e= %.10f\n", ecc);
            #endif
            e_tmp = (e_anomaly > 2.0 * M_PI - 1e-3) ? e_anomaly - 2.0 * M_PI : e_anomaly;
            if(e_tmp < 1e-3)
            {
                printf("[Elliptic orbit] |E| < 1e-3\n");
                e_tmp *= e_tmp;
                ecc_const = (jx*jx + jy*jy + jz*jz)/(m0 * a * (1 + ecc));
                // falta 1 término
                cos_const = -0.5 * e_tmp * (1 - e_tmp / 12.0 * (1 - e_tmp / 30.0));

                r_const = ecc_const + cos_const;
                v_const = w / (ecc_const - ecc * cos_const);
                #ifdef DEBUG_KEPLER
                printf("[Elliptic orbit] r_const: %.10f\n", r_const);
                printf("[Elliptic orbit] v_const: %.10f\n", v_const);
                #endif

            }
        }

        // New position
        *rx =   ax * r_const + bx * sin_e;
        *ry =   ay * r_const + by * sin_e;
        *rx =   ax * r_const + bx * sin_e;

        // New velocity
        *vx = (-ax * sin_e + bx * cos_e) * v_const;
        *vy = (-ay * sin_e + by * cos_e) * v_const;
        *vz = (-az * sin_e + bz * cos_e) * v_const;

        #ifdef DEBUG_KEPLER
        printf("[Elliptic orbit] New Position: (%.10f,%.10f,%.10f)\n", *rx, *ry, *rz);
        printf("[Elliptic orbit] New Velocity: (%.10f,%.10f,%.10f)\n", *vx, *vy, *vz);
        #endif

    }
    else
    {
        /*
         * Hyperbolic or Parabolic  Orbit case
         */

        // calculate eccentric anomaly e at t+dt
        e_anomaly = (a + r) / (ecc * a);
        if(e_anomaly < 1.0)
            e_anomaly = 0.0;
        else if((*rx * *vx + *ry * *vy + *rz * *vz) < 0)
            e_anomaly = -acosh(e_anomaly);
        else
            e_anomaly = acosh(e_anomaly);

        m_anomaly = (ecc * sinh(e_anomaly) - e_anomaly) + dt * w;

        //if(m_anomaly >= 2 * M_PI)
        //    m_anomaly = fmod(m_anomaly, 2 * M_PI);

        // Solving Kepler's equation for Hyperbolic/Parabolic orbits
        e_anomaly = kepler(ecc, m_anomaly);
        #ifdef DEBUG_KEPLER
        printf("[Hyp/Par orbit] Mean anomaly: %.10f\n", m_anomaly);
        printf("[Hyp/Par orbit] Ecce anomaly: %.10f\n", e_anomaly);
        #endif
        cos_e = cosh(e_anomaly);
        sin_e = sinh(e_anomaly);
        v_const = w / (ecc * cos_e - 1.);


        // New position
        *rx = ax * (ecc - cos_e)  + bx * sin_e;
        *ry = ay * (ecc - cos_e)  + by * sin_e;
        *rz = az * (ecc - cos_e)  + bz * sin_e;

        // New velocity
        *vx = (-ax * sin_e + bx * cos_e) * v_const;  // direction of v only
        *vy = (-ay * sin_e + by * cos_e) * v_const;  // direction of v only
        *vz = (-az * sin_e + bz * cos_e) * v_const;  // direction of v only

        #ifdef DEBUG_KEPLER
        printf("[Hyp/Par orbit] New Position: (%.10f,%.10f,%.10f)\n", *rx,*ry,*rz);
        printf("[Hyp/Par orbit] New Velocity: (%.10f,%.10f,%.10f)\n", *vx,*vy,*vz);
        #endif

    }
    // End of the orbit-type treatment

    double r2 = (*rx * *rx + *ry * *ry + *rz * *rz);
    double v2 = (*vx * *vx + *vy * *vy + *vz * *vz);

    //j = sqrt(jx*jx + jy*jy + jz*jz);
    //rv2 = (*rx * *vx + *ry * *vy + *rz * *vz);
    //v = j / (r * v * sin(acos(rv2/(r * v))));

    // total motion relative to fix central mass
    h_p_r[i].x = h_p_r[0].x + (*rx);
    h_p_r[i].y = h_p_r[0].y + (*ry);
    h_p_r[i].z = h_p_r[0].z + (*rz);

    h_p_v[i].x = h_p_v[0].x + (*vx);
    h_p_v[i].y = h_p_v[0].y + (*vy);
    h_p_v[i].z = h_p_v[0].z + (*vz);

    // Force contribution of the central mass on a particle
    // and its derivates


    double r2inv = 1/(r2);
    double mr3inv = -G * m0 * r2inv * sqrt(r2inv);

    // Acceleration (a)
    // a = -G * m0 * \vec{r} / r^{3}
    double a0x =  (*rx) * mr3inv;
    double a0y =  (*ry) * mr3inv;
    double a0z =  (*rz) * mr3inv;

    double vr = (*vx)*(*rx) + (*vy)*(*ry) + (*vz)*(*rz);

    double a1x = mr3inv * ( (*vx) - 3 * r2inv * vr * (*rx));
    double a1y = mr3inv * ( (*vy) - 3 * r2inv * vr * (*ry));
    double a1z = mr3inv * ( (*vz) - 3 * r2inv * vr * (*rz));


    double ra0 = ((*rx)*a0x + (*ry)*a0y + (*rz)*a0z);

    double a2x = mr3inv * (a0x - 3.0 * r2inv * (vr * ( 2.0 * (*vx) - 5.0 * vr * (*rx) * r2inv)));
    double a2y = mr3inv * (a0y - 3.0 * r2inv * (vr * ( 2.0 * (*vy) - 5.0 * vr * (*ry) * r2inv)));
    double a2z = mr3inv * (a0z - 3.0 * r2inv * (vr * ( 2.0 * (*vz) - 5.0 * vr * (*rz) * r2inv)));

    double va = ((*vx)*a0x + (*vy)*a0y + (*vz)*a0z);
    double ra1 = ((*rx)*a1x + (*ry)*a1y + (*rz)*a1z);

    double a3x = mr3inv * ( a1x - 3.0 * r2inv * ( 3.0 * vr * a0x + 3.0
                * (v2 + ra0) * ((*vx) - 5.0 * vr * r2inv * (*rx))
                + (3.0 * va + ra1) * (*rx) + (vr * vr * r2inv)
                * (-15.0 * (*vx) + 35.0 * vr * r2inv * (*rx))));

    double a3y = mr3inv * ( a1y - 3.0 * r2inv * ( 3.0 * vr * a0y + 3.0
                * (v2 + ra0) * ((*vy) - 5.0 * vr * r2inv * (*ry))
                + (3.0 * va + ra1) * (*ry) + (vr * vr * r2inv)
                * (-15.0 * (*vy) + 35.0 * vr * r2inv * (*ry))));

    double a3z = mr3inv * ( a1z - 3.0 * r2inv * ( 3.0 * vr * a0z + 3.0
                * (v2 + ra0) * ((*vz) - 5.0 * vr * r2inv * (*rz))
                + (3.0 * va + ra1) * (*rz) + (vr * vr * r2inv)
                * (-15.0 * (*vz) + 35.0 * vr * r2inv * (*rz))));


    h_a[i].x = a0x;
    h_a[i].y = a0y;
    h_a[i].z = a0z;

    h_a1[i].x = a1x;
    h_a1[i].y = a1y;
    h_a1[i].z = a1z;

    h_a2[i].x = a2x;
    h_a2[i].y = a2y;
    h_a2[i].z = a2z;

    h_a3[i].x = a3x;
    h_a3[i].y = a3y;
    h_a3[i].z = a3z;

    //printf("a0 = (%.10f, %.10f, %.10f)\n", a0x, a0y, a0z);
    //printf("a1 = (%.10f, %.10f, %.10f)\n", a1x, a1y, a1z);
    //printf("a2 = (%.10f, %.10f, %.10f)\n", a2x, a2y, a2z);
    //printf("a3 = (%.10f, %.10f, %.10f)\n", a3x, a3y, a3z);
}

double solve_kepler(double m_anomaly, double ecc)
{
    // Solving Kepler's equation
    double e_new = ecc > 0.8 ? M_PI : m_anomaly;
    double d = 0;
    int e_iter = 0;

    while(fabs(d) > DEL_E)
    {
        d = e_new - ecc * sin(e_new) - m_anomaly;
        if(e_iter-1 >= KEPLER_ITE)
            break;

        e_new -= d / (1.0 - ecc * cos(e_new));
        e_iter++;
    }
    return e_new;
}

/*
 * Following function taken from
 * http://www.projectpluto.com/kepler.htm
 */
double kepler(const double ecc, double mean_anom)
{
    #ifdef DEBUG_KEPLER
    printf("[DEBUG] kepler()\n");
    #endif
    double curr, err;
    int is_negative = 0, n_iter = 0;

    if(!mean_anom)
        return(0.);

   // if(ecc < .3)     /* low-eccentricity formula from Meeus,  p. 195 */
   // {
   //    curr = atan2(sin(mean_anom), cos(mean_anom) - ecc);
   //     /* one correction step,  and we're done */
   //     err = curr - ecc * sin(curr) - mean_anom;
   //     curr -= err / (1. - ecc * cos(curr));
   //     //printf("(1) %.10f\n", curr);
   //     return(curr);
   // }

    if(mean_anom < 0.)
    {
        mean_anom = -mean_anom;
        is_negative = 1;
    }

    curr = mean_anom;
    if((ecc > .8 && mean_anom < M_PI / 3.) || ecc > 1.)    /* up to 60 degrees */
    {
        double trial = mean_anom / fabs(1. - ecc);

        if(trial * trial > 6. * fabs(1. - ecc))   /* cubic term is dominant */
        {
            if(mean_anom < M_PI)
                trial = pow(6. * mean_anom, 1.0/3.0);
            else        /* hyperbolic w/ 5th & higher-order terms predominant */
                trial = asinh(mean_anom / ecc);
        }
        curr = trial;
    }

    double thresh = DEL_E_HYP * mean_anom;
    curr = log(2.*mean_anom/ecc+1.8); // taken from Burkardt & Danby, CeMec 31 (1983), 317-328
    double curr_abs = fabs(curr);
    err = (ecc * sinh(curr) - curr) - mean_anom;
    while(fabs(err) > thresh)
    {
        //n_iter++;
        //curr -= err / (ecc * cosh( curr) - 1.);
        //err = ecc * sinh( curr) - curr - mean_anom;

        n_iter++;
        if(curr_abs < .72 && ecc < 1.1)
        {
            // [e * sinh(E) - E] / E << 1, and/or e * cosh(E) - 1 << 1
            // so don't calculate it directly
            double curr2 = curr * curr;
            // relative error when omitting nth order term needs to be smaller than resolution 1.e-15:
            // .5 * E^2 > 1e15 * E^n/n!, i.e. E < (n!/2e15)^(1/(n-2))
            // n = 16: E < .72, n = 10: E < .08
         if(curr_abs > .08)
                curr -= err / ((ecc - 1) * cosh(curr) +
                        (((((((_1_16_15 * curr2 + 1.)
                        * _1_14_13 * curr2 + 1.)
                        * _1_12_11 * curr2 + 1.)
                        * _1_10_9 * curr2 + 1.)
                        * _1_8_7 * curr2 + 1.)
                        * _1_6_5 * curr2 + 1.)
                        * _1_4_3 * curr2 + 1.)
                        * .5 *curr2);
            else
                curr -= err / ((ecc - 1) * cosh(curr) +
                        ((((_1_10_9 * curr2 + 1.)
                        * _1_8_7 * curr2 + 1.)
                        * _1_6_5 * curr2 + 1.)
                        * _1_4_3 * curr2 + 1.)
                        * .5 *curr2);
            curr2 = curr * curr;
            curr_abs = fabs(curr);

            if(curr_abs > .08)
                err = ((ecc - 1) * sinh(curr) +
                        (((((((_1_17_16 * curr2 + 1.)
                        * _1_15_14 * curr2 + 1.)
                        * _1_13_12 * curr2 + 1.)
                        * _1_11_10 * curr2 + 1.)
                        * _1_9_8 * curr2 + 1.)
                        * _1_7_6 * curr2 + 1.)
                        * .05 * curr2 + 1.)
                        * _1_3_2 * curr2 * curr) - mean_anom;
            else
                err = ((ecc - 1) * sinh(curr) +
                        ((((_1_11_10 * curr2 + 1.)
                        * _1_9_8 * curr2 + 1.)
                        * _1_7_6 * curr2 + 1.)
                        * .05 * curr2 + 1.)
                        * _1_3_2 * curr2 * curr) - mean_anom;
        }
        else
        {
            curr -= err / (ecc * cosh(curr) - 1.);
            err = (ecc * sinh(curr) - curr) - mean_anom;
        }

        if(n_iter > KEPLER_ITE) // amended
        {
            #ifdef DEBUG_KEPLER
            fprintf(stderr,
                "[Warning] Aborting hyperbolic kepler solution after %d iterations, keeping error of %e (e=%e, M=%e, E_=%1.12e, sinh(E_)=%1.12e)\n",
                KEPLER_ITE, err,
                ecc, mean_anom, curr, sinh(curr));
            #endif
            break;
        }
    }

    return(is_negative ? -curr : curr);
}
