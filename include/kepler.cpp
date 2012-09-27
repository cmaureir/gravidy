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
                       double m,   double m0,  double dt, int i)
{

/*
     * Orbital parameters
     */
    double jx, jy, jz;
    double ex, ey, ez;
    double a, b, w;
    double ax,ay,az,bx,by,bz;
    double ecc, e_anomaly, m_anomaly;
    double r_const, v_const;
    double r = sqrt(*rx * *rx + *ry * *ry + *rz * *rz);
    double v;
    double cos_e, sin_e, e_tmp, ecc_const, cos_const;

    printf("Starting with the values:\n \
            Position vector: (%.10f,%.10f,%.10f)\n \
            Velocity vector: (%.10f,%.10f,%.10f)\n \
            m: %.10f\n m0: %.10f\n",*rx,*ry,*rz,*vx,*vy,*vz,m,m0);

    // Angular momentum
    // j = r x v
    jx = *ry * *vz + *rz * *vy;
    jy = *rz * *vx + *rx * *vz;
    jz = *rx * *vy + *ry * *vx;

    printf("Angular momentum vector (j): %.10f %.10f %.10f\n", jx, jy, jz);
    // Runge-Lenz-vector
    // e = { (v x j) / (G * m) }  - { r / |r| }
    ex = *vy * jz + jy * *vz;
    ey = *vx * jz + jx * *vz;
    ez = *vx * jy + jx * *vy;

    ex = ex/(G * m0) - *rx/r;
    ey = ey/(G * m0) - *ry/r;
    ez = ez/(G * m0) - *rz/r;

    printf("Runge-Lenz-vector (e): %.10f %.10f %.10f\n", ex, ey, ez);

    // Eccentricity
    ecc = sqrt(ex*ex + ey*ey + ez*ez);
    printf("Eccentricity (ecc): %.10f\n", ecc);

    // Semi-major axis
    // a = ( j * j ) / (G * m * | 1 - ecc^2 | )
    a = (jx*jx + jy*jy + jz*jz) / (G * m0 * fabs(1 - (ex*ex + ey*ey + ez*ez)));
    printf("Semi-major axis (a): %.10f\n", a);


    // Frequency of the orbit
    // w = sqrt(( G * m )/ a^3 )
    w = sqrt((G * m0) / (a*a*a));
    printf("Frequency of the orbit (w): %.10f\n", w);

    // Semi-major vector
    ax = a * ex/ecc;
    ay = a * ey/ecc;
    az = a * ez/ecc;
    printf("Semi-major vector: %.10f %.10f %.10f\n", ax, ay, az);

    // Semi-minor vector
    bx = (jy * ez + ey * jz);
    by = (jx * ez + ex * jz);
    bz = (jx * ey + ex * jy);

    double b_vmag = sqrt(bx*bx + by*by + bz*bz);
    bx *= b / b_vmag;
    by *= b / b_vmag;
    bz *= b / b_vmag;

    printf("Semi-minor vector: (%.10f) %.10f %.10f %.10f\n", b_vmag, bx, by, bz);

    // Semi-minor axis
    // b = a * sqrt ( 1 - ecc^2 )
    b = a * sqrt(fabs(1 - ecc*ecc))/b_vmag;
    printf("Semi-minor axis (b): %.10f\n", b);


    // End of the variables calculation.

    if(ecc < 1)
    {

    //    std::cout << "Elliptical orbit, with e=" << ecc << std::endl;

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

        // TO DO: Why?
        while(m_anomaly >= TWOPI)
            m_anomaly -= TWOPI;

        e_anomaly = solve_kepler(m_anomaly, ecc);

        cos_e = cos(e_anomaly);
        sin_e = sin(e_anomaly);

        r_const = cos_e - ecc;
        v_const = w / (1.0 - ecc * cos_e);

        // TO DO: Why?
        if(ecc > 0.99)
        {
            e_tmp = (e_anomaly > 2.0 * M_PI - 1e-3) ? e_anomaly - 2.0 * M_PI : e_anomaly;
            if(e_tmp < 1e-3)
            {
                e_tmp *= e_tmp;
                ecc_const = (jx*jx + jy*jy + jz*jz)/(m0 * a * (1 + ecc));
                cos_const = -0.5 * e_tmp * (1 - e_tmp / 12.0 * (1 - e_tmp / 30.0));
                r_const = ecc_const + cos_const;
                v_const = w / (ecc_const - ecc * cos_const);
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

    //printf("New Position vector: (%.10f,%.10f,%.10f)\n \
            New Velocity vector: (%.10f,%.10f,%.10f)\n", *rx,*ry,*rz,*vx,*vy,*vz);

    }
    else
    {
        // Hyperbola or Parabola orbit
        //std::cerr << "Nothing to do here!" << std::endl;
        // calculate eccentric anomaly e at t+dt
        std::cout << "---------> Hyperbolic or Parabolic orbit, with e=" << ecc << std::endl;
        e_anomaly = (a + r) / (ecc * a);
        printf("e_anomaly (1): %.10f\n", e_anomaly);
        if(e_anomaly < 1.0)
            e_anomaly = 0.0;
        else if((*rx * *vx + *ry * *vy + *rz * *vz) < 0)
            e_anomaly = -acosh(e_anomaly);
        else
            e_anomaly = acosh(e_anomaly);

        printf("e_anomaly (2): %.10f\n", e_anomaly);
        printf("kepler args: %.10f %.10f\n", ecc, ecc * sinh(e_anomaly) - e_anomaly + dt * w);

        e_anomaly = kepler(ecc, ecc * sinh(e_anomaly) - e_anomaly + dt * w);
        printf("e_anomaly (3): %.10f\n", e_anomaly);
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

    printf("New Position vector: (%.10f,%.10f,%.10f)\n \
            New Velocity vector: (%.10f,%.10f,%.10f)\n", *rx,*ry,*rz,*vx,*vy,*vz);
    }

    r = sqrt(*rx * *rx + *ry * *ry + *rz * *rz);
    v = sqrt(*vx * *vx + *vy * *vy + *vz * *vz);
    v = sqrt(jx*jx + jy*jy + jz*jz) / (r * v * sin(acos((*rx * *vx + *ry * *vy + *rz * *vz)/(r * v))));

    // total motion relative to fix central mass
    h_p_r[i].x = h_p_r[0].x + *rx;
    h_p_r[i].y = h_p_r[0].y + *ry;
    h_p_r[i].z = h_p_r[0].z + *rz;

    h_p_v[i].x = h_p_v[0].x + *vx;
    h_p_v[i].y = h_p_v[0].y + *vy;
    h_p_v[i].z = h_p_v[0].z + *vz;

//
//    if(out_a != NULL)
//    {
//        _1_r2 = 1. / scal_prod(r_, r_);
//        afact = - m_c * _1_r2 * sqrt(_1_r2);
//        //printf("4  %e %e %e\n", *(out_a), *(out_a+1), *(out_a+2));
//        for(i = 0; i < DIMENSIONS; i++)
//            out_a[i] = afact * r_[i];
//            if(out_a_ != NULL)
//            {
//                v_r_ = scal_prod(v_, r_);
//                for(i = 0; i < DIMENSIONS; i++)
//                    out_a_[i] = afact * (v_[i] - 3 * _1_r2 * v_r_ * r_[i]);
//                    if(out_a_2 != NULL)
//                    {
//                        v_v_ = scal_prod(v_, v_);
//                        r_a_ = scal_prod(r_, out_a);
//                        for(i = 0; i < DIMENSIONS; i++)
//                            out_a_2[i] = afact * (out_a[i] - 3. * _1_r2 * (v_r_ * (2. * v_[i] - 5. * v_r_ * r_[i] * _1_r2)
//                                         + (v_v_ + r_a_) * r_[i]));
//                        if(out_a_3 != NULL)
//                        {
//                            v_a_  = scal_prod(v_, out_a);
//                            r_a__  = scal_prod(r_, out_a_);
//                            for(i = 0; i < DIMENSIONS; i++)
//                                out_a_3[i] = afact * (out_a_[i]
//                                            - 3. * _1_r2 * (3. * v_r_ * out_a[i]
//                                            + 3. * (v_v_ + r_a_)
//                                            * (v_[i] - 5. * v_r_ * _1_r2 * r_[i])
//                                            + (3. * v_a_ + r_a__) * r_[i]
//                                            + v_r_ * v_r_ * _1_r2
//                                            * (-15. * v_[i] + 35. * v_r_ * _1_r2 * r_[i])));
//                        }
//                    }
//            }
//   }
//
}

double solve_kepler(double m_anomaly, double ecc)
{

// Solving Kepler's equation
double e_new = ecc > 0.8 ? PI : m_anomaly;
double d;
int e_iter = 0;

    for(d = e_new - ecc * sin(e_new) - m_anomaly;
        fabs(d) > DEL_E;
        d = e_new - ecc * sin(e_new) - m_anomaly)
    {
        if(e_iter-1 >= KEPLER_ITE)
        {
            break;
        }
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
    double curr, err;
    int is_negative = 0, n_iter = 0;

    if(!mean_anom)
        return(0.);

    if(ecc < .3)     /* low-eccentricity formula from Meeus,  p. 195 */
    {
        curr = atan2(sin(mean_anom), cos(mean_anom) - ecc);
        /* one correction step,  and we're done */
        err = curr - ecc * sin(curr) - mean_anom;
        curr -= err / (1. - ecc * cos(curr));
        printf("(1) %.10f\n", curr);
        return(curr);
    }

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
    if(ecc < 1.)
    {
        err = (curr - ecc * sin(curr)) - mean_anom;
        while(fabs(err) > thresh)
        {
            n_iter++;
            curr -= err / (1. - ecc * cos(curr));
            err = (curr - ecc * sin(curr)) - mean_anom;

            if(n_iter > KEPLER_ITE) // amended
            {
                fprintf(stderr,
                    "#### ! Aborting kepler solution after %d iterations, keeping error of %e (e=%e, M=%e, E_=%e) ####\n",
                    KEPLER_ITE, err,
                    ecc, mean_anom, curr);
                break;
            }
        }
    }
    else
    {
        curr = log(2.*mean_anom/ecc+1.8); // taken from Burkardt & Danby, CeMec 31 (1983), 317-328
        double curr_abs = fabs(curr);
        err = (ecc * sinh(curr) - curr) - mean_anom;
        while(fabs(err) > thresh)
        {
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
                fprintf(stderr,
                    "### Aborting hyperbolic kepler solution after %d iterations, keeping error of %e (e=%e, M=%e, E_=%1.12e, sinh(E_)=%1.12e)\n",
                    KEPLER_ITE, err,
                    ecc, mean_anom, curr, sinh(curr));

                break;
            }
        }
     }

    return(is_negative ? -curr : curr);
}
