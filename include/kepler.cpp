#include "kepler.hpp"

void kepler_prediction(double *rx, double *ry, double *rz,
                       double *vx, double *vy, double *vz,
                       double m,   double m0,  double dt)
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
    double cos_e, sin_e, e_tmp, ecc_const, cos_const;

    // Angular momentum
    jx = *ry * *vz + *vy * *rz;
    jy = *rx * *vz + *vx * *rz;
    jz = *rx * *vy + *vx * *ry;

    // Runge-Lenz-vector
    ex = *vy * jz + jy * *vz;
    ey = *vx * jz + jx * *vz;
    ez = *vx * jy + jx * *vy;

    ex = ex/(G * m0) - *rx/r;
    ey = ey/(G * m0) - *ry/r;
    ez = ez/(G * m0) - *rz/r;

    // Semi-major axis
    a = (jx*jx + jy*jy + jz*jz) / (G * m0 * fabs(1 - (ex*ex + ey*ey + ez*ez)));
    // Semi-minor axis
    b = a * sqrt(1-(ex*ex + ey*ey + ez*ez));

    // Frequency of the orbit
    w = sqrt((G * m0) / (a*a*a));

    // Eccentricity
    ecc = sqrt(ex*ex + ey*ey + ez*ez);

    // Semi-major vector
    ax = a * ex/ecc;
    ay = a * ey/ecc;
    az = a * ez/ecc;

    // Semi-minor vector
    bx = (jy * ez + ey * jz);
    by = (jx * ez + ex * jz);
    bz = (jx * ey + ex * jy);
    double b_vmag = sqrt(bx*bx + by*by + bz*bz);
    bx *= b / b_vmag;
    by *= b / b_vmag;
    bz *= b / b_vmag;

    /*
     * Working with elliptical orbits
     * TO DO: Case with Parabola or Hyperbola
     */
    if(ecc < 1)
    {

        std::cout << "Elliptical orbit, with e=" << ecc << std::endl;

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

        e_anomaly = e_new;

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

    }
    else
    {
        // Hyperbola or Parabola orbit
        std::cerr << "Nothing to do here!" << std::endl;
        //// TO DO
        //// calculate eccentric anomaly e at t+dt
        //e = (a + v_abs(r_)) / (ecc * a);
        //if(e < 1.0) e = .0;
        //else if(scal_prod(r_, v_) < 0) e = -acosh(e);
        //else e = acosh(e);

        //e = kepler(ecc, ecc * sinh(e) - e + dt * omega);
        //cosp = cosh(e);
        //sinp = sinh(e);
        //de_dt = omega / (ecc * cosp - 1.);
        //for(i = 0; i < DIMENSIONS; i++)
        //{
        //    r_[i] =   a_[i] * (ecc - cosp)  + b_[i] * sinp;  // new location
        //    v_[i] = (-a_[i] * sinp          + b_[i] * cosp) * de_dt;  // direction of v only
        //}
    }

//  // get |v_| from j_ = r_ x v_
//    v = v_abs(v_);
//    r = v_abs(r_);
//    v = v_abs(j_) / (r * v * sin(acos(scal_prod(r_, v_)/ (r * v))));
//
//    for(i = 0; i < DIMENSIONS; i++)
//    {
//        //v_[i] *= v;
//        // total motion relative to fix central mass
//        p1->xp[i] = p0->xp[i] + r_[i];
//        p1->vp[i] = p0->vp[i] + v_[i];
//    }
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

