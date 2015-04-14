#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>

const double D_TIME_MIN = 1.1920928955078125e-07;
//const double D_TIME_MIN = 9.5367431640625e-07;
const double D_TIME_MAX = 0.125;
//const double DT_FIX = 0.0015625;

const double e2 = 1e-08;

typedef struct
{
    double m;
    double t, dt, old_dt;
    double rx, ry, rz;
    double vx, vy, vz;
    double prx, pry, prz;
    double pvx, pvy, pvz;
    double iprx, ipry, iprz;
    double ipvx, ipvy, ipvz;
    double a0x, a0y, a0z;
    double a1x, a1y, a1z;
    double oa0x, oa0y, oa0z;
    double oa1x, oa1y, oa1z;
    double a2x, a2y, a2z;
    double a3x, a3y, a3z;
} particle;

void prediction(double ITIME);
void evaluation();
void correction(double ITIME, bool check);
void force_calculation(particle &pi, particle &pj);
void print_all();
void next_c_time(double &c_time);
double get_magnitude(double x, double y, double z);
double get_timestep_normal(particle p);
void init_dt(double &CTIME);
double get_energy();
void save_old();

// Global structure
particle p[2];
double ETA_B = 0;

int main(int argc, char *argv[])
{
    double c_time = 1e6;

    p[0] = {0.0};
    p[1] = {0.0};

    /****************************************/
    // Masses m1 = m2
    //p[0].m  = 0.5;
    //p[1].m  = 0.5;
    // Position
    //double r_cte = 0.005; // a = 1, e = 0.99
    //double r_cte = 0.05; // a = 1, e = 0.9
    //double r_cte = 0.4; // a = 1, e = 0.2

    // Velocity
    //double vcte = 7.053367989832939; // a = 1, e = 0.99
    //double vcte = 2.17944947; // a = 1, e = 0.9
    //double vcte = 0.61237244; // a = 1, e = 0.2
    //p[0].rx =  -r_cte;
    //p[0].prx = -r_cte;
    //p[1].rx =   r_cte;
    //p[1].prx =  r_cte;

    //p[0].vy =   vcte;
    //p[0].pvy =  vcte;
    //p[1].vy =  -vcte;
    //p[1].pvy = -vcte;
    /*****************************************/

    //// Masses m1 = 2 m2
    //p[0].m  = 0.666666666666666;
    //p[1].m  = 0.333333333333333;

    //// e = 0.2
    //p[0].rx =  -0.26666667; p[0].prx = -0.26666667;
    //p[1].rx =   0.53333333; p[1].prx =  0.53333333;

    //p[0].vy =   0.40824829; p[0].pvy =  0.40824829;
    //p[1].vy =  -0.81649658; p[1].pvy = -0.81649658;

    //// e = 0.9
    //p[0].rx =  -0.03333333; p[0].prx = -0.03333333;
    //p[1].rx =   0.06666667; p[1].prx =  0.06666667;

    //p[0].vy =   1.45296631; p[0].pvy =  1.45296631;
    //p[1].vy =  -2.90593263; p[1].pvy = -2.90593263;

    //// e = 0.99
    //p[0].rx =  -0.003333333; p[0].prx = -0.003333333;
    //p[1].rx =   0.006666667; p[1].prx =  0.006666667;

    //p[0].vy =   4.70224533; p[0].pvy =  4.70224533;
    //p[1].vy =  -9.40449065; p[1].pvy = -9.40449065;

    /*****************************************/

    // Masses m1 = 10 m2
    //p[0].m  = 0.9090909090909091;
    //p[1].m  = 0.09090909090909091;

    //// e = 0.2
    //p[0].rx = -0.07272727; p[0].prx = -0.07272727;
    //p[1].rx = 0.72727273; p[1].prx = 0.72727273;

    //p[0].vy = 0.11134044; p[0].pvy = 0.11134044;
    //p[1].vy =  -1.11340443; p[1].pvy = -1.11340443;

    //// e = 0.9
    //p[0].rx =  -0.00909091; p[0].prx = -0.00909091;
    //p[1].rx =  0.09090909; p[1].prx = 0.09090909;

    //p[0].vy =  0.39626354; p[0].pvy = 0.39626354;
    //p[1].vy =  -3.9626354; p[1].pvy = -3.9626354;

    //// e = 0.99
    //p[0].rx =  -0.00090909; p[0].prx = -0.00090909;
    //p[1].rx =  0.00909091; p[1].prx = 0.00909091;

    //p[0].vy =  1.28243054; p[0].pvy = 1.28243054;
    //p[1].vy =  -12.82430544; p[1].pvy = -12.82430544;

    /*****************************************/

    // Masses m1 = 100 m2
    //p[0].m  = 0.9900990099009901;
    //p[1].m  = 0.009900990099009901;

    //// e = 0.2
    //p[0].rx =  -0.00792079; p[0].prx = -0.00792079;
    //p[1].rx =  0.79207921; p[1].prx = 0.79207921;

    //p[0].vy =  0.01212619; p[0].pvy = 0.01212619;
    //p[1].vy =  -1.21261868; p[1].pvy = -1.21261868;

    // e = 0.9
    //p[0].rx =  -0.0009901; p[0].prx = -0.0009901;
    //p[1].rx =  0.0990099; p[1].prx = 0.0990099;

    //p[0].vy =  0.04315742; p[0].pvy = 0.04315742;
    //p[1].vy =  -4.31574153; p[1].pvy = -4.31574153;

    //// e = 0.99
    //p[0].rx =  -9.90099010e-05; p[0].prx = -9.90099010e-05;
    //p[1].rx =  0.00990099; p[1].prx = 0.00990099;

    //p[0].vy =  0.13967065; p[0].pvy = 0.13967065;
    //p[1].vy =  -13.96706533; p[1].pvy = -13.96706533;
    /*****************************************/

    // Masses m1 = 10^6 m2
    //p[0].m  = 0.999999000001;
    //p[1].m  = 9.99999000001e-07;

    //// e = 0.2
    //p[0].rx =  -7.99999200e-07; p[0].prx = -7.99999200e-07;
    //p[1].rx =  0.7999992; p[1].prx = 0.7999992;

    //p[0].vy =  1.22474365e-06; p[0].pvy = 1.22474365e-06;
    //p[1].vy =  -1.22474365; p[1].pvy = -1.22474365;

    // e = 0.9
    //p[0].rx =  -9.99999000e-08; p[0].prx = -9.99999000e-08;
    //p[1].rx =  0.0999999; p[1].prx = 0.0999999;

    //p[0].vy =  4.35889458e-06; p[0].pvy = 4.35889458e-06;
    //p[1].vy =  -4.35889458; p[1].pvy = -4.35889458;

    //// e = 0.99
    //p[0].rx =  -9.99999000e-09; p[0].prx = -9.99999000e-09;
    //p[1].rx =  0.00999999; p[1].prx = 0.00999999;

    //p[0].vy =  1.41067219e-05; p[0].pvy = 1.41067219e-05;
    //p[1].vy =  -14.10672187; p[1].pvy = -14.10672187;
    /*****************************************/

    p[0].m = 9.765625000000000e-04;
    p[1].m = 9.765625000000000e-04;

    p[0].rx = 9.850490359358810e-05;
    p[0].prx = 9.850490359358810e-05;
    p[0].ry = -7.263274200450531e-04;
    p[0].pry = -7.263274200450531e-04;
    p[0].rz = -1.747488448700651e-03;
    p[0].prz = -1.747488448700651e-03;
    p[0].vx = 3.740929334219047e-01;
    p[0].pvx = 3.740929334219047e-01;
    p[0].vy = -1.857405959958897e-01;
    p[0].pvy = -1.857405959958897e-01;
    p[0].vz = 2.883360826020155e-01;
    p[0].pvz = 2.883360826020155e-01;

    p[1].rx = -9.850490359358810e-05;
    p[1].prx = -9.850490359358810e-05;
    p[1].ry = 7.263274200450809e-04;
    p[1].pry = 7.263274200450809e-04;
    p[1].rz = 1.747488448700651e-03;
    p[1].prz = 1.747488448700651e-03;
    p[1].vx = -3.740929334219047e-01;
    p[1].pvx = -3.740929334219047e-01;
    p[1].vy = 1.857405959958897e-01;
    p[1].pvy = 1.857405959958897e-01;
    p[1].vz = -2.883360826020155e-01;
    p[1].pvz = -2.883360826020155e-01;

    /*****************************************/

    double m  = p[0].m + p[1].m;
    double G  = 1.0;

    // Setup: end

    // First to get initial a0 and a1
    evaluation();
    init_dt(c_time);
    double ini_e = get_energy();
    long long int ite = 0;

    double mu = p[0].m * p[1].m / (p[0].m + p[1].m);
    double rx = p[1].rx - p[0].rx;
    double ry = p[1].ry - p[0].ry;
    double rz = p[1].rz - p[0].rz;

    double vx = p[1].vx - p[0].vx;
    double vy = p[1].vy - p[0].vy;
    double vz = p[1].vz - p[0].vz;

    double jx = ry * vz - rz * vy;
    double jy = rz * vx - rx * vz;
    double jz = rx * vy - ry * vx;

    double rr = sqrt(rx*rx + ry*ry + rz*rz);
    double vv = sqrt(vx*vx + vy*vy + vz*vz);
    double jj = sqrt(jx*jx + jy*jy + jz*jz);

    double j2 = jj*jj;
    double v2 = vv*vv;
    double m1m2 = p[0].m * p[1].m;

    double mu_std = G * m;
    double espec = v2 * 0.5 - mu_std/rr;
    double semimajor_ini = -mu_std / (2*espec);
    double ecc_ini = sqrt(1.0+2.0*espec*j2/(mu_std*mu_std));

    //std::cout << "mu_std: " << mu_std << std::endl;
    //std::cout << "rr: " << rr << std::endl;
    //std::cout << "vv: " << vv << std::endl;
    //std::cout << "jj: " << jj << std::endl;
    //std::cout << "a: " << semimajor_ini << std::endl;
    //std::cout << "espec: " << espec << std::endl;
    //std::cout << "ecc: " << ecc_ini << std::endl;
    //std::cout << "energy: " << ini_e << std::endl;
    printf("mu_std: %.15e\n", mu_std);
    printf("rr: %.15e\n", rr);
    printf("vv: %.15e\n", vv);
    printf("jj: %.15e\n", jj);
    printf("a: %.15e\n", semimajor_ini);
    printf("espec: %.15e\n", espec);
    printf("ecc: %.15e\n", ecc_ini);
    printf("energy: %.15e\n", ini_e);

    double ini_time = omp_get_wtime();
    //while (c_time < 2002 * 2 * 3.1415)
    c_time = 1.132813096046448e-01;
    while (c_time < 10 * 2 * 3.1415)
    {
        prediction(c_time);
        save_old();
        // P(EC)^n, n = 3
        for (int i = 0; i < (int)atoi(argv[1]); i++)
        {
            evaluation();
            correction(c_time, true);
        }

        for (int i = 0; i < 2; i++)
        {
            // Saving old timestep
            p[i].old_dt = p[i].dt;
            // Getting new timestep
            double first_dt  = get_timestep_normal(p[i]);
            double new_dt = 0.0;

            // Check even or odd timestep
            if ((int)(c_time / p[i].dt) % 2 == 0)
            {
                //std::cout << "Even" << std::endl;
                // We double it
                p[i].dt = p[i].old_dt * 2;
                evaluation();
                correction(c_time, false);
                double last_dt  = get_timestep_normal(p[i]);
                double avg = 0.5 * (first_dt + last_dt);

                if ( 2 * p[i].old_dt <= avg ||
                         p[i].old_dt <= avg &&
                     2 * p[i].old_dt >= avg)
                {
                    //std::cout << "Doubling it" << std::endl;
                    new_dt = p[i].dt;
                }
                else // If not, we keep the same value and try again
                {
                    p[i].dt = p[i].old_dt;
                    evaluation();
                    correction(c_time, false);
                    double last_dt  = get_timestep_normal(p[i]);
                    double avg = 0.5 * (first_dt + last_dt);

                    if ( 2 * p[i].old_dt <= avg ||
                             p[i].old_dt <= avg &&
                         2 * p[i].old_dt >= avg)
                    {
                        //std::cout << "Keeping the same" << std::endl;
                        new_dt = p[i].dt;
                    }
                    else
                    {
                        if ( p[i].old_dt * 0.5 <= avg &&
                             p[i].old_dt > avg)
                        {
                            //std::cout << "Halving it" << std::endl;
                            new_dt = p[i].dt * 0.5;
                        }
                        else
                        {
                            //std::cout << "...in any other case we halve" << std::endl;
                            new_dt = p[i].dt * 0.5;
                        }
                    }
                }
            }
            else
            {
                evaluation();
                correction(c_time, false);
                double last_dt  = get_timestep_normal(p[i]);
                double avg = 0.5 * (first_dt + last_dt);

                if ( 2 * p[i].old_dt <= avg ||
                         p[i].old_dt <= avg &&
                     2 * p[i].old_dt >= avg)
                {
                    //std::cout << "Keeping the same" << std::endl;
                    new_dt = p[i].dt;
                }
                else
                {
                    if ( p[i].old_dt * 0.5 <= avg &&
                         p[i].old_dt > avg)
                    {
                        //std::cout << "Halving it" << std::endl;
                        new_dt = p[i].dt * 0.5;
                    }
                    else
                    {
                        //std::cout << "...in any other case we halve" << std::endl;
                        new_dt = p[i].dt * 0.5;
                    }
                }
            }

            if (new_dt > D_TIME_MAX)
                new_dt = D_TIME_MAX;

            if (new_dt < D_TIME_MIN)
                new_dt = D_TIME_MIN;

            p[i].dt = new_dt;

            p[i].t = c_time;

            p[i].rx = p[i].prx;
            p[i].ry = p[i].pry;
            p[i].rz = p[i].prz;

            // Correcting velocity
            p[i].vx = p[i].pvx;
            p[i].vy = p[i].pvy;
            p[i].vz = p[i].pvz;

        }
        if (ite == 0)
            std::cout << "2" << std::endl;

        //if(ite%2 == 0)
        //{
            //print_all();
            //print_all();
            double end_e = get_energy();

            double mu = p[0].m * p[1].m / (p[0].m + p[1].m);
            double rx = p[1].rx - p[0].rx;
            double ry = p[1].ry - p[0].ry;
            double rz = p[1].rz - p[0].rz;

            double vx = p[1].vx - p[0].vx;
            double vy = p[1].vy - p[0].vy;
            double vz = p[1].vz - p[0].vz;


            double jx = ry * vz - rz * vy;
            double jy = rz * vx - rx * vz;
            double jz = rx * vy - ry * vx;

            double r = sqrt(rx*rx + ry*ry + rz*rz);
            double v = sqrt(vx*vx + vy*vy + vz*vz);
            double j = sqrt(jx*jx + jy*jy + jz*jz);

            double j2 = j*j;
            double v2 = v*v;
            double m1m2 = p[0].m * p[1].m;

            double binde = 0.5 * mu * v2 - m1m2/ r;

            double semimajor = -0.5 * m1m2 / binde;

            double jmax2 = semimajor * (p[0].m + p[1].m);
            double ecc = sqrt(1. - j2 / jmax2);

            printf("%.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
                    c_time/3.1415,
                    std::abs((semimajor - semimajor_ini)/semimajor_ini),
                    std::abs((ecc - ecc_ini)/ecc_ini),
                    r,
                    std::abs((end_e - ini_e)/ini_e),
                    omp_get_wtime() - ini_time,
                    end_e);
            getchar();
        //}
        next_c_time(c_time);
        ite++;

        //get_energy();
        //std::cout << c_time << std::endl;
    }
    return 0;
}

void save_old()
{
    int i;
    for (i = 0; i < 2; i++)
    {
        p[i].oa0x = p[i].a0x;
        p[i].oa0y = p[i].a0y;
        p[i].oa0z = p[i].a0z;

        p[i].oa1x = p[i].a1x;
        p[i].oa1y = p[i].a1y;
        p[i].oa1z = p[i].a1z;
    }
}

double get_energy()
{
    double pot = 0.0;
    double kin = 0.0;

    for (int i = 0; i < 2; i++)
    {
        double epot_tmp = 0.0;
        for (int j = i+1; j < 2; j++)
        {
            double rx = p[j].rx - p[i].rx;
            double ry = p[j].ry - p[i].ry;
            double rz = p[j].rz - p[i].rz;
            double r2 = rx*rx + ry*ry + rz*rz;

            epot_tmp -= (p[0].m * p[1].m) / sqrt(r2);
        }

        double vx = p[i].vx * p[i].vx;
        double vy = p[i].vy * p[i].vy;
        double vz = p[i].vz * p[i].vz;
        double v2 = vx + vy + vz;

        double ekin_tmp = 0.5 * p[i].m * v2;

        kin += ekin_tmp;
        pot += epot_tmp;
    }
    return kin + pot;
}

void init_dt(double &CTIME)
{
    double dt_min = 1e6;
    for (int i = 0; i < 2; i++)
    {
        double a2 = p[i].a0x * p[i].a0x +
                    p[i].a0y * p[i].a0y +
                    p[i].a0z * p[i].a0z;

        double j2 = p[i].a1x * p[i].a1x +
                    p[i].a1y * p[i].a1y +
                    p[i].a1z * p[i].a1z;

        double dt_0 = sqrt(a2)/sqrt(j2);

            // Finding the minimum
        if (dt_0 < dt_min)
             dt_min = dt_0;

        p[i].dt = 0.5 * D_TIME_MIN;
        p[i].t = 0.0;
    }

    ETA_B = D_TIME_MIN / (2.0 * dt_min);
    CTIME = 0.5 * D_TIME_MIN;
}

void next_c_time(double &c_time)
{
    double t1 = p[0].t + p[0].dt;
    double t2 = p[1].t + p[1].dt;

    if (t1 < t2)
        c_time = t1;
    else
        c_time = t2;
}

void prediction(double ITIME)
{
    for (int i = 0; i < 2; i++)
    {
        double dt  = ITIME - p[i].t;
        double dt2 = (dt  * dt);
        double dt3 = (dt2 * dt);

        p[i].prx = (dt3/6 * p[i].a1x) + (dt2/2 * p[i].a0x) + (dt * p[i].vx) + p[i].rx;
        p[i].pry = (dt3/6 * p[i].a1y) + (dt2/2 * p[i].a0y) + (dt * p[i].vy) + p[i].ry;
        p[i].prz = (dt3/6 * p[i].a1z) + (dt2/2 * p[i].a0z) + (dt * p[i].vz) + p[i].rz;

        p[i].pvx = (dt2/2 * p[i].a1x) + (dt * p[i].a0x) + p[i].vx;
        p[i].pvy = (dt2/2 * p[i].a1y) + (dt * p[i].a0y) + p[i].vy;
        p[i].pvz = (dt2/2 * p[i].a1z) + (dt * p[i].a0z) + p[i].vz;

        p[i].iprx = p[i].prx;
        p[i].ipry = p[i].pry;
        p[i].iprz = p[i].prz;

        p[i].ipvx = p[i].pvx;
        p[i].ipvy = p[i].pvy;
        p[i].ipvz = p[i].pvz;
    }
}

void evaluation()
{
    for (int i = 0; i < 2; i++)
    {
        p[i].a0x = 0.0;
        p[i].a0y = 0.0;
        p[i].a0z = 0.0;

        p[i].a1x = 0.0;
        p[i].a1y = 0.0;
        p[i].a1z = 0.0;

        for (int j = 0; j < 2; j++)
        {
            if(i == j) continue;
            force_calculation(p[i], p[j]);
        }
    }
}

void correction(double ITIME, bool check)
{
    for (int i = 0; i < 2; i++)
    {
        double dt1 = p[i].dt;
        double dt2 = dt1 * dt1;
        double dt3 = dt2 * dt1;
        double dt4 = dt2 * dt2;
        double dt5 = dt4 * dt1;

        // Acceleration 2nd derivate
        p[i].a2x = (-6 * (p[i].oa0x - p[i].a0x ) - dt1 * (4 * p[i].oa1x + 2 * p[i].a1x) ) / dt2;
        p[i].a2y = (-6 * (p[i].oa0y - p[i].a0y ) - dt1 * (4 * p[i].oa1y + 2 * p[i].a1y) ) / dt2;
        p[i].a2z = (-6 * (p[i].oa0z - p[i].a0z ) - dt1 * (4 * p[i].oa1z + 2 * p[i].a1z) ) / dt2;

        // Acceleration 3rd derivate
        p[i].a3x = (12 * (p[i].oa0x - p[i].a0x ) + 6 * dt1 * (p[i].oa1x + p[i].a1x) ) / dt3;
        p[i].a3y = (12 * (p[i].oa0y - p[i].a0y ) + 6 * dt1 * (p[i].oa1y + p[i].a1y) ) / dt3;
        p[i].a3z = (12 * (p[i].oa0z - p[i].a0z ) + 6 * dt1 * (p[i].oa1z + p[i].a1z) ) / dt3;

        if (check)
        {
            // Correcting position
            p[i].prx = p[i].iprx + (dt4/24)*p[i].a2x + (dt5/120)*p[i].a3x;
            p[i].pry = p[i].ipry + (dt4/24)*p[i].a2y + (dt5/120)*p[i].a3y;
            p[i].prz = p[i].iprz + (dt4/24)*p[i].a2z + (dt5/120)*p[i].a3z;

            // Correcting velocity
            p[i].pvx = p[i].ipvx + (dt3/6)*p[i].a2x + (dt4/24)*p[i].a3x;
            p[i].pvy = p[i].ipvy + (dt3/6)*p[i].a2y + (dt4/24)*p[i].a3y;
            p[i].pvz = p[i].ipvz + (dt3/6)*p[i].a2z + (dt4/24)*p[i].a3z;
        }
    }
}

double get_magnitude(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}

double get_timestep_normal(particle p)
{
    // Calculating a_{1,i}^{(2)} = a_{0,i}^{(2)} + dt * a_{0,i}^{(3)}
    double ax1_2 = p.a2x + p.dt * p.a3x;
    double ay1_2 = p.a2y + p.dt * p.a3y;
    double az1_2 = p.a2z + p.dt * p.a3z;

    // |a_{1,i}|
    double abs_a1 = get_magnitude(p.a0x, p.a0y, p.a0z);
    // |j_{1,i}|
    double abs_j1 = get_magnitude(p.a1x, p.a1y, p.a1z);
    // |j_{1,i}|^{2}
    double abs_j12  = abs_j1 * abs_j1;
    // a_{1,i}^{(3)} = a_{0,i}^{(3)} because the 3rd-order interpolation
    double abs_a1_3 = get_magnitude(p.a3x, p.a3y, p.a3z);
    // |a_{1,i}^{(2)}|
    double abs_a1_2 = get_magnitude(ax1_2, ay1_2, az1_2);
    // |a_{1,i}^{(2)}|^{2}
    double abs_a1_22  = abs_a1_2 * abs_a1_2;

    double normal_dt = sqrt(ETA_B * ((abs_a1 * abs_a1_2 + abs_j12) / (abs_j1 * abs_a1_3 + abs_a1_22)));
    return normal_dt;
}

void force_calculation(particle &pi, particle &pj)
{
    double rx = pj.prx - pi.prx;
    double ry = pj.pry - pi.pry;
    double rz = pj.prz - pi.prz;

    double vx = pj.pvx - pi.pvx;
    double vy = pj.pvy - pi.pvy;
    double vz = pj.pvz - pi.pvz;

    double r2     = rx*rx + ry*ry + rz*rz;
    double rinv   = 1.0/sqrt(r2);
    double r2inv  = rinv  * rinv;
    double r3inv  = r2inv * rinv;
    double r5inv  = r2inv * r3inv;
    double mr3inv = r3inv * pj.m;
    double mr5inv = r5inv * pj.m;

    double rv = rx*vx + ry*vy + rz*vz;

    pi.a0x += (rx * mr3inv);
    pi.a0y += (ry * mr3inv);
    pi.a0z += (rz * mr3inv);

    pi.a1x += (vx * mr3inv - (3 * rv ) * rx * mr5inv);
    pi.a1y += (vy * mr3inv - (3 * rv ) * ry * mr5inv);
    pi.a1z += (vz * mr3inv - (3 * rv ) * rz * mr5inv);
}

void print_all()
{
    for (int i = 0; i < 2; i++)
    {
        std::cout << std::setw(2) << i << " ";
        std::cout.precision(1);
        std::cout << std::scientific;
        std::cout.precision(3);
        //std::cout << std::setw(10) << p[i].t   << " ";
        //std::cout << std::setw(3) << p[i].m   << " ";
        //std::cout << std::setw(10) << p[i].dt  << " ";
        std::cout << std::setw(10) << p[i].rx  << " ";
        std::cout << std::setw(10) << p[i].ry  << " ";
        std::cout << std::setw(10) << p[i].rz  << " ";
        std::cout << std::endl;
    }
}

