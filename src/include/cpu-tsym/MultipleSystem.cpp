#include "MultipleSystem.hpp"

MultipleSystem::MultipleSystem(NbodySystem *ns, NbodyUtils *nu)
{
    this->ns = ns;
    this->nu = nu;
    this->members = 0;
    this->ETA_B = 0;
}

MultipleSystem::~MultipleSystem()
{
    // ...
}

void MultipleSystem::adjust_particles(SParticle sp)
{
    for(int i = 0; i < members; i++)
    {
        parts[i].r.x -= sp.r.x;
        parts[i].r.y -= sp.r.y;
        parts[i].r.z -= sp.r.z;

        parts[i].p.r[0] = parts[i].r.x;
        parts[i].p.r[1] = parts[i].r.y;
        parts[i].p.r[2] = parts[i].r.z;

        parts[i].v.x -= sp.v.x;
        parts[i].v.y -= sp.v.y;
        parts[i].v.z -= sp.v.z;

        parts[i].p.v[0] = parts[i].v.x;
        parts[i].p.v[1] = parts[i].v.y;
        parts[i].p.v[2] = parts[i].v.z;

        //parts[i].f.a[0] -= sp.f.a[0];
        //parts[i].f.a[1] -= sp.f.a[1];
        //parts[i].f.a[2] -= sp.f.a[2];

        //parts[i].f.a1[0] -= sp.f.a1[0];
        //parts[i].f.a1[1] -= sp.f.a1[1];
        //parts[i].f.a1[2] -= sp.f.a1[2];
    }
}

void MultipleSystem::add_particle(int id)
{
    if (members < MAX_MEMBERS)
    {
        // New particle
        MParticle new_member;

        new_member.id = id;
        new_member.r  = ns->h_r[id];
        new_member.v  = ns->h_v[id];
        new_member.f  = ns->h_f[id];
        new_member.a2 = ns->h_a2[id];
        new_member.a3 = ns->h_a3[id];
        new_member.t  = ns->h_t[id];
        new_member.dt = ns->h_dt[id];

        new_member.p.r[0] = new_member.r.x;
        new_member.p.r[1] = new_member.r.y;
        new_member.p.r[2] = new_member.r.z;

        new_member.p.v[0] = new_member.v.x;
        new_member.p.v[1] = new_member.v.y;
        new_member.p.v[2] = new_member.v.z;

        new_member.p.m    = new_member.r.w;

        parts[members] = new_member;
        members += 1;
    }
    else
    {
        std::cerr << "Reaching MAX_MEMBERS in a MultipleSystem" << std::endl;
    }
}

void MultipleSystem::print()
{
    for (int i = 0; i < members; i++)
    {
        MParticle part = parts[i];

        std::cout << std::setw(2) << i << " ";
        std::cout.precision(1);
        std::cout << std::scientific;
        std::cout.precision(15);
        std::cout << std::setw(10) << part.t   << " ";
        //std::cout << std::setw(3)  << part.r.w << " ";
        std::cout << std::setw(10) << part.dt  << " ";
        std::cout << std::setw(10) << part.r.x  << " ";
        std::cout << std::setw(10) << part.r.y  << " ";
        std::cout << std::setw(10) << part.r.z  << " ";
        std::cout << std::setw(10) << part.v.x  << " ";
        std::cout << std::setw(10) << part.v.y  << " ";
        std::cout << std::setw(10) << part.v.z  << " ";
        //std::cout << std::setw(10) << part.p.r[0] << " ";
        //std::cout << std::setw(10) << part.p.r[1] << " ";
        //std::cout << std::setw(10) << part.p.r[2] << " ";
        //std::cout << std::setw(10) << part.p.v[0] << " ";
        //std::cout << std::setw(10) << part.p.v[1] << " ";
        //std::cout << std::setw(10) << part.p.v[2] << " ";
        std::cout << std::setw(10) << part.f.a[0] << " ";
        std::cout << std::setw(10) << part.f.a[1] << " ";
        std::cout << std::setw(10) << part.f.a[2] << " ";
        //std::cout << std::setw(10) << part.f.a1[0] << " ";
        //std::cout << std::setw(10) << part.f.a1[1] << " ";
        //std::cout << std::setw(10) << part.f.a1[2] << " ";
        //std::cout << std::setw(10) << part.old_f.a[0] << " ";
        //std::cout << std::setw(10) << part.old_f.a[1] << " ";
        //std::cout << std::setw(10) << part.old_f.a[2] << " ";
        //std::cout << std::setw(10) << part.old_f.a1[0] << " ";
        //std::cout << std::setw(10) << part.old_f.a1[1] << " ";
        //std::cout << std::setw(10) << part.old_f.a1[2] << " ";
        std::cout << std::endl;
    }
}

void MultipleSystem::next_itime(double &CTIME)
{
    double new_time = 1e6;

    for (int i = 0; i < members; i++)
    {
        MParticle part = parts[i];
        double t = part.t + part.dt;
        if (t < new_time)
            new_time = t;
    }
    CTIME = new_time;
}

SParticle MultipleSystem::get_center_of_mass(MParticle p1, MParticle p2)
{

    SParticle pcm;

    pcm.r.w   = (p1.r.w + p2.r.w);
    double minv = 1.0/pcm.r.w;

    // Calculate CoM position
    pcm.r.x = ((p1.r.x * p1.r.w) + (p2.r.x * p2.r.w))*minv;
    pcm.r.y = ((p1.r.y * p1.r.w) + (p2.r.y * p2.r.w))*minv;
    pcm.r.z = ((p1.r.z * p1.r.w) + (p2.r.z * p2.r.w))*minv;

    // Calculate CoM velocity
    pcm.v.x = ((p1.v.x * p1.r.w) + (p2.v.x * p2.r.w))*minv;
    pcm.v.y = ((p1.v.y * p1.r.w) + (p2.v.y * p2.r.w))*minv;
    pcm.v.z = ((p1.v.z * p1.r.w) + (p2.v.z * p2.r.w))*minv;

    // Calculate CoM acceleration
    pcm.f.a[0] = ((p1.f.a[0] * p1.r.w) + (p2.f.a[0] * p2.r.w))*minv;
    pcm.f.a[1] = ((p1.f.a[1] * p1.r.w) + (p2.f.a[1] * p2.r.w))*minv;
    pcm.f.a[2] = ((p1.f.a[2] * p1.r.w) + (p2.f.a[2] * p2.r.w))*minv;

    // Calculate CoM acceleration first derivative
    pcm.f.a1[0] = ((p1.f.a1[0] * p1.r.w) + (p2.f.a1[0] * p2.r.w))*minv;
    pcm.f.a1[1] = ((p1.f.a1[1] * p1.r.w) + (p2.f.a1[1] * p2.r.w))*minv;
    pcm.f.a1[2] = ((p1.f.a1[2] * p1.r.w) + (p2.f.a1[2] * p2.r.w))*minv;

    //pcm.old = pcm.f;

    return pcm;
}

void MultipleSystem::prediction(double CTIME)
{

    for (int i = 0; i < members; i++)
    {
        // Accesing system particle
        MParticle part = parts[i];

        double dt  = CTIME - part.t;
        double dt2 = (dt  * dt);
        double dt3 = (dt2 * dt);

        double dt2c = dt2/2.0;
        double dt3c = dt3/6.0;

        // Temporal predictor variable
        Predictor pp;

        // Access particle information only once
        Forces ff  = part.f;
        double4 vv = part.v;
        double4 rr = part.r;

        pp.r[0] = (dt3c * ff.a1[0]) + (dt2c * ff.a[0]) + (dt * vv.x) + rr.x;
        pp.r[1] = (dt3c * ff.a1[1]) + (dt2c * ff.a[1]) + (dt * vv.y) + rr.y;
        pp.r[2] = (dt3c * ff.a1[2]) + (dt2c * ff.a[2]) + (dt * vv.z) + rr.z;

        pp.v[0] = (dt2c * ff.a1[0]) + (dt * ff.a[0]) + vv.x;
        pp.v[1] = (dt2c * ff.a1[1]) + (dt * ff.a[1]) + vv.y;
        pp.v[2] = (dt2c * ff.a1[2]) + (dt * ff.a[2]) + vv.z;

        pp.m = rr.w;

        // Update real predictor
        parts[i].p = pp;

        // Saving initial value
        parts[i].p0 = pp;
    }
}

void MultipleSystem::force_calculation(Predictor pi, Predictor pj, Forces &fi)
{
    double rx = pj.r[0] - pi.r[0];
    double ry = pj.r[1] - pi.r[1];
    double rz = pj.r[2] - pi.r[2];

    double vx = pj.v[0] - pi.v[0];
    double vy = pj.v[1] - pi.v[1];
    double vz = pj.v[2] - pi.v[2];

    double r2   = rx*rx + ry*ry + rz*rz;
    double rinv = 1.0/sqrt(r2);

    double rv = rx*vx + ry*vy + rz*vz;

    double r2inv  = rinv  * rinv;
    double r3inv  = r2inv * rinv;
    double r5inv  = r2inv * r3inv;
    double mr3inv = r3inv * pj.m;
    double mr5inv = r5inv * pj.m;

    fi.a[0] += (rx * mr3inv);
    fi.a[1] += (ry * mr3inv);
    fi.a[2] += (rz * mr3inv);

    fi.a1[0] += (vx * mr3inv - (3 * rv ) * rx * mr5inv);
    fi.a1[1] += (vy * mr3inv - (3 * rv ) * ry * mr5inv);
    fi.a1[2] += (vz * mr3inv - (3 * rv ) * rz * mr5inv);
}

void MultipleSystem::save_old()
{
    for (int i = 0; i < members; i++)
    {
        parts[i].old_f = parts[i].f;
    }
}

void MultipleSystem::evaluation(int *nb_list)
{
    int i, j;
    //#pragma omp parallel for private(i,j)
    for (i = 0; i < members; i++)
    {
        // Initialization of all the struct to zero
        Forces ff = {0};
        Predictor pi = parts[i].p;

        //#pragma omp parallel for
        for (j = 0; j < members; j++)
        {
            if(i == j) continue;
            force_calculation(pi, parts[j].p, ff);
        }

        parts[i].f = ff;
    }
    // Perturbers
    //for (i = 0; i < members; i++)
    //{
    //    Forces ff = {0};
    //    Predictor pi = parts[i].p;

    //    #pragma omp parallel for
    //    for (j = 0; j < ns->h_f[parts[0].id].nb; j++)
    //    {
    //        if(i == j) continue;
    //        force_calculation(pi, ns->h_p[nb_list[j]], ff);
    //    }

    //    for (j = 0; j < ns->n; j++)
    //    {
    //        if(i == j) continue;
    //        if ( j != parts[0].id)
    //            force_calculation(pi, ns->h_p[j], ff);
    //    }

    //    printf("WHOLE %.15e %.15e %.15e\n", ff.a[0], ff.a[1], ff.a[2]);

    //    // Write on the Forces array
    //    parts[i].f.a[0] += ff.a[0];
    //    parts[i].f.a[1] += ff.a[1];
    //    parts[i].f.a[2] += ff.a[2];

    //    parts[i].f.a1[0] += ff.a1[0];
    //    parts[i].f.a1[1] += ff.a1[1];
    //    parts[i].f.a1[2] += ff.a1[2];
    //}
}

void MultipleSystem::correction(double CTIME, bool check)
{

    for (int i = 0; i < members; i++)
    {
        // Accesing system particle
        MParticle part = parts[i];

        double dt1 = part.dt;
        double dt2 = dt1 * dt1;
        double dt3 = dt2 * dt1;
        double dt4 = dt2 * dt2;
        double dt5 = dt4 * dt1;

        double dt3c = dt3/6.0;
        double dt4c = dt4/24.0;
        double dt5c = dt5/120.0;

        double dt2inv = 1.0/dt2;
        double dt3inv = 1.0/dt3;

        // Load of the data that we will use only once, and not in every line
        Forces fi   =  part.f;
        Forces oldi =  part.old_f;
        Predictor p0 = part.p0;
        double3 a2i;
        double3 a3i;

        // Acceleration 2nd derivate
        a2i.x = (-6 * (oldi.a[0] - fi.a[0]) - dt1 * (4 * oldi.a1[0] + 2 * fi.a1[0])) * dt2inv;
        a2i.y = (-6 * (oldi.a[1] - fi.a[1]) - dt1 * (4 * oldi.a1[1] + 2 * fi.a1[1])) * dt2inv;
        a2i.z = (-6 * (oldi.a[2] - fi.a[2]) - dt1 * (4 * oldi.a1[2] + 2 * fi.a1[2])) * dt2inv;

        // Acceleration 3rd derivate
        a3i.x = (12 * (oldi.a[0] - fi.a[0] ) + 6 * dt1 * (oldi.a1[0] + fi.a1[0])) * dt3inv;
        a3i.y = (12 * (oldi.a[1] - fi.a[1] ) + 6 * dt1 * (oldi.a1[1] + fi.a1[1])) * dt3inv;
        a3i.z = (12 * (oldi.a[2] - fi.a[2] ) + 6 * dt1 * (oldi.a1[2] + fi.a1[2])) * dt3inv;

        if (check)
        {
            // Correcting position
            parts[i].p.r[0] = p0.r[0] + dt4c * a2i.x + dt5c * a3i.x;
            parts[i].p.r[1] = p0.r[1] + dt4c * a2i.y + dt5c * a3i.y;
            parts[i].p.r[2] = p0.r[2] + dt4c * a2i.z + dt5c * a3i.z;

            // Correcting velocity
            parts[i].p.v[0] = p0.v[0] + dt3c * a2i.x + dt4c * a3i.x;
            parts[i].p.v[1] = p0.v[1] + dt3c * a2i.y + dt4c * a3i.y;
            parts[i].p.v[2] = p0.v[2] + dt3c * a2i.z + dt4c * a3i.z;
        }

        // Write on the arrays
        parts[i].a2 = a2i;
        parts[i].a3 = a3i;

    }

}

double MultipleSystem::get_timestep_normal(MParticle p)
{
    // Calculating a_{1,i}^{(2)} = a_{0,i}^{(2)} + dt * a_{0,i}^{(3)}
    double ax1_2 = p.a2.x + p.dt * p.a3.x;
    double ay1_2 = p.a2.y + p.dt * p.a3.y;
    double az1_2 = p.a2.z + p.dt * p.a3.z;

    // |a_{1,i}|
    double abs_a1 = nu->get_magnitude(p.f.a[0], p.f.a[1], p.f.a[2]);
    // |j_{1,i}|
    double abs_j1 = nu->get_magnitude(p.f.a1[0], p.f.a1[1], p.f.a1[2]);
    // |j_{1,i}|^{2}
    double abs_j12  = abs_j1 * abs_j1;
    // a_{1,i}^{(3)} = a_{0,i}^{(3)} because the 3rd-order interpolation
    double abs_a1_3 = nu->get_magnitude(p.a3.x, p.a3.y, p.a3.z);
    // |a_{1,i}^{(2)}|
    double abs_a1_2 = nu->get_magnitude(ax1_2, ay1_2, az1_2);
    // |a_{1,i}^{(2)}|^{2}
    double abs_a1_22  = abs_a1_2 * abs_a1_2;

    double normal_dt = sqrt(ETA_B * ((abs_a1 * abs_a1_2 + abs_j12) / (abs_j1 * abs_a1_3 + abs_a1_22)));
    return normal_dt;
}

void MultipleSystem::init_timestep()
{
    double dt_min = 1e6;
    if (members < 2)
    {
        std::cout << "[Error] We cannot initialize the timesteps of a system with only 1 member" << std::endl;
    }
    else
    {
        for (int i = 0; i < members; i++)
        {
            // Accesing system particle
            MParticle p = parts[i];

            // |a|
            double a = sqrt(p.f.a[0] * p.f.a[0] +
                            p.f.a[1] * p.f.a[1] +
                            p.f.a[2] * p.f.a[2]);

            // |a1|
            double a1 = sqrt(p.f.a1[0] * p.f.a1[0] +
                             p.f.a1[1] * p.f.a1[1] +
                             p.f.a1[2] * p.f.a1[2]);

            // `dt` for every member [S. Konstantinidis 2010] eq. 22
            double dt_0 = a/a1;

            // Finding the minimum
            if (dt_0 < dt_min)
                dt_min = dt_0;

            parts[i].dt = 0.5 * D_TIME_MIN;
        }

        // ETA_B for the binary system
        ETA_B = D_TIME_MIN / (2.0 * dt_min);
    }
}

void MultipleSystem::get_orbital_elements()
{
    // Calculate orbital elements only for binary systems
    if (members == 2)
    {
        MParticle p = parts[0];
        MParticle q = parts[1];

        //double mu = (q.r.w * p.r.w) / (q.r.w + p.r.w);

        double rx = q.r.x - p.r.x;
        double ry = q.r.y - p.r.y;
        double rz = q.r.z - p.r.z;

        double vx = q.v.x - p.v.x;
        double vy = q.v.y - p.v.y;
        double vz = q.v.z - p.v.z;

        double jx = ry * vz - rz * vy;
        double jy = rz * vx - rx * vz;
        double jz = rx * vy - ry * vx;

        double rr = sqrt(rx*rx + ry*ry + rz*rz);
        double vv = sqrt(vx*vx + vy*vy + vz*vz);
        double jj = sqrt(jx*jx + jy*jy + jz*jz);

        double j2 = jj*jj;
        double v2 = vv*vv;
        //double m1m2 = q.r.w * p.r.w;

        double mu_std = G * (q.r.w + p.r.w);
        double espec = v2 * 0.5 - mu_std/rr;
        double semimajor_ini = -mu_std / (2*espec);
        double ecc_ini = sqrt(1.0+2.0*espec*j2/(mu_std*mu_std));

        std::cout << "mu_std: " << mu_std << std::endl;
        std::cout << "rr: " << rr << std::endl;
        std::cout << "vv: " << vv << std::endl;
        std::cout << "jj: " << jj << std::endl;
        std::cout << "a: " << semimajor_ini << std::endl;
        std::cout << "espec: " << espec << std::endl;
        std::cout << "ecc: " << ecc_ini << std::endl;
    }
}

double MultipleSystem::get_energy()
{
    double pot = 0.0;
    double kin = 0.0;

    for (int i = 0; i < members; i++)
    {
        MParticle pi = parts[i];
        double epot_tmp = 0.0;
        for (int j = i+1; j < members; j++)
        {

            MParticle pj = parts[j];
            double rx = pj.r.x - pi.r.x;
            double ry = pj.r.y - pi.r.y;
            double rz = pj.r.z - pi.r.z;
            double r2 = rx*rx + ry*ry + rz*rz;

            epot_tmp -= (pj.r.w * pi.r.w) / sqrt(r2);
        }

        double vx = pi.v.x * pi.v.x;
        double vy = pi.v.y * pi.v.y;
        double vz = pi.v.z * pi.v.z;
        double v2 = vx + vy + vz;

        double ekin_tmp = 0.5 * pi.r.w * v2;

        kin += ekin_tmp;
        pot += epot_tmp;
    }

    //printf("BB: K = %.15e | U = %.15e | Tot = %.15e\n", kin, pot, kin+pot);
    return kin + pot;
}

void MultipleSystem::update_timestep(double CTIME)
{
    for (int i = 0; i < members; i++)
    {
        // Saving old timestep
        MParticle part = parts[i];

        part.old_dt = part.dt;
        // Getting new timestep
        double first_dt  = get_timestep_normal(part);
        double new_dt = 0.0;

        // Check even or odd timestep
        if ((int)(CTIME / part.dt) % 2 == 0)
        {
            //std::cout << "Even" << std::endl;
            // We double it
            part.dt = part.old_dt * 2;
            evaluation(NULL);
            correction(CTIME, false);
            double last_dt  = get_timestep_normal(part);
            double avg = 0.5 * (first_dt + last_dt);

            if ( (2 * part.old_dt <= avg ||
                     part.old_dt <= avg) &&
                 2 * part.old_dt >= avg )
            {
                //std::cout << "Doubling it" << std::endl;
                new_dt = part.dt;
            }
            else // If not, we keep the same value and try again
            {
                part.dt = part.old_dt;
                evaluation(NULL);
                correction(CTIME, false);
                double last_dt  = get_timestep_normal(part);
                double avg = 0.5 * (first_dt + last_dt);

                if ( (2 * part.old_dt <= avg ||
                         part.old_dt <= avg) &&
                     2 * part.old_dt >= avg )
                {
                    //std::cout << "Keeping the same" << std::endl;
                    new_dt = part.dt;
                }
                else
                {
                    if ( part.old_dt * 0.5 <= avg &&
                         part.old_dt > avg)
                    {
                        //std::cout << "Halving it" << std::endl;
                        new_dt = part.dt * 0.5;
                    }
                    else
                    {
                        //std::cout << "...in any other case we halve" << std::endl;
                        new_dt = part.dt * 0.5;
                    }
                }
            }
        }
        else
        {
            evaluation(NULL);
            correction(CTIME, false);
            double last_dt  = get_timestep_normal(part);
            double avg = 0.5 * (first_dt + last_dt);

            if ( (2 * part.old_dt <= avg ||
                     part.old_dt <= avg) &&
                 2 * part.old_dt >= avg )
            {
                //std::cout << "Keeping the same" << std::endl;
                new_dt = part.dt;
            }
            else
            {
                if ( part.old_dt * 0.5 <= avg &&
                     part.old_dt > avg )
                {
                    //std::cout << "Halving it" << std::endl;
                    new_dt = part.dt * 0.5;
                }
                else
                {
                    //std::cout << "...in any other case we halve" << std::endl;
                    new_dt = part.dt * 0.5;
                }
            }
        }

        if (new_dt > D_TIME_MAX)
            new_dt = D_TIME_MAX;

        if (new_dt < D_TIME_MIN)
            new_dt = D_TIME_MIN;

        part.dt = new_dt;
        part.t = CTIME;
        part.r.x = part.p.r[0];
        part.r.y = part.p.r[1];
        part.r.z = part.p.r[2];

        // Correcting velocity
        part.v.x = part.p.v[0];
        part.v.y = part.p.v[1];
        part.v.z = part.p.v[2];

        // Writing information
        parts[i] = part;
    }
}
