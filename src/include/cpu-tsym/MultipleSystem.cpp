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

    return pcm;
}

// A new particle will be created always with two members
MParticle MultipleSystem::get_new_particle(MParticle p1, MParticle p2)
{

    MParticle newp;

    newp.r.w  = p1.r.w + p2.r.w;

    //newp.r.x = p2.r.x - p1.r.x;
    //newp.r.y = p2.r.y - p1.r.y;
    //newp.r.z = p2.r.z - p1.r.z;

    //double r = sqrt(newp.r.x * newp.r.x +
    //                newp.r.y * newp.r.y +
    //                newp.r.z * newp.r.z);

    ////double rinv = 1.0/r;
    ////double r3inv = rinv*rinv*rinv;

    //newp.v.x = p2.v.x - p1.v.x;
    //newp.v.y = p2.v.y - p1.v.y;
    //newp.v.z = p2.v.z - p1.v.z;

    ////newp.a.x = -(newp.m * newp.r.x)*r3inv + (p2.f.a[0] - p1.f.a[0]);
    ////newp.a.y = -(newp.m * newp.r.y)*r3inv + (p2.f.a[1] - p1.f.a[1]);
    ////newp.a.z = -(newp.m * newp.r.z)*r3inv + (p2.f.a[2] - p1.f.a[2]);

    //// TODO: Which a1, a2, a3, t, dt, we use?

    return newp;
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

    //double r2   = rx*rx + ry*ry + rz*rz + ns->e2;
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
    int i;
    for (i = 0; i < members; i++)
    {
        parts[i].old_f = parts[i].f;
    }
}

void MultipleSystem::evaluation()
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

        // Write on the Forces array
        parts[i].f = ff;
    }
}

void MultipleSystem::correction(double CTIME)
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

        // Correcting position
        parts[i].p.r[0] = p0.r[0] + dt4c * a2i.x + dt5c * a3i.x;
        parts[i].p.r[1] = p0.r[1] + dt4c * a2i.y + dt5c * a3i.y;
        parts[i].p.r[2] = p0.r[2] + dt4c * a2i.z + dt5c * a3i.z;

        // Correcting velocity
        parts[i].p.v[0] = p0.v[0] + dt3c * a2i.x + dt4c * a3i.x;
        parts[i].p.v[1] = p0.v[1] + dt3c * a2i.y + dt4c * a3i.y;
        parts[i].p.v[2] = p0.v[2] + dt3c * a2i.z + dt4c * a3i.z;

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

void MultipleSystem::init_timestep(double &CTIME)
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

            //parts[i].dt = 0.0009765625;
            parts[i].dt = 0.5 * D_TIME_MIN;
            //parts[i].t = 0.0;
        }

        // ETA_B for the binary system
        ETA_B = D_TIME_MIN / (2.0 * dt_min);
        //ETA_B = 0.01;
        //CTIME = parts[0].t + 0.5 * D_TIME_MIN;
        CTIME = parts[0].t + D_TIME_MIN;

        // First timestep is the minimum
        //std::cout << std::scientific;
        //std::cout << "Initial dt " << 0.5 * D_TIME_MIN << std::endl;
    }
}

void MultipleSystem::update_information(double c_time)
{

    for (int i = 0; i < members; i++)
    {
        MParticle part = parts[i];
        // New time
        part.t = c_time;

        // Updating the real pos/vel to the last corrected one.
        parts[i].r.x = part.p.r[0];
        parts[i].r.y = part.p.r[1];
        parts[i].r.z = part.p.r[2];

        parts[i].v.x = part.p.v[0];
        parts[i].v.y = part.p.v[1];
        parts[i].v.z = part.p.v[2];

        double normal_dt  = get_timestep_normal(part);
        if (normal_dt < part.dt)
        {
            if (0.5 * part.dt < D_TIME_MIN)
                part.dt = D_TIME_MIN;
            else
                part.dt *= 0.5;
        }
        else if (normal_dt > 2*part.dt && (fmod(part.t, 2*part.dt) == 0))
        {
            if (2*part.dt < D_TIME_MAX)
                part.dt *= 2;
            else
                part.dt = D_TIME_MAX;
        }
        parts[i].dt = part.dt;
        parts[i].t = part.t;
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
    return kin + pot;
}

//void MultipleSystem::update_timestep(double CTIME)
//{
//    for (int i = 0; i < members; i++)
//    {
//        MParticle p = parts[i];
//
//        // |a|
//        double a = sqrt(p.f.a[0] * p.f.a[0] +
//                        p.f.a[1] * p.f.a[1] +
//                        p.f.a[2] * p.f.a[2]);
//
//        // |a1|
//        double a1 = sqrt(p.f.a1[0] * p.f.a1[0] +
//                         p.f.a1[1] * p.f.a1[1] +
//                         p.f.a1[2] * p.f.a1[2]);
//
//        // `dt` for every member [S. Konstantinidis 2010] eq. 22
//        double dt_b = ETA_B * a/a1;
//        double new_dt;
//        double dt_avg;
//
//        // Check *even* or *odd* with the previous dt
//
//        //std::cout << "t: " << CTIME << ", dt: " << parts[i].dt << std::endl;
//        if (fmod((CTIME/parts[i].dt),2) == 0)
//        {
//            // Even
//
//            // We double the old `dt` value
//            new_dt = 2.0 * parts[i].dt;
//
//            // `dt` average
//            dt_avg = 0.5 * (new_dt + dt_b);
//
//            if (new_dt  < dt_avg && 0.5 * dt_avg < new_dt)
//            {
//                // We use twice the value
//                new_dt *= 2.0;
//                //std::cout << "even, mutiplied by 2" << std::endl;
//            }
//            else
//            {
//                new_dt = new_dt;
//                //std::cout << "even, more steps" << std::endl;
//                if (new_dt  < dt_avg && 0.5 * dt_avg < new_dt)
//                {
//                    // We use the same
//                    new_dt = new_dt;
//                    //std::cout << "even, we use the same" << std::endl;
//                }
//                else
//                {
//                    // new_dt/2.0
//                    new_dt *= 0.5;
//                    //std::cout << "even, divided by 2" << std::endl;
//                }
//            }
//        }
//        else
//        {
//            // Odd
//
//            // We keep the old `dt` value
//            new_dt = parts[i].dt;
//
//            // `dt` average
//            dt_avg = 0.5 * (new_dt + dt_b);
//
//            if (new_dt  < dt_avg && 0.5 * dt_avg < new_dt)
//            {
//                // We use the same
//                new_dt = new_dt;
//                //std::cout << "odd, we use the same" << std::endl;
//            }
//            else
//            {
//                // new_dt/2.0
//                new_dt *= 0.5;
//                //std::cout << "odd, divided by 2" << std::endl;
//            }
//        }
//        parts[i].dt = new_dt;
//
//        //std::cout << i << " dt_b=" << dt_b << " dt=" << parts[i].dt << std::endl;
//    }
//}
