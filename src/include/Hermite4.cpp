#include "Hermite4.hpp"

Hermite4::Hermite4(NbodySystem *ns, Logger *logger, NbodyUtils *nu)
{
    this->ns     = ns;
    this->logger = logger;
    this->nu     = nu;
    alloc_arrays_host();
    init_data();
}

Hermite4::~Hermite4()
{
    free_arrays_host();
}

int Hermite4::find_particles_to_move(double ITIME)
{

    int j = 0;
    for (int i = FIRST_PARTICLE; i < ns->n; i++)
    {
        ns->h_move[i] = -1;
        // Getting larger mass of the system
        if (ns->h_r[i].w > ns->max_mass)
            ns->max_mass = ns->h_r[i].w;

        double tmp_time = ns->h_t[i] + ns->h_dt[i];
        if(std::fabs(ITIME - tmp_time) < 2*std::numeric_limits<double>::epsilon())
        {
            ns->h_move[j] = i;
            j++;

            ns->h_f[i].nb = 0;
        }
    }
    return j;
}

void Hermite4::next_integration_time(double &CTIME)
{
    // Initial number as the maximum
    CTIME = ns->h_t[FIRST_PARTICLE] + ns->h_dt[FIRST_PARTICLE];
    for (int i = FIRST_PARTICLE + 1; i < ns->n; i++)
    {
        double time = ns->h_t[i] + ns->h_dt[i];
        if(time < CTIME)
        {
            CTIME = time;
        }
    }
}

void Hermite4::init_dt(double &CTIME, float ETA, double ITIME)
{
    // Aarseth initial timestep
    // dt_{i} = ETA_S * sqrt( (|a|) / (|j|) )
    double tmp_dt;
    for (int i = FIRST_PARTICLE; i < ns->n; i++)
    {
        double a2 = nu->get_magnitude(ns->h_f[i].a[0],
                                      ns->h_f[i].a[1],
                                      ns->h_f[i].a[2]);
        double j2 = nu->get_magnitude(ns->h_f[i].a1[0],
                                      ns->h_f[i].a1[1],
                                      ns->h_f[i].a1[2]);
        tmp_dt = ETA * (a2/j2);

        // Adjusting to block timesteps
        // to the nearest-lower power of two
        int exp = (int)(std::ceil(log(tmp_dt)/log(2.0))-1);
        tmp_dt = pow(2,exp);
        //2 << (exp - 1);

        if (tmp_dt < D_TIME_MIN)
            tmp_dt = D_TIME_MIN;
        else if (tmp_dt > D_TIME_MAX)
            tmp_dt = D_TIME_MAX;

        ns->h_dt[i] = tmp_dt;
        //ns->h_dt[i] = D_TIME_MIN;
        //ns->h_t[i] = 0.0;
        ns->h_t[i] = ITIME;

        // Obtaining the first integration time
        if(tmp_dt < CTIME)
            CTIME = tmp_dt;
    }
}

void Hermite4::save_old_acc_jrk(int nact)
{
    //#pragma omp parallel for
    for (int k = 0; k < nact; k++)
    {
        int i = ns->h_move[k];
        ns->h_old[i].a[0]  = ns->h_f[i].a[0];
        ns->h_old[i].a[1]  = ns->h_f[i].a[1];
        ns->h_old[i].a[2]  = ns->h_f[i].a[2];

        ns->h_old[i].a1[0] = ns->h_f[i].a1[0];
        ns->h_old[i].a1[1] = ns->h_f[i].a1[1];
        ns->h_old[i].a1[2] = ns->h_f[i].a1[2];
    }

}

void Hermite4::alloc_arrays_host()
{
    //int d4_size = ns->n * sizeof(double4);
    //int d1_size = ns->n * sizeof(double);
    //int i1_size = ns->n * sizeof(int);
    //int ff_size = ns->n * sizeof(Forces);
    //int pp_size = ns->n * sizeof(Predictor);

    ns->h_f        = new Forces[ns->n];
    ns->h_a2       = new double3[ns->n];
    ns->h_a3       = new double3[ns->n];
    ns->h_old      = new Forces[ns->n];
    ns->h_t        = new double[ns->n];
    ns->h_dt       = new double[ns->n];
    ns->h_move     = new int[ns->n];
    ns->h_p        = new Predictor[ns->n];
    ns->h_i        = new Predictor[ns->n];
    ns->h_ekin     = new double[ns->n];
    ns->h_epot     = new double[ns->n];
    ns->h_r_sphere = new double[ns->n];

}

void Hermite4::free_arrays_host()
{
    delete ns->h_f;
    delete ns->h_a2;
    delete ns->h_a3;
    delete ns->h_old;
    delete ns->h_t;
    delete ns->h_dt;
    delete ns->h_move;
    delete ns->h_p;
    delete ns->h_i;
    delete ns->h_ekin;
    delete ns->h_epot;
    delete ns->h_r_sphere;
}

void Hermite4::init_data()
{
    double3 empty = {0.0, 0.0, 0.0};
    ns->m_g = 0.0;
    for (int i = 0; i < ns->n; i++) {

        double mass = ns->h_r[i].w;

        ns->h_p[i].m    = mass;

        ns->h_p[i].r[0] = ns->h_r[i].x;
        ns->h_p[i].r[1] = ns->h_r[i].y;
        ns->h_p[i].r[2] = ns->h_r[i].z;

        ns->h_p[i].v[0] = ns->h_v[i].x;
        ns->h_p[i].v[1] = ns->h_v[i].y;
        ns->h_p[i].v[2] = ns->h_v[i].z;

        ns->h_f[i].a[0]  = 0.0;
        ns->h_f[i].a[1]  = 0.0;
        ns->h_f[i].a[2]  = 0.0;

        ns->h_f[i].a1[0] = 0.0;
        ns->h_f[i].a1[1] = 0.0;
        ns->h_f[i].a1[2] = 0.0;

        ns->h_f[i].nb = 0;

        ns->h_a2[i]      = empty;
        ns->h_a3[i]      = empty;

        ns->h_old[i].a[0]  = 0.0;
        ns->h_old[i].a[1]  = 0.0;
        ns->h_old[i].a[2]  = 0.0;

        ns->h_old[i].a1[0] = 0.0;
        ns->h_old[i].a1[1] = 0.0;
        ns->h_old[i].a1[2] = 0.0;

        ns->h_t[i]       = 0.0;
        ns->h_dt[i]      = 0.0;

        ns->h_move[i]    = 0;

        ns->h_r_sphere[i] = 0.0;

        // Heaviest star
        if (mass > ns->m_g)
            ns->m_g = mass;

    }
}
