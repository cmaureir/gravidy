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
    for (int i = 0; i < ns->n; i++)
    {
        ns->h_move[i] = -1;
        double tmp_time = ns->h_t[i] + ns->h_dt[i];
        if(std::fabs(ITIME - tmp_time) < 2*std::numeric_limits<double>::epsilon())
        {
            ns->h_move[j] = i;
            j++;
        }
    }
    return j;
}

void Hermite4::next_integration_time(double &ATIME)
{
    // Big number to find the minimum
    //ATIME = 1.0e10;
    ATIME = ns->h_t[0] + ns->h_dt[0];
    for (int i = 1; i < ns->n; i++)
    {
        double time = ns->h_t[i] + ns->h_dt[i];
        if(time < ATIME)
        {
            ATIME = time;
        }
    }
}

void Hermite4::init_dt(double &ATIME)
{
    // Aarseth initial timestep
    // dt_{i} = ETA_S * sqrt( (|a|) / (|j|) )
    double tmp_dt;
    for (int i = 0; i < ns->n; i++)
    {
        double a2 = nu->get_magnitude(ns->h_f[i].a[0],
                                      ns->h_f[i].a[1],
                                      ns->h_f[i].a[2]);
        double j2 = nu->get_magnitude(ns->h_f[i].a1[0],
                                      ns->h_f[i].a1[1],
                                      ns->h_f[i].a1[2]);
        tmp_dt = ETA_S * (a2/j2);

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
        ns->h_t[i] = 0.0;

        // Obtaining the first integration time
        if(tmp_dt < ATIME)
            ATIME = tmp_dt;
    }
}

void Hermite4::save_old_acc_jrk(int nact)
{
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
    int d4_size = ns->n * sizeof(double4);
    int d1_size = ns->n * sizeof(double);
    int i1_size = ns->n * sizeof(int);
    int ff_size = ns->n * sizeof(Forces);
    int pp_size = ns->n * sizeof(Predictor);

    ns->h_f    = new Forces[ff_size];
    ns->h_a2   = new double4[d4_size];
    ns->h_a3   = new double4[d4_size];
    ns->h_old  = new Forces[ff_size];
    ns->h_t    = new double[d1_size];
    ns->h_dt   = new double[d1_size];
    ns->h_move = new int[i1_size];
    ns->h_p    = new Predictor[pp_size];
    ns->h_i    = new Predictor[pp_size];
    ns->h_ekin = new double[d1_size];
    ns->h_epot = new double[d1_size];

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
}

void Hermite4::init_data()
{
    double4 empty = {0.0, 0.0, 0.0, 0.0};
    for (int i = 0; i < ns->n; i++) {

        ns->h_p[i].m    = ns->h_r[i].w;

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
    }
}
