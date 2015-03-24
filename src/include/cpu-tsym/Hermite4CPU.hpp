#ifndef HERMITE4CPU_HPP
#define HERMITE4CPU_HPP
#include "../Hermite4.hpp"
#include "MultipleSystem.hpp"

typedef struct
{
    int id_a;
    int id_b;
} binary_id;

class Hermite4CPU : public Hermite4 {
    public:
        Hermite4CPU(NbodySystem *ns, Logger *logger, NbodyUtils *nu);
        ~Hermite4CPU();

        Forces *h_fn;
        Forces *h_fn_old;
        double *h_tn;
        int **nb_list;
        int nb_number;
        int *ghosts;

        void force_calculation(Predictor pi, Predictor pj, Forces &fi,
                               int i, int j, double hi);

        void init_acc_jrk(Predictor *p, Forces *f, double *r_sphere);

        void update_acc_jrk(int nact, int *move, double *r_sphere,
                            Predictor *p, Forces *f);

        void predicted_pos_vel(double itime, double *t, double4 *r, double4 *v,
                               Forces *f, Predictor *p);

        void correction_pos_vel(double itime, int nact, double *dt, double *t,
                                int *move, Predictor *p, Forces *f, Forces *old,
                                double3 *a2, double3 *a3, double4 *r, double4 *v);

        void update_neighbour_radius();

        void print_nb(double itime, int **nb_list, Forces *f, int n, Predictor *p,
                      double *r_sphere);

        bool get_close_encounters(double itime, int **nb_list, Forces *f, int n,
                                  Predictor *p, double *r_sphere, std::vector<binary_id> &pairs, int nact);

        void integration();


        void create_ghost_particle(MultipleSystem ms);

        void multiple_systems_integration(std::vector<MultipleSystem> ms, double ITIME);
};
#endif
