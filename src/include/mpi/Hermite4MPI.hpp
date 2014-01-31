#ifndef HERMITE4MPI_HPP
#define HERMITE4MPI_HPP
#include "../Hermite4.hpp"

#define MPI_NUM_SLAVES 32

class Hermite4MPI : public Hermite4 {
    public:
        Hermite4MPI(NbodySystem *ns, Logger *logger, NbodyUtils *nu,
                    int rank, int nproc);
        ~Hermite4MPI();

        int rank;
        int nproc;
        int tag;
        MPI_Status   status;
        MPI_Datatype f_type;
        MPI_Op       f_op;

        int chunk_size;
        int chunk_begin;
        int chunk_end;
        Forces *tmp_f;

        void alloc_slaves_memory(int rank);
        void force_calculation(Predictor pi, Predictor pj, Forces &fi);
        void init_acc_jrk();
        void update_acc_jrk(int nact);
        void predicted_pos_vel(double ITIME);
        void correction_pos_vel(double ITIME, int nact);
        void integration();
};

#endif
