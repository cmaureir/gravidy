#ifndef HERMITE4_HPP
#define HERMITE4_HPP
#include "utils/Logger.hpp"
#include "utils/NbodyUtils.hpp"
#include <limits>

class Hermite4 {
    public:
        Hermite4(NbodySystem *ns, Logger *logger, NbodyUtils *nu);
        ~Hermite4();

        NbodySystem *ns;
        Logger      *logger;
        NbodyUtils  *nu;

        int  find_particles_to_move(double ITIME);
        void next_integration_time(double &ATIME);
        void init_dt(double &ATIME, float ETA);
        void save_old_acc_jrk(int nact);
        void alloc_arrays_host();
        void free_arrays_host();
        void init_data();

        /** Virtual methods to be implemented by the different versions **/

        /// @brief This method contains the whole integration process,
        ///        using the predictor-corrector scheme
        virtual void integration() {}

        /// @brief The prediction step method, which extrapolates the
        ///        positions and velocities of ALL the particles.
        /// @param itime Current integration time
        /// @param t Array containing the particle Times.
        /// @param r Array containing the particle Positions.
        /// @param v Array containing the particle Velocities.
        /// @param f Array containing the particle Forces.
        /// @param p Array containing the particle Predicted position and velocity.
        virtual void predicted_pos_vel(double itime,
                                       double *t,
                                       double4 *r,
                                       double4 *v,
                                       Forces *f,
                                       Predictor *p) {}

        /// @brief The prediction step method, which extrapolates the
        ///        positions and velocities of ALL the particles.
        /// @param itime Current integration time
        /// @param nact Amount of active particles
        /// @param dt Array containing the particle Timesteps.
        /// @param t Array containing the particle Times.
        /// @param move Array containing the ID for the active particles.
        /// @param p Array containing the particle Predicted position and velocity.
        /// @param f Array containing the particle Forces.
        /// @param old Array containing the particle previous Forces.
        /// @param a2 Array containing the particle Acceleration 2nd derivative.
        /// @param a3 Array containing the particle Acceleration 3rd derivative.
        /// @param r Array containing the particle Positions.
        /// @param v Array containing the particle Velocities.
        virtual void correction_pos_vel(double itime,
                                        int nact,
                                        double *dt,
                                        double *t,
                                        int *move,
                                        Predictor *p,
                                        Forces *f,
                                        Forces *old,
                                        double3 *a2,
                                        double3 *a3,
                                        double4 *r,
                                        double4 *v) {}

        virtual void force_calculation(Predictor pi, Predictor pj, Forces &fi) {}
        virtual void init_acc_jrk() {}
        virtual void update_acc_jrk(int nact) {}

};

#endif // HERMITE4_HPP
