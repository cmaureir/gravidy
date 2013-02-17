#include "extra_utils.hpp"

void print_all(int limit, float ITIME)
{
//    for (int i = 0; i < limit; i++) {
//        printf("%.10f %6d %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
//                ITIME,
//                i,
//                h_r[i].x, h_r[i].y, h_r[i].z,
//                h_v[i].x, h_v[i].y, h_v[i].z,
//                h_a[i].x, h_a[i].y, h_a[i].z,
//                h_a1[i].x, h_a1[i].y, h_a1[i].z,
//                h_dt[i]);
//    }
}

void print_positions(int limit)
{
    for (int i = 0; i < limit; i++) {
        printf("%6d %.10f %.10f %.10f\n", i, h_r[i].x, h_r[i].y, h_r[i].z );
    }
}

void print_velocities(int limit)
{
    for (int i = 0; i < limit; i++) {
        printf("%6d %.10f %.10f %.10f\n", i, h_v[i].x, h_v[i].y, h_v[i].z );
    }
}
void print_accelerations(int limit)
{
//    for (int i = 0; i < limit; i++) {
//        printf("%6d %.10f %.10f %.10f\n", i, h_a[i].x, h_a[i].y, h_a[i].z );
//    }

}
void print_jrks(int limit)
{
//    for (int i = 0; i < limit; i++) {
//        printf("%6d %.10f %.10f %.10f\n", i, h_a1[i].x, h_a1[i].y, h_a1[i].z );
//    }
}

void print_accjrk(int limit)
{
//    for (int i = 0; i < limit; i++) {
//        printf("%6d %4.10f %4.10f %4.10f %4.10f %4.10f %4.10f\n",
//                i,
//                h_a[i].x, h_a[i].y, h_a[i].z,
//                h_a1[i].x, h_a1[i].y, h_a1[i].z
//                );
//    }
}
void print_masses(int limit)
{
    for (int i = 0; i < limit; i++) {
        printf("%6d %.10f\n", i, h_m[i] );
    }
}

void print_times(int limit, float itime)
{
    int exp = 0;
    for (int i = 0; i < limit; i++) {
        exp = (int)std::ceil(log(h_dt[i])/log(2));
        printf("%.10f %6d %.10f %2d %.10f\n", itime, i, h_dt[i], exp, h_t[i]);
    }
}

// Print old
void print_old(int limit)
{
    for (int i = 0; i < limit; i++) {
        printf("%6d %.10f %.10f %.10f %.10f %.10f %.10f\n",
                i,
                h_old_a[i].x, h_old_a[i].y, h_old_a[i].z,
                h_old_a1[i].x, h_old_a1[i].y, h_old_a1[i].z);
    }
}

void print_predicted(int limit)
{
//    for (int i = 0; i < limit; i++) {
//        printf("%6d %.10f %.10f %.10f %.10f %.10f %.10f\n",
//                i,
//                h_p_r[i].x, h_p_r[i].y, h_p_r[i].z,
//                h_p_v[i].x, h_p_v[i].y, h_p_v[i].z);
//    }
}

void print_movement(int limit, int total, float ITIME)
{
    printf("%.6f ", ITIME);
    for (int i = 0; i < limit; i++)
    {
        int value = 0;
        for(int j = 0; j < total ; j++)
        {
            int k = h_move[j];
            if (k == i)
                value = 1;

        }
        printf("%d ", value);
    }
    printf("\n");
}


void print_particle(int i)
{
//    printf("%5d %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
//            i,
//            h_r[i].x, h_r[i].y, h_r[i].z, h_v[i].x, h_v[i].y, h_v[i].z,
//            h_a[i].x, h_a[i].y, h_a[i].z, h_a1[i].x, h_a1[i].y, h_a1[i].z,
//            h_dt[i]);

}



double get_magnitude(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}
