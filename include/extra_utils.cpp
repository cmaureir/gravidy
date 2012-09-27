#include "extra_utils.hpp"

void print_all(int limit)
{
    cout << setw(6)  << "id";
    cout << setw(15) << "rx";
    cout << setw(15) << "ry";
    cout << setw(15) << "rz";
    cout << setw(15) << "vx";
    cout << setw(15) << "vy";
    cout << setw(15) << "vz";
    cout << setw(15) << "ax";
    cout << setw(15) << "ay";
    cout << setw(15) << "az";
    cout << setw(15) << "jx";
    cout << setw(15) << "jy";
    cout << setw(15) << "jz" << endl;
    for (int i = 0; i < limit; i++) {
        printf("%6d %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f\n",
                i,
                h_r[i].x, h_r[i].y, h_r[i].z,
                h_v[i].x, h_v[i].y, h_v[i].z,
                h_a[i].x, h_a[i].y, h_a[i].z,
                h_j[i].x, h_j[i].y, h_j[i].z,
                h_dt[i]);
    }
}

void print_positions(int limit)
{
    for (int i = 0; i < limit; i++) {
        printf("%6d %.14f %.14f %.14f\n", i, h_r[i].x, h_r[i].y, h_r[i].z );
    }
}

void print_velocities(int limit)
{
    for (int i = 0; i < limit; i++) {
        printf("%6d %.14f %.14f %.14f\n", i, h_v[i].x, h_v[i].y, h_v[i].z );
    }
}
void print_accelerations(int limit)
{
    for (int i = 0; i < limit; i++) {
        printf("%6d %.14f %.14f %.14f\n", i, h_a[i].x, h_a[i].y, h_a[i].z );
    }

}
void print_jrks(int limit)
{
    cout << left;
    cout << setw(10)  << "id";
    cout << setw(22) << "jx";
    cout << setw(22) << "jy";
    cout << setw(22) << "jz" << endl;
    for (int i = 0; i < limit; i++) {
        cout << setw(10)  << i;
        cout << setw(22) << h_j[i].x;
        cout << setw(22) << h_j[i].y;
        cout << setw(22) << h_j[i].z << endl;
    }
}
void print_masses(int limit)
{
    cout << left;
    cout << setw(10)  << "id";
    cout << setw(22) << "m" << endl;
    for (int i = 0; i < limit; i++) {
        cout << setw(10)  << i;
        cout << setw(22) << h_m[i] << endl;
    }
}

void print_times(int limit)
{
    for (int i = 0; i < limit; i++) {
        printf("%6d %.15f %.15f %.15f\n", i, h_dt[i], h_t[i], h_dt[i] + h_t[i]);
    }
}

// Print old
void print_old(int limit)
{
    cout << left;
    cout << setw(10)  << "id";
    cout << setw(22) << "old_ax";
    cout << setw(22) << "old_ay";
    cout << setw(22) << "old_az";
    cout << setw(22) << "old_jx";
    cout << setw(22) << "old_jy";
    cout << setw(22) << "old_jz" << endl;
    for (int i = 0; i < limit; i++) {
        cout << setw(10)  << i;
        cout << setw(22) << h_old_a[i].x;
        cout << setw(22) << h_old_a[i].y;
        cout << setw(22) << h_old_a[i].z;
        cout << setw(22) << h_old_j[i].x;
        cout << setw(22) << h_old_j[i].y;
        cout << setw(22) << h_old_j[i].z << endl;
    }
}

void print_predicted(int limit)
{
    for (int i = 0; i < limit; i++) {
        printf("%6d %.14f %.14f %.14f %.14f %.14f %.14f\n",
                i,
                h_p_r[i].x, h_p_r[i].y, h_p_r[i].z,
                h_p_v[i].x, h_p_v[i].y, h_p_v[i].z);
    }
}

double magnitude(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}
