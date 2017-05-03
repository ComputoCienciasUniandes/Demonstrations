#include "init.h"
#include "omp.h"
#include <stdio.h>
#include <stdlib.h>

void init_conditions(int n, DOUBLE m, DOUBLE g, DOUBLE epsilon, DOUBLE tolerance, int energy_print, int threads)
{
    N = n;
    M = m;
    G = g;
    EPSILON = epsilon;
    TOLERANCE = tolerance;
    print_energy = energy_print;
    if (threads >= 1)
    {
        omp_set_num_threads(threads);
    }
}

void allocate_space(void)
{
    pos_x = malloc(N*sizeof(DOUBLE));
    pos_y = malloc(N*sizeof(DOUBLE));
    pos_z = malloc(N*sizeof(DOUBLE));
    speed_x = malloc(N*sizeof(DOUBLE));
    speed_y = malloc(N*sizeof(DOUBLE));
    speed_z = malloc(N*sizeof(DOUBLE));
    acc_x = malloc(N*sizeof(DOUBLE));
    acc_y = malloc(N*sizeof(DOUBLE));
    acc_z = malloc(N*sizeof(DOUBLE));
}

void init_from_ram(DOUBLE *x, DOUBLE *y, DOUBLE *z, DOUBLE *vx, DOUBLE *vy, DOUBLE *vz)
{
    allocate_space();
    int i = 0;
    for(i = 0; i < N; i++)
    {
        pos_x[i] = x[i];
        pos_y[i] = y[i];
        pos_z[i] = z[i];
        speed_x[i] = vx[i];
        speed_y[i] = vy[i];
        speed_z[i] = vz[i];
    }
}

int *where(int *n, int *pos, DOUBLE *min_bound, DOUBLE *max_bound)
{
    int i = 0, j = 0;
    DOUBLE x, y, z;
    int *the_ones = malloc(*n*sizeof(int));

    for(i = 0; i<*n; i++)
    {
        x = pos_x[pos[i]];
        y = pos_y[pos[i]];
        z = pos_z[pos[i]];
        if((x >= min_bound[0]) && (x < max_bound[0]))
        {
            if((y >= min_bound[1]) && (y < max_bound[1]))
            {
                if((z >= min_bound[2]) && (z < max_bound[2]))
                {
                    the_ones[j] = pos[i];
                    j++;
                }
            }
        }
    }
    the_ones = realloc(the_ones, j*sizeof(int));
    *n = j;
    return the_ones;
}
