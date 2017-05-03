#define DOUBLE double

int N, print_energy;
DOUBLE M, G, EPSILON, TOLERANCE, energy;
DOUBLE *pos_x, *pos_y, *pos_z, *speed_x, *speed_y, *speed_z;
DOUBLE *acc_x, *acc_y, *acc_z;

void allocate_space(void);
int *where(int *n, int *pos, DOUBLE *min_bound, DOUBLE *max_bound);
void init_from_ram(DOUBLE *x, DOUBLE *y, DOUBLE *z, DOUBLE *vx, DOUBLE *vy, DOUBLE *vz);
void init_conditions(int n, DOUBLE m, DOUBLE g, DOUBLE epsilon, DOUBLE tolerance, int energy_print, int threads);
