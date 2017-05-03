typedef struct box_str
{
    DOUBLE *points;
    DOUBLE *lower_bound;
    DOUBLE *upper_bound;
    DOUBLE *center_of_mass;
    DOUBLE coordinate_size;
    DOUBLE box_half_size;
    DOUBLE mass_unit;
    DOUBLE mass;
    int number_of_subs;
    int number_of_points;

    struct box_str *subBoxes;
} box;

box *init_tree(void);
void clean_tree(box *tree);
void bounds(DOUBLE *low_bounds, DOUBLE *size);
DOUBLE *calc_center_of_mass(int n, int *pos);
box *create_box(int n, int *points,  DOUBLE *lb, DOUBLE cs);
