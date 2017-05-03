#include "init.h"
#include "box.h"
#include "omp.h"
#include <stdlib.h>

void clean_tree(box *tree)
{
    int j = 0;
    int n = tree->number_of_subs;
    box temp;
    #pragma omp parallel for private(temp)
    for(j = 0; j<n; j++)
    {
        temp = tree->subBoxes[j];
        clean_tree(&temp);
    }
    free(tree->points);
    free(tree->lower_bound);
    free(tree->upper_bound);
    free(tree->center_of_mass);
    if(n > 0)
    {
        free(tree->subBoxes);
    }
}

box *init_tree(void)
{
    int n = N, i;
    int *points = malloc(sizeof(int)*n);
    DOUBLE *low_bounds = malloc(sizeof(DOUBLE)*3);
    DOUBLE rank = 0;
    #pragma omp parallel for
    for(i=0; i<n; i++)
    {
        points[i] = i;
    }
    bounds(low_bounds, &rank);
    return create_box(n, points, low_bounds, rank);;
}

void bounds(DOUBLE *low_bounds, DOUBLE *size)
{
    int i;
    DOUBLE x_min, x_max, x;
    DOUBLE y_min, y_max, y;
    DOUBLE z_min, z_max, z;

    x_min = x_max = pos_x[0];
    y_min = y_max = pos_y[0];
    z_min = z_max = pos_z[0];
    DOUBLE *rank = malloc(3*sizeof(DOUBLE));
    #pragma omp parallel for private(x, y, z)
    for(i = 0; i<N; i++)
    {
        x = pos_x[i];
        y = pos_y[i];
        z = pos_z[i];
        if(x < x_min)
        {
            x_min = x;
        }
        else if(x > x_max)
        {
            x_max = x;
        }
        if(y < y_min)
        {
            y_min = y;
        }
        else if(y > y_max)
        {
            y_max = y;
        }
        if(z < z_min)
        {
            z_min = z;
        }
        else if(z > z_max)
        {
            z_max = z;
        }
    }

    low_bounds[0] = x_min;
    low_bounds[1] = y_min;
    low_bounds[2] = z_min;
    rank[0] = x_max - x_min;
    rank[1] = y_max - y_min;
    rank[2] = z_max - z_min;
    if ((rank[0] >= rank[1]) && (rank[0] >= rank[2]))
    {
        *size = rank[0];
    }
    else if ((rank[1] >= rank[0]) && (rank[1] >= rank[2]))
    {
        *size = rank[1];
    }
    else if ((rank[2] >= rank[0]) && (rank[2] >= rank[1]))
    {
        *size = rank[2];
    }
    free(rank);
}

DOUBLE *calc_center_of_mass(int n, int *pos)
{
    int i = 0;
    DOUBLE *cm = malloc(3*sizeof(DOUBLE));
    cm[0] = 0; cm[1] = 0; cm[2] = 0;
    for(i = 0; i<n; i++)
    {
        cm[0] += pos_x[pos[i]];
        cm[1] += pos_y[pos[i]];
        cm[2] += pos_z[pos[i]];
    }
    cm[0] *= 1.0/n;
    cm[1] *= 1.0/n;
    cm[2] *= 1.0/n;
    return cm;
}

box *create_box(int n, int *points, DOUBLE *lb, DOUBLE cs)
{
    box *current_box = malloc(sizeof(box));
    current_box-> points = malloc(n*sizeof(DOUBLE));
    current_box-> lower_bound = malloc(3*sizeof(DOUBLE));
    current_box-> upper_bound = malloc(3*sizeof(DOUBLE));
    current_box-> center_of_mass = malloc(3*sizeof(DOUBLE));
    current_box-> coordinate_size = cs;
    current_box-> box_half_size = 0.5*cs;
    current_box-> mass = M*n;
    current_box-> number_of_points = n;

    DOUBLE *cm = calc_center_of_mass(n, points);
    current_box-> center_of_mass[0] = cm[0];
    current_box-> center_of_mass[1] = cm[1];
    current_box-> center_of_mass[2] = cm[2];
    current_box-> lower_bound[0] = lb[0];
    current_box-> lower_bound[1] = lb[1];
    current_box-> lower_bound[2] = lb[2];
    current_box-> upper_bound[0] = lb[0] + cs;
    current_box-> upper_bound[1] = lb[1] + cs;
    current_box-> upper_bound[2] = lb[2] + cs;

    int j = 0, k;
    DOUBLE *min_bound = malloc(3*sizeof(DOUBLE));
    DOUBLE *max_bound = malloc(3*sizeof(DOUBLE));
    int *sub_pos;
    if(n > 1)
    {
        current_box-> subBoxes = malloc(8*sizeof(box));
        k = n;
        min_bound[0] = lb[0];
        min_bound[1] = lb[1];
        min_bound[2] = lb[2];
        max_bound[0] = lb[0] + current_box-> box_half_size;
        max_bound[1] = lb[1] + current_box-> box_half_size;
        max_bound[2] = lb[2] + current_box-> box_half_size;
        sub_pos = where(&k, points, min_bound, max_bound);
        if (k > 0)
        {
            current_box-> subBoxes[j] = *create_box(k, sub_pos, min_bound, current_box-> box_half_size);
            j ++;
        }
        k = n;
        min_bound[0] = lb[0] + current_box-> box_half_size;
        min_bound[1] = lb[1];
        min_bound[2] = lb[2];
        max_bound[0] = lb[0] + cs;
        max_bound[1] = lb[1] + current_box-> box_half_size;
        max_bound[2] = lb[2] + current_box-> box_half_size;
        sub_pos = where(&k, points, min_bound, max_bound);
        if (k > 0)
        {
            current_box-> subBoxes[j] = *create_box(k, sub_pos, min_bound, current_box-> box_half_size);
            j ++;
        }
        k = n;
        min_bound[0] = lb[0];
        min_bound[1] = lb[1] + current_box-> box_half_size;
        min_bound[2] = lb[2];
        max_bound[0] = lb[0] + current_box-> box_half_size;
        max_bound[1] = lb[1] + cs;
        max_bound[2] = lb[2] + current_box-> box_half_size;
        sub_pos = where(&k, points, min_bound, max_bound);
        if (k > 0)
        {
            current_box-> subBoxes[j] = *create_box(k, sub_pos, min_bound, current_box-> box_half_size);
            j ++;
        }
        k = n;
        min_bound[0] = lb[0];
        min_bound[1] = lb[1];
        min_bound[2] = lb[2] + current_box-> box_half_size;
        max_bound[0] = lb[0] + current_box-> box_half_size;
        max_bound[1] = lb[1] + current_box-> box_half_size;
        max_bound[2] = lb[2] + cs;
        sub_pos = where(&k, points, min_bound, max_bound);
        if (k > 0)
        {
            current_box-> subBoxes[j] = *create_box(k, sub_pos, min_bound, current_box-> box_half_size);
            j ++;
        }
        k = n;
        min_bound[0] = lb[0] + current_box-> box_half_size;
        min_bound[1] = lb[1] + current_box-> box_half_size;
        min_bound[2] = lb[2];
        max_bound[0] = lb[0] + cs;
        max_bound[1] = lb[1] + cs;
        max_bound[2] = lb[2] + current_box-> box_half_size;
        sub_pos = where(&k, points, min_bound, max_bound);
        if (k > 0)
        {
            current_box-> subBoxes[j] = *create_box(k, sub_pos, min_bound, current_box-> box_half_size);
            j ++;
        }
        k = n;
        min_bound[0] = lb[0] + current_box-> box_half_size;
        min_bound[1] = lb[1];
        min_bound[2] = lb[2] + current_box-> box_half_size;
        max_bound[0] = lb[0] + cs;
        max_bound[1] = lb[1] + current_box-> box_half_size;
        max_bound[2] = lb[2] + cs;
        sub_pos = where(&k, points, min_bound, max_bound);
        if (k > 0)
        {
            current_box-> subBoxes[j] = *create_box(k, sub_pos, min_bound, current_box-> box_half_size);
            j ++;
        }
        k = n;
        min_bound[0] = lb[0];
        min_bound[1] = lb[1] + current_box-> box_half_size;
        min_bound[2] = lb[2] + current_box-> box_half_size;
        max_bound[0] = lb[0] + current_box-> box_half_size;
        max_bound[1] = lb[1] + cs;
        max_bound[2] = lb[2] + cs;
        sub_pos = where(&k, points, min_bound, max_bound);
        if (k > 0)
        {
            current_box-> subBoxes[j] = *create_box(k, sub_pos, min_bound, current_box-> box_half_size);
            j ++;
        }
        k = n;
        min_bound[0] = lb[0] + current_box-> box_half_size;
        min_bound[1] = lb[1] + current_box-> box_half_size;
        min_bound[2] = lb[2] + current_box-> box_half_size;
        max_bound[0] = lb[0] + cs;
        max_bound[1] = lb[1] + cs;
        max_bound[2] = lb[2] + cs;
        sub_pos = where(&k, points, min_bound, max_bound);
        if (k > 0)
        {
            current_box-> subBoxes[j] = *create_box(k, sub_pos, min_bound, current_box-> box_half_size);
            j ++;
        }
        if((j!=8) && (j>0))
        {
            current_box->subBoxes = realloc(current_box->subBoxes, j*sizeof(box));
        }
    }
    current_box-> number_of_subs = j;
    free(min_bound);
    free(max_bound);
    free(cm);
    free(points);
    return current_box;
}
