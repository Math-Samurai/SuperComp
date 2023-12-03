#include <stdint.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

enum POSITION
{
    OUTSIDE,
    INSIDE,
    INTERSECT
};

struct Point
{
    double x;
    double y;
} typedef point;

void free_matrix(double **A, int32_t row_count)
{
    if (!A)
    {
        return;
    }

    for (int32_t i = 0; i < row_count; ++i)
    {
        free(A[i]);
    }

    free(A);
}

enum POSITION calc_point_position(const point *p)
{
    if (p -> x <= 2 || p -> y <= -3 * (p -> x) + 9)
    {
        return INSIDE;
    }
    
    return OUTSIDE;
}

enum POSITION calc_position(const point *point1, const point *point2)
{
    enum POSITION pos1 = calc_point_position(point1);
    enum POSITION pos2 = calc_point_position(point2);
    if (pos1 != pos2)
    {
        return INTERSECT;
    }
    else if (pos1 == INSIDE)
    {
        return INSIDE;
    }
    else
    {
        return OUTSIDE;
    }
}

double *multiply_matrix_by_vector(double **mx, double *vector, int32_t size)
{
    double *res = (double *)malloc(sizeof(double) * size);
    
#pragma omp parallel
{
    int id = omp_get_thread_num();
    int numt = omp_get_num_threads();
    for (int32_t i = id; i < size; i += numt)
    {
        double tmp = 0;
        for (int32_t j = 0; j < size; ++j)
        {
            tmp += mx[i][j] * vector[j];
        }

        res[i] = tmp;
    }
}
    return res;
}

double *multiply_vector_by_number(double *vector, int32_t size, double num)
{
    double *res = (double *)malloc(sizeof(double) * size);
#pragma omp parallel
{
    int id = omp_get_thread_num();
    int numt = omp_get_num_threads();
    for (int32_t i = id; i < size; i += numt)
    {
        res[i] = vector[i] * num;
    }
}
    return res;
}

double *add_vectors(double *vec1, double *vec2, int32_t size)
{
    double *res = (double *)malloc(sizeof(double) * size);
#pragma omp parallel
{
    int id = omp_get_thread_num();
    int numt = omp_get_num_threads();
    for (int32_t i = id; i < size; i += numt)
    {
        res[i] = vec1[i] + vec2[i];
    }
}
    return res;
}

double *subtract_vectors(double *vec1, double *vec2, int32_t size)
{
    double *res = (double *)malloc(sizeof(double) * size);
#pragma omp parallel
{
    int id = omp_get_thread_num();
    int numt = omp_get_num_threads();
    for (int32_t i = id; i < size; i += numt)
    {
        res[i] = vec1[i] - vec2[i];
    }
}
    return res;
}

double scalar_product(double *vec1, double *vec2, int32_t size)
{
    double res = 0;
    for (size_t i = 0; i < size; ++i)
    {
        res += vec1[i] * vec2[i];
    }

    return res;
}

long double calc_snorm(double *vec, int32_t size)
{
    long double res = 0;
    for (int32_t i = 0; i < size; ++i)
    {
        res += vec[i] * vec[i];
    }

    return res;
}

double calc_a_coef(int32_t first_ind, int32_t second_ind, double step)
{
    point p1 = {(first_ind - 0.5) * step, (second_ind - 0.5) * step};
    point p2 = {(first_ind - 0.5) * step, (second_ind + 0.5) * step};
    if (p1.x <= 2 || p2.y <= -3 * (p2.x) + 9)
    {
        return 1;
    }
    else if (p1.y > -3 * (p1.x) + 9)
    {
        return 1 / (step * step);
    }
    else
    {
        double coef = -3 * (p1.x) + 9 - p1.y;
        return coef / step + (1 - coef / step) / (step * step);
    }
}

double calc_b_coef(int32_t first_ind, int32_t second_ind, double step)
{
    point p1 = {(first_ind - 0.5) * step, (second_ind - 0.5) * step};
    point p2 = {(first_ind + 0.5) * step, (second_ind - 0.5) * step};
    if (p2.x <= 2 || p2.y <= -3 * (p2.x) + 9)
    {
        return 1;
    }
    else if (p1.x > (p1.y - 9) / (-3))
    {
        return 1 / (step * step);
    }
    else
    {
        double coef = (p1.y - 9) / (-3) - p1.x;
        return coef / step + (1 - coef / step) / (step * step);
    }
}

int get_point_position(const point *p)
{
    if (p -> x <= 2 || p -> y <= -3 * (p -> x) + 9)
    {
        return INSIDE;
    }

    return OUTSIDE;
}

double calc_f(int32_t first_ind, int32_t second_ind, double step)
{
    point p1 = {(first_ind - 0.5) * step, (second_ind - 0.5) * step};
    point p2 = {(first_ind + 0.5) * step, (second_ind - 0.5) * step};
    point p3 = {(first_ind + 0.5) * step, (second_ind + 0.5) * step};
    point p4 = {(first_ind - 0.5) * step, (second_ind + 0.5) * step};
    int p1_pos = get_point_position(&p1);
    int p2_pos = get_point_position(&p2);
    int p3_pos = get_point_position(&p3);
    int p4_pos = get_point_position(&p4);
    if (p3_pos == INSIDE)
    {
        return 1;
    }
    else if (p1_pos == OUTSIDE)
    {
        return 0;
    }
    else if (p1_pos == INSIDE && p2_pos == INSIDE && p3_pos == OUTSIDE && p4_pos == INSIDE)
    {
        point intersection_point1 = {(p4.y - 9) / (-3), p4.y};
        point intersection_point2 = {p2.x, -3 * (p2.x) + 9};
        return (step * step - (p3.x - intersection_point1.x) * ((intersection_point2.y - p2.y) / 2)) / (step * step);
    }
    else if (p1_pos == INSIDE && p2_pos == OUTSIDE && p3_pos == OUTSIDE && p4_pos == INSIDE)
    {
        point intersection_point1 = {(p4.y - 9) / (-3), p4.y};
        point intersection_point2 = {(p1.y - 9) / 3, p1.y};
        return (((intersection_point1.x - p4.x) + (intersection_point2.x - p1.x)) / 2) / step;
    }
    else if (p1_pos == INSIDE && p2_pos == OUTSIDE && p3_pos == OUTSIDE && p4_pos == OUTSIDE)
    {
        point intersection_point1 = {p1.x, -3 * p1.x + 9};
        point intersection_point2 = {-3 * p1.x + 9, p1.y};
        return ((intersection_point1.y - p1.y) * (intersection_point2.x - p1.x) / 2) / (step * step);
    }
    else
    {
        point intersection_point1 = {p1.x, -3 * p1.x + 9};
        point intersection_point2 = {p2.x, -3 * p2.x + 9};
        return (((intersection_point1.y - p1.y) + (intersection_point2.y - p2.y)) / 2) / step;
    }
}

double **create_matrix_A(double *right, int32_t node_count_on_axis, double step)
{
    double **A = (double **)malloc(sizeof(double *) * (node_count_on_axis - 2) * (node_count_on_axis - 2));

    for (int32_t node_index = 0, row_index = 0; node_index < node_count_on_axis * node_count_on_axis; ++node_index)
    {
        int32_t row = node_index / node_count_on_axis;
        int32_t col = node_index % node_count_on_axis;
        if (node_index < node_count_on_axis || node_index >= node_count_on_axis * (node_count_on_axis - 1) || col == 0 || col == node_count_on_axis - 1)
        {
            continue;
        }

        A[row_index] = (double *)malloc(sizeof(double) * (node_count_on_axis - 2) * (node_count_on_axis - 2));
        double a1 = calc_a_coef(row + 1, col, step) / (step * step);
        double a2 = calc_a_coef(row, col, step) / (step * step);
        double b1 = calc_b_coef(row, col + 1, step) / (step * step);
        double b2 = calc_b_coef(row, col, step) / (step * step);
        right[row_index] = calc_f(row, col, step);;
        if (row_index + node_count_on_axis < (node_count_on_axis - 2) * (node_count_on_axis - 2))
        {
            A[row_index][row_index + node_count_on_axis] = a1;
        }
        if (row_index - node_count_on_axis >= 0)
        {
            A[row_index][row_index - node_count_on_axis] = a2;
        }
        if (row_index + 1 < (node_count_on_axis - 2) * (node_count_on_axis - 2))
        {
            A[row_index][row_index + 1] = b1;
        }
        A[row_index][row_index - 1] = b2;
        A[row_index][row_index] = a1 + a2 + b1 + b2;
        ++row_index;
    }

    return A;
}

double *solve_system(double **mx, double *right, int32_t size, double delta)
{
    double *diff = (double *)malloc(sizeof(double) * size);
    double *sol = NULL;
    double *prev_sol = (double *)malloc(sizeof(double) * size);
    if (!prev_sol || !diff)
    {
        return NULL;
    }

    for (int32_t i = 0; i < size; ++i)
    {
        prev_sol[i] = 0;
        diff[i] = -right[i];
    }

    double *result = multiply_matrix_by_vector(mx, diff, size);
    while (1)
    {
        double numinator = scalar_product(result, diff, size);
        double denominator = calc_snorm(result, size);

        double t = numinator / denominator;
        double *subtrahend = multiply_vector_by_number(diff, size, t);

        sol = subtract_vectors(prev_sol, subtrahend, size);

        double *sol_diff = subtract_vectors(sol, prev_sol, size);

        long double diff_norm = calc_snorm(sol_diff, size);
        if (sqrt(diff_norm) < delta)
        {
            return sol;
        }

        prev_sol = sol;
        result = multiply_matrix_by_vector(mx, sol, size);
        diff = subtract_vectors(result, right, size);
        result = multiply_matrix_by_vector(mx, diff, size);
    }
}

int main(int argc, char *argv[])
{
    const point rec_point_1 = {0, 0};
    const point rec_point_2 = {0, 3};
    const point rec_point_3 = {3, 0};
    const point rec_point_4 = {3, 3};
    const double rec_width = rec_point_4.x - rec_point_1.x;
    const double rec_height = rec_point_4.y - rec_point_1.y;
    int32_t rec_count = strtol(argv[1], NULL, 10);
    --rec_count;
    int32_t mx_size = (rec_count - 1) * (rec_count - 1);
    double *right = (double *)malloc(sizeof(double) * mx_size);
    double **A = create_matrix_A(right, rec_count + 1, rec_width / rec_count);
    double *res = solve_system(A, right, mx_size, 1E-6);
    FILE *fl = fopen("par_result", "w");
    double zero = 0;
    for (int i = 0; i < rec_count + 1; ++i)
    {
        fprintf(fl, "%0.8lf ", zero);
    }

    fprintf(fl, "\n");
    for (size_t i = 0; i < rec_count - 1; ++i)
    {
        for (size_t j = 0; j < rec_count - 1; ++j)
        {
            if (j == 0)
            {
                fprintf(fl, "%0.8lf %0.8lf ", zero, res[i * (rec_count - 1) + j]);
            }
            else if (j == rec_count - 2)
            {
                fprintf(fl, "%0.8lf %0.8lf\n", res[i * (rec_count - 1) + j], zero);
            }
            else
            {
                fprintf(fl, "%0.8lf ", res[i * (rec_count - 1) + j]);
            }
        }
        
    }

    for (int i = 0; i < rec_count + 1; ++i)
    {
        fprintf(fl, "%0.8lf ", zero);
    }
    fclose(fl);
    fclose(mx);
    fclose(b);
    return 0;
}
