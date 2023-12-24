#include <stdint.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>

int proc_count, rank;

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
    return res;
}

double *multiply_vector_by_number(double *vector, int32_t size, double num)
{
    double *res = (double *)malloc(sizeof(double) * size);
    int id = omp_get_thread_num();
    int numt = omp_get_num_threads();
    int32_t i;
    for (i = id; i < size; i += numt)
    {
        res[i] = vector[i] * num;
    }
    return res;
}

double *add_vectors(double *vec1, double *vec2, int32_t size)
{
    double *res = (double *)malloc(sizeof(double) * size);
    int id = omp_get_thread_num();
    int numt = omp_get_num_threads();
    int32_t i;
    for (i = id; i < size; i += numt)
    {
        res[i] = vec1[i] + vec2[i];
    }
    return res;
}

double *subtract_vectors(double *vec1, double *vec2, int32_t size)
{
    double *res = (double *)malloc(sizeof(double) * size);
    int id = omp_get_thread_num();
    int numt = omp_get_num_threads();
    int32_t i;
    for (i = id; i < size; i += numt)
    {
        res[i] = vec1[i] - vec2[i];
    }
    return res;
}

double scalar_product(double *vec1, double *vec2, int32_t size)
{
    double res = 0;
    int32_t i;
    for (i = 0; i < size; ++i)
    {
        res += vec1[i] * vec2[i];
    }

    return res;
}

long double calc_snorm(double *vec, int32_t size)
{
    long double res = 0;
    int32_t i;
    for (i = 0; i < size; ++i)
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

    int32_t node_index, row_index;
    for (node_index = 0, row_index = 0; node_index < node_count_on_axis * node_count_on_axis && row_index < (node_count_on_axis - 2) * (node_count_on_axis - 2); ++node_index)
    {
        int32_t row = node_index / node_count_on_axis;
        int32_t col = node_index % node_count_on_axis;
        if (node_index < node_count_on_axis || node_index >= node_count_on_axis * (node_count_on_axis - 1) || col == 0 || col == node_count_on_axis - 1)
        {
            continue;
        }

        A[row_index] = (double *)malloc(sizeof(double) * (node_count_on_axis - 2) * (node_count_on_axis - 2));
        for (int j = 0; j < (node_count_on_axis - 2) * (node_count_on_axis - 2); ++j)
        {
            A[row_index][j] = 0;
        }

        double a1 = calc_a_coef(row + 1, col, step) / (step * step);
        double a2 = calc_a_coef(row, col, step) / (step * step);
        double b1 = calc_b_coef(row, col + 1, step) / (step * step);
        double b2 = calc_b_coef(row, col, step) / (step * step);
        right[row_index] = calc_f(row, col, step);
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
        if (row_index - 1 >= 0)
        {
            A[row_index][row_index - 1] = b2;
        }
        A[row_index][row_index] = a1 + a2 + b1 + b2;
        ++row_index;
    }

    return A;
}

int get_vec_size(int size, int rank)
{
    int vec_size;
    if (size % proc_count == 0 || rank != proc_count - 1)
    {
        return size / proc_count;
    }
    
    return size / proc_count + size % proc_count;
}

double *solve_system(double **mx, double *right, int32_t size, double delta)
{
    int32_t vec_size = get_vec_size(size, rank);
    int i, j;
    if (!rank)
    {
        double *diff = (double *)malloc(sizeof(double) * size);
        double *result = (double *)malloc(sizeof(double) * size);
        double *sol = (double *)malloc(sizeof(double) * size);
        double *prev_sol = (double *)malloc(sizeof(double) * size);
        double *subtrahend = (double *)malloc(sizeof(double) * size);

        for (i = 0; i < size; ++i)
        {
            prev_sol[i] = 0;
            diff[i] = -right[i];
        }


        while (1)
        {
            for (i = 1; i < proc_count; ++i)
            {
                MPI_Send(diff, size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }

            for (i = 0; i < vec_size; ++i)
            {
                result[i] = 0;
                for (j = 0; j < size; ++j)
                {
                    result[i] += mx[i][j] * diff[j];
                }
            }

            MPI_Status stat;
            for (i = 1; i < proc_count; ++i)
            {
                MPI_Recv(result + i * (size / proc_count), get_vec_size(size, i), MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &stat);
            }
            
            double numinator = 0;
            double denominator = 0;
            double *tmp = (double *)malloc(sizeof(double) * 2);
            for (i = 0; i < vec_size; ++i)
            {
                numinator += result[i] * diff[i];
                denominator += result[i] * result[i];
            }
            
            for (i = 1; i < proc_count; ++i)
            {
                MPI_Recv(tmp, 2, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &stat);
                numinator += tmp[0];
                denominator += tmp[1];
            }
            
            double t = numinator / denominator;
            for (i = 1; i < proc_count; ++i)
            {
                MPI_Send(&t, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
            }

            for (i = 0; i < vec_size; ++i)
            {
                subtrahend[i] = t * diff[i];
            }

            for(i = 0; i < vec_size; ++i)
            {
                sol[i] = prev_sol[i] - subtrahend[i];
            }

            for (i = 1; i < proc_count; ++i)
            {
                MPI_Send(prev_sol + i * (size / proc_count), get_vec_size(size, i), MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
            }

            for (i = 1; i < proc_count; ++i)
            {
                MPI_Recv(sol + i * (size / proc_count), get_vec_size(size, i), MPI_DOUBLE, i, 5, MPI_COMM_WORLD, &stat);
            }

            printf("W");
            double norm = 0;
            for (i = 0; i < vec_size; ++i)
            {
                norm += (sol[i] - prev_sol[i]) * (sol[i] - prev_sol[i]);
            }

            for (i = 1; i < proc_count; ++i)
            {
                double tmp;
                MPI_Recv(&tmp, 1, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, &stat);
                norm += tmp;
            }

            int code = 0;
            if (sqrt(norm) < delta)
            {
                code = 1;
                for (i = 1; i < proc_count; ++i)
                {
                    MPI_Send(&code, 1, MPI_INT, i, 7, MPI_COMM_WORLD);
                }

                return sol;
            }

            for (i = 1; i < proc_count; ++i)
            {
                MPI_Send(&code, 1, MPI_INT, i, 7, MPI_COMM_WORLD);
            }

            prev_sol = sol;
            for (i = 1; i < proc_count; ++i)
            {
                MPI_Send(sol, size, MPI_DOUBLE, i, 8, MPI_COMM_WORLD);
            }

            for (i = 0; i < vec_size; ++i)
            {
                result[i] = 0;
                for (j = 0; j < size; ++j)
                {
                    result[i] += mx[i][j] * sol[j];
                }

                diff[i] = result[i] - right[i];
            }

            for (i = 1; i < proc_count; ++i)
            {
                MPI_Recv(diff + i * (size / proc_count), get_vec_size(size, i), MPI_DOUBLE, i, 9, MPI_COMM_WORLD, &stat);
            }
        }
    }
    else
    {
        double *diff_buf = (double *)malloc(sizeof(double) * size);
        double *result_buf = (double *)malloc(sizeof(double) * vec_size);
        double *subtrahend = (double *)malloc(sizeof(double) * vec_size);
        double *prev_sol = (double *)malloc(sizeof(double) * vec_size);
        double *sol = (double *)malloc(sizeof(double) * vec_size);
        double *right = (double *)malloc(sizeof(double) * vec_size);
        double *full_sol = (double *)malloc(sizeof(double) * size);
        int start_point = rank * (size / proc_count);
        

        MPI_Status stat;
        while (1)
        {
            MPI_Recv(diff_buf, size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &stat);

            for (i = start_point; i < start_point + vec_size; ++i)
            {
                result_buf[i - start_point] = 0;
                for (j = 0; j < size; ++j)
                {
                    result_buf[i - start_point] += mx[i][j] * diff_buf[j];
                }
            }

            MPI_Send(result_buf, vec_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            double *nums = (double *)malloc(sizeof(double) * 2);
            nums[0] = 0;
            nums[1] = 0;
            for (i = start_point; i < start_point + vec_size; ++i)
            {
                nums[0] += result_buf[i - start_point] * diff_buf[i];
                nums[1] += result_buf[i - start_point] * result_buf[i - start_point];
            }

            MPI_Send(nums, 2, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
            double t;
            MPI_Recv(&t, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &stat);
            for (i = start_point; i < start_point + vec_size; ++i)
            {
                subtrahend[i - start_point] = t * diff_buf[i];
            }

            MPI_Recv(prev_sol, vec_size, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &stat);
            for (i = 0; i < vec_size; ++i)
            {
                sol[i] = prev_sol[i] - subtrahend[i];
            }
            MPI_Send(sol, vec_size, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
            double norm = 0;
            for (i = 0; i < vec_size; ++i)
            {
                norm += (sol[i] - prev_sol[i]) * (sol[i] - prev_sol[i]);
            }

            MPI_Send(&norm, 1, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD);
            int code;
            MPI_Recv(&code, 1, MPI_INT, 0, 7, MPI_COMM_WORLD, &stat);
            if (code == 1)
            {
                return NULL;
            }

            MPI_Recv(full_sol, size, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD, &stat);
            for (i = start_point; i < start_point + vec_size; ++i)
            {
                result_buf[i - start_point] = 0;
                for (j = 0; j < size; ++j)
                {
                    result_buf[i - start_point] += mx[i][j] * full_sol[j];
                }
                diff_buf[i] = result_buf[i - start_point] - right[i];
            }

            MPI_Send(diff_buf + start_point, vec_size, MPI_DOUBLE, 0, 9, MPI_COMM_WORLD);
        }
    }

    /*double *result = multiply_matrix_by_vector(mx, diff, size);
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
    }*/
}

int main(int argc, char *argv[])
{
    const double rec_width = 3;
    const double rec_height = 3;
    int32_t rec_count = strtol(argv[1], NULL, 10);
    --rec_count;
    int32_t mx_size = (rec_count - 1) * (rec_count - 1);
    double *right = (double *)malloc(sizeof(double) * mx_size);
    double **A = create_matrix_A(right, rec_count + 1, rec_width / rec_count);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_count);
    double *res = solve_system(A, right, mx_size, 1E-6);
    if (!rank)
    {
        FILE *fl = fopen("mpi_result", "w");
        double zero = 0;
        for (int i = 0; i < rec_count + 1; ++i)
        {
            fprintf(fl, "%lf ", zero);
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
    }
    MPI_Finalize();
    return 0;
}
