#include <stdint.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#define NODE_COUNT_ON_AXIS 160
#define NODE_COUNT 24964

#define OUTSIDE 0
#define INSIDE 1
#define INTERSECT 2

int calc_point_position(double x, double y)
{
    if (x <= 2 || y <= -3 * x + 9)
    {
        return INSIDE;
    }
    
    return OUTSIDE;
}

int calc_position(double x1, double y1, double x2, double y2)
{
    int pos1 = calc_point_position(x1, y1);
    int pos2 = calc_point_position(x2, y2);
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
    
    int i, j;
    for (i = 0; i < size; ++i)
    {
        double tmp = 0;
        for (j = 0; j < size; ++j)
        {
            tmp += mx[i][j] * vector[j];
        }

        res[i] = tmp;
    }
    return res;
}

double calc_a_coef(int32_t first_ind, int32_t second_ind, double step)
{
    double x1 = (first_ind - 0.5) * step;
    double y1 = (second_ind - 0.5) * step;
    double x2 = (first_ind - 0.5) * step;
    double y2 = (second_ind + 0.5) * step;
    if (x1 <= 2 || y2 <= -3 * x2 + 9)
    {
        return 1;
    }
    else if (y1 > -3 * x1 + 9)
    {
        return 1 / (step * step);
    }
    else
    {
        double coef = -3 * (x1) + 9 - y1;
        return coef / step + (1 - coef / step) / (step * step);
    }
}

double calc_b_coef(int32_t first_ind, int32_t second_ind, double step)
{
    double x1 = (first_ind - 0.5) * step;
    double y1 = (second_ind - 0.5) * step;
    double x2 = (first_ind + 0.5) * step;
    double y2 = (second_ind - 0.5) * step;
    if (x2 <= 2 || y2 <= -3 * x2 + 9)
    {
        return 1;
    }
    else if (x1 > (y1 - 9) / (-3))
    {
        return 1 / (step * step);
    }
    else
    {
        double coef = (y1 - 9) / (-3) - x1;
        return coef / step + (1 - coef / step) / (step * step);
    }
}

int get_point_position(double x, double y)
{
    if (x <= 2 || y <= -3 * x + 9)
    {
        return INSIDE;
    }

    return OUTSIDE;
}

double calc_f(int32_t first_ind, int32_t second_ind, double step)
{
    double ix1, iy1, ix2, iy2;
    double x1 = (first_ind - 0.5) * step;
    double y1 = (second_ind - 0.5) * step;
    double x2 = (first_ind + 0.5) * step;
    double y2 = (second_ind - 0.5) * step;
    double x3 = (first_ind + 0.5) * step;
    double y3 = (second_ind + 0.5) * step;
    double x4 = (first_ind - 0.5) * step;
    double y4 = (second_ind + 0.5) * step;
    int p1_pos = get_point_position(x1, y1);
    int p2_pos = get_point_position(x2, y2);
    int p3_pos = get_point_position(x3, y3);
    int p4_pos = get_point_position(x4, y4);
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
        ix1 = (y4 - 9) / (-3);
        iy1 = y4;
        ix2 = x2;
        iy2 = -3 * (x2) + 9;
        return (step * step - (x3 - ix1) * ((iy2 - y2) / 2)) / (step * step);
    }
    else if (p1_pos == INSIDE && p2_pos == OUTSIDE && p3_pos == OUTSIDE && p4_pos == INSIDE)
    {
        ix1 = (y4 - 9) / (-3);
        iy1 = y4;
        ix2 = (y1 - 9) / 3;
        iy2 = y1;
        return (((ix1 - x4) + (ix2 - x1)) / 2) / step;
    }
    else if (p1_pos == INSIDE && p2_pos == OUTSIDE && p3_pos == OUTSIDE && p4_pos == OUTSIDE)
    {
        ix1 = x1;
        iy1 = -3 * x1 + 9;
        ix2 = -3 * x1 + 9;
        iy2 = y1;
        return ((iy1 - y1) * (ix2 - x1) / 2) / (step * step);
    }
    else
    {
        ix1 = x1;
        iy1 = -3 * x1 + 9;
        ix2 = x2;
        iy2 = -3 * x2 + 9;
        return (((iy1 - y1) + (iy2 - y2)) / 2) / step;
    }
}

int main()
{
    double delta = 1E-6;
    double step = 3.0 / NODE_COUNT_ON_AXIS;
#pragma dvm array distribute[block][]
    double mx[NODE_COUNT][NODE_COUNT];
#pragma dvm array align([i] with mx[][i])
    double right[NODE_COUNT];
#pragma dvm region
{
#pragma dvm parallel([i][j] on mx[i][j])
    for (int i = 0; i < NODE_COUNT; i++)
    {
        for (int j = 0; j < NODE_COUNT; j++)
        {
            double coef = 0;
            if (j == i + NODE_COUNT_ON_AXIS)
            {
                coef = calc_a_coef(i / NODE_COUNT_ON_AXIS + 1, i % NODE_COUNT_ON_AXIS, step);
	    }
            else if (j == i - NODE_COUNT_ON_AXIS)
            {
                coef = calc_a_coef(i / NODE_COUNT_ON_AXIS, i % NODE_COUNT_ON_AXIS, step);
            }
            else if (j == i - 1)
            {
                coef = calc_b_coef(i / NODE_COUNT_ON_AXIS, i % NODE_COUNT_ON_AXIS, step);
	    }
            else if (j == i + 1)
            {
                coef = calc_b_coef(i / NODE_COUNT_ON_AXIS, i % NODE_COUNT_ON_AXIS + 1, step);
            }
            else if (i == j)
            {
                coef = calc_a_coef(i / NODE_COUNT_ON_AXIS + 1, i % NODE_COUNT_ON_AXIS, step) + calc_a_coef(i / NODE_COUNT_ON_AXIS, i % NODE_COUNT_ON_AXIS, step) +
                    calc_b_coef(i / NODE_COUNT_ON_AXIS, i % NODE_COUNT_ON_AXIS, step) + calc_b_coef(i / NODE_COUNT_ON_AXIS, i % NODE_COUNT_ON_AXIS + 1, step);
            }

            mx[i][j] = coef / (step * step);
        }
    }
}

#pragma dvm region
{
#pragma dvm parallel([i] on mx[][i])
    for (int i = 0; i < NODE_COUNT; i++)
    {
        right[i] = calc_f(i / NODE_COUNT_ON_AXIS, i % NODE_COUNT_ON_AXIS, step);
    }
}

#pragma dvm array align([i] with mx[i][])
    double result[NODE_COUNT];
#pragma dvm array align([i] with mx[][i])
    double diff[NODE_COUNT];
#pragma dvm array align([i] with result[i])
    double prev_sol[NODE_COUNT];
#pragma dvm array align([i] with result[i])
    double sol[NODE_COUNT];
#pragma dvm array align([i] with result[i])
    double sol_diff[NODE_COUNT];
#pragma dvm array align([i] with result[i])
    double subtrahend[NODE_COUNT];

#pragma dvm region
{
#pragma dvm parallel([i] on mx[i][])
    for (int i = 0; i < NODE_COUNT; i++)
    {
        prev_sol[i] = 0;
    }
}

#pragma dvm region
{
#pragma dvm parallel([i] on mx[][i])
    for (int i = 0; i < NODE_COUNT; i++)
    {
        diff[i] = -right[i];
    }
}
    while (1)
    {
#pragma dvm region
{
#pragma dvm parallel([i] on mx[i][])
        for (int i = 0; i < NODE_COUNT; i++)
        {
            double num = 0;
            for (int j = 0; j < NODE_COUNT; j++)
            {
                num += mx[i][j] * diff[j];
            }

            result[i] = num;
        };
}
    
        double scalar_product = 0;
#pragma dvm region
{
#pragma dvm parallel([i] on mx[i][]) reduction(sum(scalar_product))
        for (int i = 0; i < NODE_COUNT; i++)
        {
            scalar_product += diff[i] * result[i];
        }
}

        double snorm = 0;
#pragma dvm region
{
#pragma dvm parallel([i] on mx[i][]) reduction(sum(snorm))
        for (int i = 0; i < NODE_COUNT; i++)
        {
            snorm += result[i] * result[i];
        }
}

        double t = scalar_product / snorm;
#pragma dvm region
{
#pragma dvm parallel([i] on mx[i][])
        for (int i = 0; i < NODE_COUNT; i++)
        {
            subtrahend[i] = diff[i] * t;
        }
}

#pragma dvm region
{
#pragma dvm parallel([i] on mx[i][])
        for (int i = 0; i < NODE_COUNT; i++)
        {
            sol_diff[i] = sol[i] - prev_sol[i];
        }
}

        snorm = 0;
#pragma dvm region
{
#pragma dvm parallel ([i] on mx[i][]) reduction (sum(snorm))
        for (int i = 0; i < NODE_COUNT; i++)
        {
            snorm += sol_diff[i] * sol_diff[i];
        }
}

        if (snorm * snorm < delta * delta)
        {
            break;
        }

#pragma dvm region
{
#pragma dvm parallel([i] on mx[i][])
        for (int i = 0; i < NODE_COUNT; i++)
        {
            prev_sol[i] = sol[i];
        }
}

#pragma dvm region
{
#pragma dvm parallel([i] on mx[i][])
        for (int i = 0; i < NODE_COUNT; i++)
        {
            double num = 0;
            for (int j = 0; j < NODE_COUNT; j++)
            {
                num += mx[i][j] * sol[j];
            }

            result[i] = num;
        }
}

#pragma dvm region
{
#pragma dvm parallel([i] on mx[i][])
        for (int i = 0; i < NODE_COUNT; i++)
        {
            diff[i] = result[i] - right[i];
        }
}
    }
    return 0;
}
