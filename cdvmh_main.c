#include <stdint.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

#define NODE_COUNT_ON_AXIS 160
#define NODE_COUNT 24964

#define OUTSIDE 0
#define INSIDE 1
#define INTERSECT 2

#pragma dvm array distribute[block][]
double mx[NODE_COUNT][NODE_COUNT];
#pragma dvm array align([i] with mx[][i])
double right[NODE_COUNT];
#pragma dvm array align([i] with mx[i][])
double result[NODE_COUNT];
#pragma dvm array align([i] with mx[][i])
double diff[NODE_COUNT];
#pragma dvm array align([i] with result[i])
double prev_sol[NODE_COUNT];
#pragma dvm array align([i] with result[i])
double sol[NODE_COUNT];
#pragma dvm array align([i] with result[i])
double subtrahend[NODE_COUNT];
#pragma dvm array align([i] with result[i])
double sol_diff[NODE_COUNT];

int main()
{
    int max = 100;
    int it = 0;
    double delta = 1E-6;
    double step = 3.0 / NODE_COUNT_ON_AXIS;
#pragma dvm region
{
#pragma dvm parallel([i][j] on mx[i][j])
    for (int i = 0; i < NODE_COUNT; i++)
    {
        for (int j = 0; j < NODE_COUNT; j++)
        {
            double coef = 0;
            double x1, x2, y1, y2;
            if (j == i + NODE_COUNT_ON_AXIS || i == j)
            {
                double first_ind = i / NODE_COUNT_ON_AXIS + 1;
                double second_ind = i % NODE_COUNT_ON_AXIS;
                double x1 = (first_ind - 0.5) * step;
                double y1 = (second_ind - 0.5) * step;
                double x2 = (first_ind - 0.5) * step;
                double y2 = (second_ind + 0.5) * step;
                if (x1 <= 2 || y2 <= -3 * x2 + 9)
                {
                    coef += 1;
                }
                else if (y1 > -3 * x1 + 9)
                {
                    coef += 1 / (step * step);
                }
                else
                {
                    double tmp = -3 * (x1) + 9 - y1;
                    coef += tmp / step + (1 - tmp / step) / (step * step);
                }
	        }
            else if (j == i - NODE_COUNT_ON_AXIS || i == j)
            {
                double first_ind = i / NODE_COUNT_ON_AXIS;
                double second_ind = i % NODE_COUNT_ON_AXIS;
                double x1 = (first_ind - 0.5) * step;
                double y1 = (second_ind - 0.5) * step;
                double x2 = (first_ind - 0.5) * step;
                double y2 = (second_ind + 0.5) * step;
                if (x1 <= 2 || y2 <= -3 * x2 + 9)
                {
                    coef += 1;
                }
                else if (y1 > -3 * x1 + 9)
                {
                    coef += 1 / (step * step);
                }
                else
                {
                    double tmp = -3 * (x1) + 9 - y1;
                    coef += tmp / step + (1 - tmp / step) / (step * step);
                }
            }
            else if (j == i - 1 || i == j)
            {
                double first_ind = i / NODE_COUNT_ON_AXIS;
                double second_ind = i % NODE_COUNT_ON_AXIS;
                double x1 = (first_ind - 0.5) * step;
                double y1 = (second_ind - 0.5) * step;
                double x2 = (first_ind + 0.5) * step;
                double y2 = (second_ind - 0.5) * step;
                if (x2 <= 2 || y2 <= -3 * x2 + 9)
                {
                    coef += 1;
                }
                else if (x1 > (y1 - 9) / (-3))
                {
                    coef += 1 / (step * step);
                }
                else
                {
                    double tmp = (y1 - 9) / (-3) - x1;
                    coef += tmp / step + (1 - tmp / step) / (step * step);
                }
	        }
            else if (j == i + 1 || i == j)
            {
                double first_ind = i / NODE_COUNT_ON_AXIS;
                double second_ind = i % NODE_COUNT_ON_AXIS + 1;
                double x1 = (first_ind - 0.5) * step;
                double y1 = (second_ind - 0.5) * step;
                double x2 = (first_ind + 0.5) * step;
                double y2 = (second_ind - 0.5) * step;
                if (x2 <= 2 || y2 <= -3 * x2 + 9)
                {
                    coef += 1;
                }
                else if (x1 > (y1 - 9) / (-3))
                {
                    coef += 1 / (step * step);
                }
                else
                {
                    double tmp = (y1 - 9) / (-3) - x1;
                    coef += tmp / step + (1 - tmp / step) / (step * step);
                }
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
	double ix1, iy1, ix2, iy2;
        double first_ind = i / NODE_COUNT_ON_AXIS;
        double second_ind = i % NODE_COUNT_ON_AXIS;
        double x1 = (first_ind - 0.5) * step;
        double y1 = (second_ind - 0.5) * step;
        double x2 = (first_ind + 0.5) * step;
        double y2 = (second_ind - 0.5) * step;
        double x3 = (first_ind + 0.5) * step;
        double y3 = (second_ind + 0.5) * step;
        double x4 = (first_ind - 0.5) * step;
        double y4 = (second_ind + 0.5) * step;
        int p1_pos, p2_pos, p3_pos, p4_pos;
	if (x1 <= 2 || y1 <= -3 * x1 + 9)
	{
            p1_pos = INSIDE;
	}
	else
	{
            p1_pos = OUTSIDE;
	}

	if (x2 <= 2 || y2 <= -3 * x2 + 9)
	{
	    p2_pos = INSIDE;	
	}
        else
        {
	    p2_pos = OUTSIDE;
	}

	if (x3 <= 2 || y3 <= -3 * x3 + 9)
	{
	    p3_pos = INSIDE;
	}
	else
	{
	    p3_pos = OUTSIDE;
	}

	if (x4 <= 2 || y4 <= -3 * x4 + 9)
	{
	    p4_pos = INSIDE;
	}
	else
	{
	    p4_pos = OUTSIDE;
	}

        if (p3_pos == INSIDE)
        {
            right[i] = 1;
        }
        else if (p1_pos == OUTSIDE)
        {
            right[i] = 0;
        }
        else if (p1_pos == INSIDE && p2_pos == INSIDE && p3_pos == OUTSIDE && p4_pos == INSIDE)
        {
            ix1 = (y4 - 9) / (-3);
            iy1 = y4;
            ix2 = x2;
            iy2 = -3 * (x2) + 9;
            right[i] = (step * step - (x3 - ix1) * ((iy2 - y2) / 2)) / (step * step);
        }
        else if (p1_pos == INSIDE && p2_pos == OUTSIDE && p3_pos == OUTSIDE && p4_pos == INSIDE)
        {
            ix1 = (y4 - 9) / (-3);
            iy1 = y4;
            ix2 = (y1 - 9) / 3;
            iy2 = y1;
            right[i] = (((ix1 - x4) + (ix2 - x1)) / 2) / step;
        }
        else if (p1_pos == INSIDE && p2_pos == OUTSIDE && p3_pos == OUTSIDE && p4_pos == OUTSIDE)
        {
            ix1 = x1;
            iy1 = -3 * x1 + 9;
            ix2 = -3 * x1 + 9;
            iy2 = y1;
            right[i] = ((iy1 - y1) * (ix2 - x1) / 2) / (step * step);
        }
        else
        {
            ix1 = x1;
            iy1 = -3 * x1 + 9;
            ix2 = x2;
            iy2 = -3 * x2 + 9;
            right[i] = (((iy1 - y1) + (iy2 - y2)) / 2) / step;
        }
    }
}

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
        }
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
        for (int i = 0; i < NODE_COUNT; ++i)
        {
            sol[i] = prev_sol[i] - subtrahend[i];
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

        double n = 0;
#pragma dvm region
{
#pragma dvm parallel ([i] on mx[i][]) reduction(sum(n))
        for (int i = 0; i < NODE_COUNT; i++)
        {
            n += sol_diff[i] * sol_diff[i];
        }
}
        if (it > max)
        {
            break;
        }
	++it;

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
