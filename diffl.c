// diffl.cpp : Defines the exported functions for the DLL application.
//

//#include "stdafx.h"

// diffusion.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"

//#include "stddef.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "diffl.h"
#include <string.h>

//#define simpleExport extern __declspec(dllexport)

struct _doubleArray {
    double* m_elements;
    int     m_size;
};

typedef struct _doubleArray doubleArray;

//------------------------------------------------------------------
int allocate_double_array(doubleArray* newArray_p, int inSize)
{
    double* array_p;
    size_t totalSize;

    totalSize = sizeof(double) * inSize;
    array_p = (double*)malloc(totalSize);
    memset(array_p, 0, totalSize);
    if(array_p)
    {
        newArray_p->m_elements = array_p;
        newArray_p->m_size = inSize;
        return 1;
    }
    return 0;
}
//------------------------------------------------------------------
int init_array(doubleArray* inArray_p, doubleArray* xArray_p, double (*start_f)(double x))
{
    int i;
    if(inArray_p)
    {
        for(i = 0; i < inArray_p->m_size; i++)
        {
            inArray_p->m_elements[i] = (*start_f)(xArray_p->m_elements[i]);
        }
        return 1;
    }
    return 0;
}
//------------------------------------------------------------------
int init_x_arr(doubleArray* xArray, double max, double min)
{
    int i;
    double dx;
    dx = (max - min)/(xArray->m_size);
    for(i = 0; i < xArray->m_size; i++)
    {
        xArray->m_elements[i] = min + i * dx;
    }
    return 1;
}
//------------------------------------------------------------------
double value_at(doubleArray* array_p, int index)
{
    return array_p->m_elements[periodic_index_at(array_p->m_size, index)];
}
//------------------------------------------------------------------
void set_value_at(doubleArray* array_p, int index, double value)
{
    array_p->m_elements[periodic_index_at(array_p->m_size, index)] = value;
}
//------------------------------------------------------------------
int periodic_index_at(int arrSize, int index)
{
    if(0 < index && index < arrSize)
    {
        return index;
    }
    else if ( index < 0 && fabs(index) < arrSize)
    {
        return index + arrSize;
    }
    else if ( arrSize < index && index < 2 * arrSize)
    {
        return index - arrSize;
    }
    return 0;
}
//------------------------------------------------------------------
//******************************************************************
//------------------------------------------------------------------
double singe_time_step(double f, double dfdt, double dt)
{
    return f + dfdt * dt;
}
//------------------------------------------------------------------
int do_time_step(doubleArray* valArray_p, doubleArray* dArray_p, double dt)
{
    int i;
    double dval, fval, newval;
    if (valArray_p->m_size != dArray_p->m_size)
    {
        return 0;
    }
    for(i = 0; i < valArray_p->m_size; i++)
    {
        fval = value_at(valArray_p, i);
        dval = value_at(dArray_p, i);
        newval = singe_time_step(fval, dval, dt);
        set_value_at(valArray_p, i, newval);
        //valArray_p->m_elements[i] = singe_time_step(valArray_p->m_elements[i],
        //                                            dArray_p->m_elements[i],
        //                                            dt);
    }
    return 1;
}
//------------------------------------------------------------------
double central_second_deriv(double fdown, double f, double fup, double dx)
{
    return (fup-2*f+fdown)/(dx*dx);
}
//------------------------------------------------------------------

void laplace(doubleArray* array_p, doubleArray* toArray_p, double dx)
{
    double fdown, f, fup, newval;
    int i, arraySize;
    arraySize = array_p->m_size;

    for(i = 0; i < arraySize; i++)
    {
        fdown = value_at(array_p, i - 1);
        f = value_at(array_p, i);
        fup = value_at(array_p, i + 1);
        newval = central_second_deriv(fdown, f, fup, dx);
        set_value_at(toArray_p, i, newval);
    }
}
//------------------------------------------------------------------
void cube(doubleArray* array_p, doubleArray* toArray_p)
{
    int i;
    double val, cv;

    for(i = 0; i < array_p->m_size ; i++)
    {
        val = value_at(array_p, i);
        cv = val*val*val;
        set_value_at(toArray_p, i, cv);
    }
}
//------------------------------------------------------------------
void multiply_scalar(double s, doubleArray* toArray_p)
{
    int i;
    double val;

    for(i = 0; i < toArray_p->m_size ; i++)
    {
        val = s * value_at(toArray_p, i);
        set_value_at(toArray_p, i, val);
    }
}
//------------------------------------------------------------------
void sum_array(doubleArray* array_p, doubleArray* toArray_p)
{
    int i;
    double val_one, val_two, sum;

    for(i = 0; i < array_p->m_size ; i++)
    {
        val_one = value_at(array_p, i);
        val_two = value_at(toArray_p, i);
        sum = val_one + val_two;
        set_value_at(toArray_p, i, sum);
    }
}
//------------------------------------------------------------------
int do_diffusion(doubleArray* valArray_p,
                doubleArray* dArray_p,
                double TMAX,
                double dt,
                double dx)
{
    double tTotal = 0;
    while(tTotal < TMAX)
    {
        laplace(valArray_p, dArray_p, dx);
        do_time_step(valArray_p, dArray_p, dt);
        tTotal += dt;
    }
    return 1;
}
//------------------------------------------------------------------
int do_cahn_hilliard(doubleArray* valArray_p,
                    doubleArray* dArray_p,
                    doubleArray* cubeArray_p,
                    doubleArray* lptermArray_p,
                    doubleArray* sumArray_p,
                    double TMAX,
                    double dt,
                    double dx)
{
    double tTotal = 0;
    while(tTotal < TMAX)
    {
        laplace(valArray_p, lptermArray_p, dx);
        cube(valArray_p, cubeArray_p);
        multiply_scalar(-1, cubeArray_p);

        sum_array(valArray_p, sumArray_p);
        sum_array(lptermArray_p, sumArray_p);
        sum_array(cubeArray_p, sumArray_p);

        laplace(sumArray_p, dArray_p, dx);
        multiply_scalar(-1, dArray_p);

        do_time_step(valArray_p, dArray_p, dt);
        tTotal += dt;
    }
    return 1;
}
//------------------------------------------------------------------
int do_phase_field_crystal(doubleArray* valArray_p,
                            doubleArray* dArray_p,
                            doubleArray* cubeArray_p,
                            doubleArray* lapOneArray_p,
                            doubleArray* lapTwoArray_p,
                            doubleArray* sumArray_p,
                            double epsilon,
                            double TMAX,
                            double dt,
                            double dx)
{
    double tTotal = 0;
    while(tTotal < TMAX)
    {
        laplace(valArray_p, lapOneArray_p, dx);
        laplace(lapOneArray_p, lapTwoArray_p, dx);
        cube(valArray_p, cubeArray_p);

        multiply_scalar((1 - epsilon), valArray_p);
        multiply_scalar(2, lapOneArray_p);

        sum_array(valArray_p, sumArray_p);
        sum_array(lapOneArray_p, sumArray_p);
        sum_array(lapTwoArray_p, sumArray_p);
        sum_array(cubeArray_p, sumArray_p);

        laplace(sumArray_p, dArray_p, dx);

        do_time_step(valArray_p, dArray_p, dt);
        tTotal += dt;
    }
    return 1;
}
//------------------------------------------------------------------
void init_grid(doubleArray* valArray_p,
               doubleArray* xArray_p,
               double (*f_p)(double x),
               int gSize,
               int xMax,
               int xMin)
{


    init_x_arr(xArray_p, xMax, xMin);
    init_array(valArray_p, xArray_p, (*f_p));
}
//------------------------------------------------------------------
void copy_result(doubleArray* valArray_p, int gSize, double* out_p)
{
    int i;
    for (i = 0; i < gSize; i++)
    {
        out_p[i] = value_at(valArray_p, i);
    }
}
//------------------------------------------------------------------
double calc_dx(doubleArray* xArray_p)
{
    return value_at(xArray_p, 1) - value_at(xArray_p, 0);
}
//------------------------------------------------------------------
simpleExport int diffusion(double (*start_val_f)(double x),
                           int gSize,
                           double xMax,
                           double xMin,
                           double dt,
                           double TMAX,
                           double* out_p)
{
    int i;
    double dx;

    doubleArray valArray, xArray, dArray;
    allocate_double_array(&valArray, gSize);
    allocate_double_array(&xArray, gSize);
    allocate_double_array(&dArray, gSize);

    init_grid(&valArray, &xArray, (*start_val_f), gSize, xMax, xMin);
    dx = calc_dx(&xArray);
    do_diffusion(&valArray, &dArray, TMAX, dt, dx);
    copy_result(&valArray, gSize, out_p);

    return 7;
}
//------------------------------------------------------------------
simpleExport int xvalues(double (*start_val_f)(double x),
                           int gSize,
                           double xMax,
                           double xMin,
                           double* out_p)
{

    doubleArray valArray, xArray;
    allocate_double_array(&valArray, gSize);
    allocate_double_array(&xArray, gSize);
    init_grid(&valArray, &xArray, (*start_val_f), gSize, xMax, xMin);
    copy_result(&xArray, gSize, out_p);

    return 7;
}
//------------------------------------------------------------------
simpleExport int cahnhilliard(double (*start_val_f)(double x),
                               int gSize,
                               double xMax,
                               double xMin,
                               double dt,
                               double TMAX,
                               double* out_p)
{
    int i;
    double dx;

    doubleArray valArray, xArray, dArray, cubeArray, lptermArray, sumArray;

    allocate_double_array(&valArray, gSize);
    allocate_double_array(&xArray, gSize);
    allocate_double_array(&dArray, gSize);
    allocate_double_array(&cubeArray, gSize);
    allocate_double_array(&lptermArray, gSize);
    allocate_double_array(&sumArray, gSize);

    init_grid(&valArray, &xArray, (*start_val_f), gSize, xMax, xMin);
    dx = calc_dx(&xArray);
    do_cahn_hilliard(&valArray,
                     &dArray,
                     &cubeArray,
                     &lptermArray,
                     &sumArray,
                     TMAX,
                     dt,
                     dx);
    copy_result(&valArray, gSize, out_p);

    return 7;
}
//------------------------------------------------------------------
simpleExport int phasefieldcrystal(double (*start_val_f)(double x),
                                    int gSize,
                                    double xMax,
                                    double xMin,
                                    double dt,
                                    double TMAX,
                                    double epsilon,
                                    double* out_p)
{
    int i;
    double dx;

    doubleArray valArray, xArray, dArray, cubeArray, lapOneArray, lapTwoArray, sumArray;

    allocate_double_array(&valArray, gSize);
    allocate_double_array(&xArray, gSize);
    allocate_double_array(&dArray, gSize);
    allocate_double_array(&cubeArray, gSize);
    allocate_double_array(&lapOneArray, gSize);
    allocate_double_array(&lapTwoArray, gSize);
    allocate_double_array(&sumArray, gSize);

    init_grid(&valArray, &xArray, (*start_val_f), gSize, xMax, xMin);
    dx = calc_dx(&xArray);
    do_phase_field_crystal(&valArray,
                             &dArray,
                             &cubeArray,
                             &lapOneArray,
                             &lapTwoArray,
                             &sumArray,
                             epsilon,
                             TMAX,
                             dt,
                             dx);
    copy_result(&valArray, gSize, out_p);

    return 7;
}