#ifndef EVAL_CIRCLE_H_
#define EVAL_CIRCLE_H_

#include <stddef.h>
#include <math.h>
#include <stdio.h>

typedef struct 
{
        double CenterX;
        double CenterY;
        double Radius;
} Circle;

typedef struct 
{
        double *x;
        double *y;
        size_t n;
} ExperimentalDots;

// compute partial derrivative by a
double gradA
    (
        double a,              //--> Current value of a variable on this particular step
        double b,              //--> Current value of b variable on this particular step
        double R,              //--> Current value of R variable on this particular step
        ExperimentalDots *Dots //--> Experimental data wrapper
        );
// compute partial derrivative bu b
double gradB
    (
        double a,              //--> Current value of a variable on this particular step
        double b,              //--> Current value of b variable on this particular step
        double R,              //--> Current value of R variable on this particular step
        ExperimentalDots *Dots //--> Experimental data wrapper
        );
// compute partial derrivative by R
double gradR
    (
        double a,              //--> Current value of a variable on this particular step
        double b,              //--> Current value of b variable on this particular step
        double R,              //--> Current value of R variable on this particular step
        ExperimentalDots *Dots //--> Experimental data wrapper
    );
// do estimation of circle that approximate experimental data
void estimate_center_of_circle
(
        const double x[],      //--> Array contains X coordinates of experimenal dots (input)
        const double y[],      //--> Array contains Y coordinates of experimenal dots (input)
        size_t n,              //--> Number of experimental dots (input)
        double *lse,           //--> Least Squares Error estimation after approximation (output)
        int *niter,            //--> Number of iterations (output) 
        Circle *crcl           //--> In this structure result will be written (result)
);

#endif






