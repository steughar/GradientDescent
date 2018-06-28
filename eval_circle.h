#ifndef EVAL_CIRCLE_H_
#define EVAL_CIRCLE_H_

#include <stddef.h>
#include <math.h>

typedef struct 
{
  double CenterX;
  double CenterY;
  double Radius;
  double phi;     //--> Rotation angle for experimental dots adjusting
} Circle;

typedef struct 
{
  double *x;
  double *y;
  double *f0;
  size_t n;
} ExperimentalDots;

typedef struct
{
  double Q;
  double f;
} Weights;


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
// compute partial derivative by Q
double gradQ
(
    double Q,              //--> Current value of a variable on this particular step
    double f,              //--> Current value of b variable on this particular step
    ExperimentalDots *Dots,//--> Experimental data wrapper
    Circle *crcl           //--> Circle structure with approximation circle
 );
// compute partial derivative
double gradF
(
    double Q,              //--> Current value of a variable on this particular step
    double f,              //--> Current value of b variable on this particular step
    ExperimentalDots *Dots,//--> Experimental data wrapper
    Circle *crcl           //--> Circle structure with approximation circle
 );
// compute weighting function
double WeightingFunction
(
    double Q,              //--> Current Q
    double f,              //--> Current f
    ExperimentalDots *Dots,//--> Experimental dots wrapper structure
    Circle *crcl           //--> Approximation circle
 );
// find alpha
double GetAlpha
(
    double Q,              //--> Current Q
    double f,              //--> Current f
    double f0              //--> Current f0
 );




// do estimation of circle that approximate experimental data
void
estimate_center_of_circle
(
    ExperimentalDots *Dots, //--> Struct that contains experimental data
    size_t n,              //--> Number of experimental dots (input)
    double *lse,           //--> Least Squares Error estimation after approximation (output)
    int *niter,            //--> Number of iterations (output) 
    Circle *crcl           //--> In this structure result will be written (result)
 );

// take experimental dots to the origin and rotate the circle 
void
adjust_experimental_data
(
    ExperimentalDots *Dots, //--> Experimental dots wrapper struct
    Circle *crcl            //--> Approximation circle structure
 );

void find_Q_factor
(
    ExperimentalDots *Dots, //--> Experimental dots wrapper struct (input)
    Weights *weights,       //--> seed value must be passed on call, result will overwrite this values 
    Circle *crcl,           //--> Approximation circle
    double *lse           //--> Optimization error (LSE)
 );

#endif






