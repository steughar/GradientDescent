#include "eval_circle.h"
#include <stdlib.h>

double gradA(double a, double b, double R, ExperimentalDots *Dots)
{
  double Result = 0;

  for (int i = 0; i < Dots->n; i++)
    {
      double x = Dots->x[i];
      double y = Dots->y[i];
      double Numerator = -2 * (x - a) * (sqrt((x - a)*(x - a) + (y - b)*(y - b)) - R);
      double Denominator = sqrt((x - a)*(x - a) + (y - b)*(y - b));
      Result += Numerator / Denominator;
    }
  
  return(Result);
}

double gradB(double a, double b, double R, ExperimentalDots *Dots)
{
  double Result = 0;

  for (int i = 0; i < Dots->n; i++)
    {
      double x = Dots->x[i];
      double y = Dots->y[i];
      double Numerator = -2 * (y - b) * (sqrt((x - a)*(x - a) + (y - b)*(y - b)) - R);
      double Denominator = sqrt((x - a)*(x - a) + (y - b)*(y - b));
      Result += Numerator / Denominator;
    }

  return(Result);
}

double gradR(double a, double b, double R, ExperimentalDots *Dots)
{
  double Result = 0;
  for (int i = 0; i < Dots->n; i++)
  {
    double x = Dots->x[i];
    double y = Dots->y[i];
    Result += -2 * (sqrt((x - a)*(x - a) + (y - b)*(y - b)) - R);
  }

  return(Result);
}

double gradF(double Q, double f, ExperimentalDots *Dots, Circle *crcl)
{
  double delta = 0.00001;
  double Result = (WeightingFunction(Q, f + delta, Dots, crcl) -
                   WeightingFunction(Q, f, Dots, crcl)) / delta;
  return(Result);
  
}

double gradQ(double Q, double f, ExperimentalDots *Dots, Circle *crcl)
{
  
  double delta = 0.00001;
  double Result = (WeightingFunction(Q + delta, f, Dots, crcl) -
                   WeightingFunction(Q, f, Dots, crcl)) / delta;
  return(Result);
}

double GetAlpha(double Q, double f, double f0)
{
  double Result = 2*atan(Q*(f / f0 - f0 / f));
  return(Result);
}

double WeightingFunction(double Q, double f, ExperimentalDots *Dots, Circle *crcl)
{
  
  double Result = 0;
  double alpha = 0;
  // dynamically declare an arrays size of n
  double *new_x = (double *)malloc(sizeof(double)*Dots->n);
  double *new_y = (double *)malloc(sizeof(double)*Dots->n);

  // generate new dots 
  for(int i = 0; i < Dots->n; i++) {
    
    alpha = GetAlpha(Q, f, Dots->f0[i]);
    new_x[i] = cos(alpha)*crcl->Radius;
    new_y[i] = sin(alpha)*crcl->Radius;
  }

  // compute error 
  for(int i = 0; i < Dots->n; i++) {
    
    Result += sqrt((Dots->x[i] - new_x[i])*(Dots->x[i] - new_x[i])
                   + (Dots->y[i] - new_y[i])*(Dots->y[i] - new_y[i]));
  }

  return(Result);
}


void estimate_center_of_circle
(
 ExperimentalDots *Dots,
 size_t n,             
 double *lse,          
 int *niter,         
 Circle *crcl        
 )
{

  double error = 0;
  double sigma = 0.02;
  double minimal_error = 999999999999;

  //-------------Preliminary estimation of seed values for algorithm-------------
  double AverageA = 0;
  double AverageB = 0;
  double AverageR = 0;

  for (int i = 0; i < n; i++)
    {
      AverageA += Dots->x[i];
      AverageB += Dots->y[i];

    }
  AverageA /= n;
  AverageB /= n;

  for (int i = 0; i < n; i++)
    {
      AverageR += sqrt((AverageA - Dots->x[i])*(AverageA - Dots->x[i]) +
                       (AverageB - Dots->y[i])*(AverageB - Dots->y[i]));
    }
  AverageR /= n;
  //------------------------------------------------------------------------------

        
  double a_bestFit = 0;
  double b_bestFit = 0;
  double R_bestFit = 0;

  // declaration of seed values
  double a = AverageA;
  double b = AverageB;
  double R = AverageR;

  //----------------Gradient Descent algorithm Count = 1000 sets of learning-------
  for (int Count = 0; Count < 1000; Count++)
    {
      error = 0;
      for (int i = 0; i < n; i++)
        {
          double x_i = Dots->x[i];
          double y_i = Dots->y[i];

          //weighting function
          error += (sqrt((x_i - a)*(x_i - a) + (y_i - b)*(y_i - b)) - R) *
            (sqrt((x_i - a)*(x_i - a) + (y_i - b)*(y_i - b)) - R);
        }
                
      if (error < minimal_error)
        {
          //update best fit result
          minimal_error = error;
          a_bestFit = a;
          b_bestFit = b;
          R_bestFit = R;
        }
      //updates weights a,b,R
      a = a - sigma * (gradA(a, b, R, Dots) / n);
      b = b - sigma * (gradB(a, b, R, Dots) / n);
      R = R - sigma * (gradR(a, b, R, Dots) / n);
    }
  //--------------------------------------------------------------------------------

  //Write final coefficients of weighting function into passed structure
  crcl->CenterX = a_bestFit;
  crcl->CenterY = b_bestFit;
  crcl->Radius  = R_bestFit;
  crcl->phi     = (AverageB - b_bestFit)/R_bestFit;

  *lse = minimal_error / n;
  *niter = 1000;
}


void
adjust_experimental_data (
    ExperimentalDots *Dots,
    Circle *crcl
)
{

  // take down experimental data to the origin
  for (int i = 0; i < Dots->n; i++) {
    Dots->x[i] -= crcl->CenterX;
    Dots->y[i] -= crcl->CenterY;
  }

  double cosPhi = cos(crcl->phi);
  double sinPhi = sin(crcl->phi);

  // rotate experimental data
  for (int i = 0; i < Dots->n; i++) {
    Dots->y[i] = Dots->y[i]*cosPhi + Dots->x[i]*sinPhi;
    Dots->x[i] = Dots->x[i]*cosPhi - Dots->y[i]*sinPhi;
  }

}

void find_Q_factor
(
    ExperimentalDots *Dots,
    Weights *weights,
    Circle *crcl,
    double *lse
)
{
  double sigma = 0.02;
  double error = 0;
  double minimal_error = 99999999;
  double Q_bestFit = 0;
  double f_bestFit = 0;
  double Q = weights->Q;
  double f = weights->f;
                 
  for(int count = 0; count < 1000; count++) {

    error = WeightingFunction(Q, f, Dots, crcl);
    
    if(error < minimal_error) {
      minimal_error = error;
      Q_bestFit = Q;
      f_bestFit = f;
    }
    
    Q -= sigma*gradQ(Q, f, Dots, crcl);
    f -= sigma*gradF(Q, f, Dots, crcl);
  }
  
  weights->Q = Q_bestFit;
  weights->f = f_bestFit;
  *lse = minimal_error;
}

    
