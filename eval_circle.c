#include "eval_circle.h"

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


void estimate_center_of_circle
(
 const double x[],   //--> Array contains X coordinates of experimenal dots
 const double y[],   //--> Array contains Y coordinates of experimenal dots
 size_t n,           //--> Number of experimental dots
 double *lse,        //--> Least Squares Error estimation after approximation
 int *niter,         //--> Number of iterations          
 Circle *crcl        //--> In this structure result will be written
 )
{
  ExperimentalDots Dots;
  Dots.x = x;
  Dots.y = y;
  Dots.n = n;

  double error = 0;
  double sigma = 0.02;
  double minimal_error = 999999999999;

  //-------------Preliminary estimation of seed values for algorithm-------------
  double AverageA = 0;
  double AverageB = 0;
  double AverageR = 0;

  for (int i = 0; i < n; i++)
    {
      AverageA += x[i];
      AverageB += y[i];

    }
  AverageA /= n;
  AverageB /= n;

  for (int i = 0; i < n; i++)
    {
      AverageR += sqrt((AverageA - Dots.x[i])*(AverageA - Dots.x[i]) +
                       (AverageB - Dots.y[i])*(AverageB - Dots.y[i]));
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
          double x_i = Dots.x[i];
          double y_i = Dots.y[i];

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
      a = a - sigma * (gradA(a, b, R, &Dots) / n);
      b = b - sigma * (gradB(a, b, R, &Dots) / n);
      R = R - sigma * (gradR(a, b, R, &Dots) / n);
    }
  //--------------------------------------------------------------------------------

  //Write final coefficients of weighting function into passed structure
  crcl->CenterX = a_bestFit;
  crcl->CenterY = b_bestFit;
  crcl->Radius  = R_bestFit;

  *lse = minimal_error / n;
  *niter = 1000;
}
