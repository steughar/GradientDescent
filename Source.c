#include "eval_circle.h"
#include <stdio.h>

int main(void)
{

#define size 120

  double ArrayOfXY[size];
  ExperimentalDots Dots;
  FILE* stream;
        
  stream = fopen("W:\\work\\LKARD\\CEX\\CEX\\40_dots.txt", "r");
  int fileCarret = fscanf(stream, "%lf", &ArrayOfXY[0]);
  int j = 1;
  
  while (fileCarret != EOF)
  {
    fileCarret = fscanf(stream, "%lf", &ArrayOfXY[j]);
    j++;
    
  }
  
  double x[size / 3];
  double y[size / 3];
  double f0[size / 3];
  int counter = 0;
  for (int i = 0; i < size; i += 3) {
    f0[counter] = ArrayOfXY[i];
    x[counter] = ArrayOfXY[i + 1];
    y[counter] = ArrayOfXY[i + 2];
    counter++;
  }
  
  double lse;
  int niter;
  Circle crcl;

  Dots.x = x;
  Dots.y = y;
  Dots.f0 = f0;
  Dots.n = counter;
  
  estimate_center_of_circle(&Dots, counter, &lse, &niter, &crcl);

  adjust_experimental_data(&Dots, &crcl);

  /*
  FILE* output;

  output = fopen("W:\\work\\LKARD\\CEX\\CEX\\40_dots_after_c.txt", "w");
  for(int i = 0; i < Dots.n; i++)
  {
    fprintf(output, "%f %f %f \n", Dots.x[i], Dots.y[i], Dots.f0[i]);
  }

  fclose(output);
  */
  Weights weights;
  // seed values for optimization
  weights.Q = 6;
  weights.f = 29000;

  
  find_Q_factor(&Dots, &weights, &crcl, &lse);
  
  return 0;
}
