#include "eval_circle.h"

int main(void)
{

#define size 120

	double ArrayOfXY[size];
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
	double f[size / 3];
	int counter = 0;
	for (int i = 0; i < size; i += 3) {
		f[counter] = ArrayOfXY[i];
		x[counter] = ArrayOfXY[i + 1];
		y[counter] = ArrayOfXY[i + 2];
		counter++;
	}

	double lse;
	int niter;
	Circle crcl;

	estimate_center_of_circle(x, y, counter, &lse, &niter, &crcl);

	return 0;
}
