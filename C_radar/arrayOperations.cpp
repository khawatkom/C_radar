#include "stdafx.h"
#include "arrayOperations.h"
#include <complex>
#include "fftw3.h"
using namespace std;

arrayOperations::arrayOperations(void)
{
}
arrayOperations::~arrayOperations(void)
{
}

complex<double> * arrayOperations::shiftArray(complex<double>* x, int dataLength, int k)
{ 
	int shiftIndex = 0;
	complex<double> *shifted;
	shifted = (complex<double>*) malloc(sizeof(complex<double>) * dataLength);
	for (int i=0; i<dataLength; i++) {
		if (i < k) {shifted[i] = complex<double>(0.0, 0.0); }
		else { shifted[i] = x[shiftIndex]; shiftIndex++;}
	}
	return shifted;
}

double * arrayOperations::complexAbs(fftw_complex *x, int dataLength)
{
	fftw_complex *in = (fftw_complex*) malloc(sizeof(fftw_complex) * dataLength);
	in = x;
	double *abs;
	abs = (double*) malloc(sizeof(double) * dataLength);
	for(int i=0; i<dataLength; i++) {abs[i] = (double) sqrt(pow(in[i][0], 2.0) + pow(in[i][1], 2.0));}
	return abs;
}



