#include <complex>
#include "math.h"
#include <cmath>
#include "fftw3.h"
using namespace std;
#pragma once

class arrayMath
{
private:

public:
	arrayMath(void);
	~arrayMath(void);


	complex<double> * shiftArray(complex<double> *x, int dataLength, int k);
	double * complexAbs(fftw_complex *x, int dataLength);
	complex<double> * resample(complex<double> *x, int dataLength, double factor);
	complex<double> * elementMultiplication(complex<double> *x, complex<double> *y, int dataLength, string type);
};

