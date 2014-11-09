#include <complex>
#include "math.h"
#include <cmath>
#include "fftw3.h"
using namespace std;
#pragma once

class arrayOperations
{
private:

public:
	arrayOperations(void);
	~arrayOperations(void);

	complex<double> * shiftArray(complex<double> *x, int dataLength, int k);
	double * complexAbs(fftw_complex *x, int dataLength);

};

