#include "stdafx.h"
#include <iostream>
#include "arrayMath.h"
#include <complex>
#include "fftw3.h"	
using namespace std;
#pragma once

class matchedFilter
{
private:

public:
	matchedFilter(void);
	~matchedFilter(void);

double * ambigutyFunction(complex<double> * referenceSignal, complex<double> * surveillenceSignal, int k, int dataLength, int nfft);

};