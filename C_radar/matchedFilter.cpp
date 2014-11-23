#include "stdafx.h"
#include <iostream>
#include "arrayMath.h"
#include "matchedFilter.h"
#include <complex>
#include "fftw3.h"	
using namespace std;

arrayMath arrayOp;
fftw_plan fft;

matchedFilter::matchedFilter(void)
{
	arrayOp = arrayMath();
}
matchedFilter::~matchedFilter(void)
{
}

complex<double> * matchedFilter::ambigutyFunction(complex<double> * referenceSignal, complex<double> * surveillenceSignal, int k, int dataLength, int nfft)
{
	// arrayMath arrayOp = arrayMath();
	complex<double> *shiftedReferenceSignal, *correlator, *correlatorResampled;
	fftw_complex *fftBlock;

	shiftedReferenceSignal = (complex<double>*) malloc(sizeof(complex<double>)*dataLength);
	correlator = (complex<double>*) malloc(sizeof(complex<double>)*dataLength);
	fftBlock = (fftw_complex*) malloc(sizeof(fftw_complex)*dataLength);


	shiftedReferenceSignal = arrayOp.shiftArray(referenceSignal, dataLength, k);
	correlator = arrayOp.elementMultiplication(shiftedReferenceSignal, surveillenceSignal, dataLength, "conj");
	
	double factor = 1/((double)dataLength/(double)nfft);
	correlatorResampled = (complex<double>*) malloc(sizeof(complex<double>)* (int) (ceil(dataLength * factor)));
	correlatorResampled = arrayOp.resample(correlator, dataLength, factor);

	return correlatorResampled;
}
