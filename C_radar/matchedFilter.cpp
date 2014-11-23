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

double * matchedFilter::ambigutyFunction(complex<double> * referenceSignal, complex<double> * surveillenceSignal, int k, int dataLength, int nfft)
{
	// arrayMath arrayOp = arrayMath();
	complex<double> *shiftedReferenceSignal, *correlator, *correlatorResampled;
	double *absFFT;
	fftw_complex *input, *output;
	fftw_plan fftBlock;

	shiftedReferenceSignal = (complex<double>*) malloc(sizeof(complex<double>)*dataLength);
	correlator = (complex<double>*) malloc(sizeof(complex<double>)*dataLength);
	input = (fftw_complex*) malloc(sizeof(fftw_complex)*nfft);
	output = (fftw_complex*) malloc(sizeof(fftw_complex)*nfft);
	absFFT = (double*) malloc(sizeof(double)*nfft);


	shiftedReferenceSignal = arrayOp.shiftArray(referenceSignal, dataLength, k);
	correlator = arrayOp.elementMultiplication(shiftedReferenceSignal, surveillenceSignal, dataLength, "conj");
	
	double factor = 1/((double)dataLength/(double)nfft);
	correlatorResampled = (complex<double>*) malloc(sizeof(complex<double>)* (int) (ceil(dataLength * factor)));
	correlatorResampled = arrayOp.resample(correlator, dataLength, factor);
	input = reinterpret_cast<fftw_complex*>(correlatorResampled);
	fftBlock = fftw_plan_dft_1d(nfft, input, output, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(fftBlock);
	absFFT = arrayOp.complexAbs(output, nfft*2);



	return absFFT;



}
