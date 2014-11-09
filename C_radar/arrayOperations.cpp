#include "stdafx.h"
#include <iostream>
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

complex<double> * arrayOperations::resample(complex<double> *x, int dataLength, double factor) 
{
	float a_real, a_imag, b_real, b_imag;
	complex<double> *resampled;
	dataLength = (int) floor(dataLength*factor);
	resampled = (complex<double>*) malloc(sizeof(complex<double>) * dataLength);
	if (factor == 1) {resampled = x;}
	if (factor > 1) {
		for(int i=0; i<dataLength; i++) {
			double p=i*(1.0/factor);
			int temp = (int) floor(p);
			float s = p-temp;
			if(temp < dataLength) {
				a_real = real(x[temp]);
				a_imag = imag(x[temp]);
			} else { a_real = 0;  a_imag = 0;}
			if(temp < dataLength-1) {
				b_real = real(x[temp+1]);
				b_imag = imag(x[temp+1]);
			} else { b_real = 0;  b_imag = 0;}
			resampled[i] = complex<double> ((1-s)*a_real + s*b_real, (1-s)*a_imag + s*b_imag);
		} 
	} else {
		for(int i=0; i<dataLength; i++) {
			double p1=i*(1.0/factor);
			double p2=(i+1)*(1.0/factor);
			if((p1 < dataLength) && ((int) floor(p1) != (int) ceil(p1))) {
				resampled[i] += (((int)ceil(p1)) - p1)*x[(int)floor(p1)];
			}
			for (int j=(int) ceil(p1); j<min((int)floor(p2), dataLength); j++) {
				resampled[i] += x[j];
			}
    		if (((int)floor(p2) < dataLength) && ((int)floor(p2) != (int)ceil(p2))) {
    			resampled[i] += (p2-((int)floor(p2)))*x[(int)floor(p2)];
    		}
			resampled[i] = complex<double> (real(resampled[i])*factor, imag(resampled[i])*factor);		
		}
		
	}
	return resampled;
}
	




