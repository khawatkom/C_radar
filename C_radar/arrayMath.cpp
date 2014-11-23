#include "stdafx.h"
#include <iostream>
#include "arrayMath.h"
#include <complex>
#include "fftw3.h"
#include <stdio.h>
#include <string.h>
using namespace std;

arrayMath::arrayMath(void)
{
}
arrayMath::~arrayMath(void)
{
}

complex<double> * arrayMath::shiftArray(complex<double>* x, int dataLength, int k)
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

double * arrayMath::complexAbs(fftw_complex *x, int dataLength)
{
	fftw_complex *in = (fftw_complex*) malloc(sizeof(fftw_complex) * dataLength);
	in = x;
	double *abs;
	abs = (double*) malloc(sizeof(double) * dataLength);
	for(int i=0; i<dataLength; i++) {abs[i] = (double) sqrt(pow(in[i][0], 2.0) + pow(in[i][1], 2.0));}
	return abs;
}

complex<double> * arrayMath::resample(complex<double> *x, int dataLength, double factor) 
{
	double a_real, a_imag, b_real, b_imag;
	complex<double> *resampled;
	int newDataLength = (int) ceil(dataLength*factor);
	resampled = (complex<double>*) malloc(sizeof(complex<double>) * newDataLength);
	if (factor == 1) {for(int i=0; i<dataLength; i++) resampled[i] = complex<double> (real(x[i]), imag(x[i]));}
	if (factor > 1) {
		for(int i=0; i<newDataLength; i++) {
			double p=i*(1.0/factor);
			int temp = (int) floor(p);
			float s = p-temp;
			if(temp < newDataLength) {
				a_real = real(x[temp]);
				a_imag = imag(x[temp]);
			} else { a_real = 0;  a_imag = 0;}
			if(temp < newDataLength-1) {
				b_real = real(x[temp+1]);
				b_imag = imag(x[temp+1]);
			} else { b_real = 0;  b_imag = 0;}
			resampled[i] = complex<double> ((1-s)*a_real + s*b_real, (1-s)*a_imag + s*b_imag);
			if (real(resampled[i]) < 10e-20 && real(resampled[i]) > -10e-20) {resampled[i] = complex<double> (0, imag(resampled[i]));} 
			if (imag(resampled[i]) < 10e-20 && imag(resampled[i]) > -10e-20) {resampled[i] = complex<double> (real(resampled[i]), 0);}
		} 
	} else {
		for(int i=0; i<newDataLength; i++) {
			double temp_real = 0, temp_imag = 0;
			double p1=i*(1.0/factor);
			double p2=(i+1)*(1.0/factor);
			if((p1 < dataLength) && ((int) floor(p1) != (int) ceil(p1))) {
				temp_real +=  (((int)ceil(p1)) - p1)*real(x[(int)floor(p1)]);
				temp_imag +=  (((int)ceil(p1)) - p1)*imag(x[(int)floor(p1)]);
				resampled[i] = complex<double> (temp_real, temp_imag);
			}
			for (int j=(int) ceil(p1); j<min((int)floor(p2), dataLength); j++) {
				temp_real += real(x[j]);
				temp_imag += imag(x[j]);
				resampled[i] = complex<double> (temp_real, temp_imag);
			}
    		if (((int)floor(p2) < dataLength) && ((int)floor(p2) != (int)ceil(p2))) {
				temp_real +=  (p2-((int)floor(p2)))*real(x[(int)floor(p2)]);
				temp_imag +=  (p2-((int)floor(p2)))*imag(x[(int)floor(p2)]);
				resampled[i] = complex<double> (temp_real, temp_imag);
    		}
			resampled[i] = complex<double> (real(resampled[i])*factor, imag(resampled[i])*factor);
			if (real(resampled[i]) < 10e-20 && real(resampled[i]) > -10e-20) {resampled[i] = complex<double> (0, imag(resampled[i]));} 
			if (imag(resampled[i]) < 10e-20 && imag(resampled[i]) > -10e-20) {resampled[i] = complex<double> (real(resampled[i]), 0);}
		}
		
	}
	return resampled;
}


complex<double> * arrayMath::elementMultiplication(complex<double> *x, complex<double> *y, int dataLength, string type) 
{
	complex<double> * result;
	result = (complex<double>*) malloc(sizeof(complex<double>)*dataLength);

	if (type == "conj") {
		for(int i=0; i<dataLength; i++)
			result[i] = x[i]*conj(y[i]);		
	}else {
		for(int i=0; i<dataLength; i++)
			result[i] = x[i]*y[i];		
	}
	return result;
}
