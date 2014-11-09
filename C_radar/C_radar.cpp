// Radar.cpp : Defines the entry point for the console application.
//



#include "stdafx.h"
#include "Windows.h"
#include <time.h>
#include <math.h>
#define _USE_MATH_DEFINES
#include <complex>
#include <cmath>
#include "fftw3.h"
#include <ctime>
#include <iostream>
#include <chrono>
#include "arrayOperations.h"

using namespace std;


int _tmain(int argc, _TCHAR* argv[])
{

	int	n = 8;
	arrayOperations arrayOp = arrayOperations();
	fftw_plan fft;

	complex<double> *x;
	fftw_complex *in, *out;
	double *absIn;
	x = (complex<double>*) malloc(sizeof(complex<double>) * n);
	in = (fftw_complex*) malloc(sizeof(fftw_complex) * n);
	out = (fftw_complex*) malloc(sizeof(fftw_complex) * n);
	absIn = (double*) malloc(sizeof(double) * n);

	// generate data
	for(int i=0; i<n; i++) x[i] = complex<double> (i, 0.0);
	for(int i=0; i<n; i++) cout << x[i] << "\n";
	cout << "End of actual data.\n";

	// shift data by k
	double factor = 1.5;
	complex<double> *resampled;
	resampled = (complex<double>*) malloc(sizeof(complex<double>) * (int) ceil(n * factor));
	resampled = arrayOp.resample(x, n, factor);
	for(int i=0; i<ceil(n*factor); i++) cout << resampled[i] << "\n";
	cout << "End of " << factor << " factorized data.\n";
	
	//// take fft of data
	//in = reinterpret_cast<fftw_complex*>(x);
	//fft = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	//fftw_execute(fft);
	//for(int i=0; i<n; i++) {printf("freq %i: %f + %f*i \n", i, out[i][0], out[i][1]);}
	//cout << "End of FFT of data.\n";



	//// take absolute value of fftw_complex data type
	//absIn = arrayOp.complexAbs(out, n);
	//for(int i=0; i<n; i++) {cout << absIn[i] << "\n";}
	//cout << "End of absolute value of FFT of data.\n";


	//free(x);
	//free(absIn);
	//fftw_destroy_plan(fft);
	//fftw_cleanup();
	return 0;


}

