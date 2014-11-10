// Radar.cpp : Defines the entry point for the console application.
//


#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
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
#include <iostream>
#include <fstream>

using namespace std;


int _tmain(int argc, _TCHAR* argv[])
{

	int	n = 512;
	arrayOperations arrayOp = arrayOperations();
	fftw_plan fft;

	complex<double> *x;
	fftw_complex *in, *out;
	double *absIn;
	x = (complex<double>*) malloc(sizeof(complex<double>) * n);
	
	
	

	// generate data
	for(int i=0; i<n; i++) x[i] = complex<double> (cos(2*M_PI*12*i/n), sin(2*M_PI*12*i/n));
	FILE * pFile;
	pFile = fopen("fft.txt", "w");
	fprintf(pFile, "clear all; close all; clc; \ndata=[");
	for(int i=0; i<n; i++) fprintf(pFile,"%3.5f + %3.5f*i\n", real(x[i]), imag(x[i]));
	fprintf(pFile, "]; \n");

	// shift data by k
	double factor = 0.5;
	complex<double> *resampled;
	resampled = (complex<double>*) malloc(sizeof(complex<double>) * (int) (ceil(n * factor)));
	resampled = arrayOp.resample(x, n, factor);
	fprintf(pFile, "\nresampled_data=[");
	for(int i=0; i<(int) (ceil(n * factor)); i++) {fprintf(pFile, "%3.5f + %3.5f*i \n", real(resampled[i]), imag(resampled[i]));}
	fprintf(pFile, "];\n");
	fclose(pFile);	

	free(resampled);
	free(x);
	fftw_cleanup();

	//// take fft of data
	//in = (fftw_complex*) malloc(sizeof(fftw_complex) * (int) (ceil(n * factor)));
	//in = reinterpret_cast<fftw_complex*>(resampled);
	//absIn = (double*) malloc(sizeof(double) * (int) (ceil(n * factor)));
	//out = (fftw_complex*) malloc(sizeof(fftw_complex) * (int) (ceil(n * factor)));
	//fft = fftw_plan_dft_1d((int) (ceil(n * factor)), in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	//fftw_execute(fft);
	//fprintf(pFile, "\n M=[");
	//for(int i=0; i<(int) (ceil(n * factor)); i++) {fprintf(pFile, "%.5f + %.5f*i \n", out[i][0], out[i][1]);}
	//fprintf(pFile, "];\n");
	//fprintf(pFile, "plot([0:%i], abs(M));", n-1);
	//fclose(pFile);
	//// take absolute value of fftw_complex data type
	//absIn = arrayOp.complexAbs(out, (int) (ceil(n * factor)));
	//for(int i=0; i<n; i++) {cout << absIn[i] << "\n";}
	//cout << "End of absolute value of FFT of data.\n";


	//free(x);
	//free(absIn);
	//fftw_destroy_plan(fft);
	//fftw_cleanup();
	return 0;


}

