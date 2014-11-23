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
#include "arrayMath.h"
#include <iostream>
#include <fstream>
#include "matchedFilter.h"

using namespace std;


int _tmain(int argc, _TCHAR* argv[])
{
	arrayMath arrayOp = arrayMath();
	matchedFilter amb = matchedFilter();
	
	const int range = 50;
	int n = 10000;
	int nfft = 1024;
	complex<double> *x, *y, **output, *test;
	FILE * pFile;
	x = (complex<double>*) malloc(sizeof(complex<double>)*n);
	y = (complex<double>*) malloc(sizeof(complex<double>)*n);
	
	int x_index = 0;
	for(int i = 0; i<n; i++) {
		x[i] = complex<double>(cos(2*M_PI*150*i/n), 0.0);
		if (i < 20) { y[i] = complex<double>(0.0,0.0);}
		else { y[i] = x[x_index]; x_index++;}
	}


	pFile = fopen("resampleddata.m", "w");
	fprintf(pFile, "clear all; clc; \ndata=[");
	output =  (complex<double>**) malloc(sizeof(complex<double>*)*range);	
	for(int k=0; k<range; k++) {
		output[k] = (complex<double>*) malloc(sizeof(complex<double>)*nfft);
		output[k] = amb.ambigutyFunction(x, y, k, n, nfft);		
		for(int i=0; i<nfft; i++) { fprintf(pFile,"%3.5f + %3.5f*i \t", real(output[k][i]), imag(output[k][i]));}
		fprintf(pFile, "\n");	
	}
	fprintf(pFile, "]; figure; surf(data); \n");
	fclose(pFile);
	free(output);


	//int	n = 8;
	//arrayMath arrayOp = arrayMath();
	//fftw_plan fft;


	//// generate data
	//for(int i=0; i<n; i++) x[i] = complex<double> (cos(2*M_PI*12*i/n), 0);
	//FILE * pFile;
	//pFile = fopen("resampleddata.m", "w");
	//fprintf(pFile, "clear all; close all; clc; \ndata=[");
	//for(int i=0; i<n; i++) fprintf(pFile,"%3.5f + %3.5f*i\n", real(x[i]), imag(x[i]));
	//fprintf(pFile, "]; \n");

	//// shift data by k
	//double factor = 0.04;
	//complex<double> *resampled;
	//resampled = (complex<double>*) malloc(sizeof(complex<double>) * (int) (ceil(n * factor)));
	//resampled = arrayOp.resample(x, n, factor);
	//fprintf(pFile, "\nresampled_data=[");
	//for(int i=0; i<(int) (ceil(n * factor)); i++) {fprintf(pFile, "%3.5f + %3.5f*i \n", real(resampled[i]), imag(resampled[i]));}
	//fprintf(pFile, "];\n");
	//fclose(pFile);	


	


	//free(resampled);
	//free(x);
	//fftw_cleanup();

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

