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
	complex<double> *x, *y;
	double **output;
	FILE * pFile;
	x = (complex<double>*) malloc(sizeof(complex<double>)*n);
	y = (complex<double>*) malloc(sizeof(complex<double>)*n);
	
	int x_index = 0;
	for(int i = 0; i<n; i++) {
		x[i] = complex<double>(cos(2*M_PI*150*i/n), 0.0);
		if (i < 20) { y[i] = complex<double>(0.0,0.0)*std::polar(1.0, 1*M_PI*30*i/n);}
		else { y[i] = x[x_index]*std::polar(1.0, 1*M_PI*30*i/n); x_index++;}
	}


	pFile = fopen("resampleddata.m", "w");
	fprintf(pFile, "clear all; clc; \ndata=[");
	output =  (double**) malloc(sizeof(double*)*range);	
	for(int k=0; k<range; k++) {
		output[k] = (double*) malloc(sizeof(double)*nfft);
		output[k] = amb.ambigutyFunction(x, y, k, n, nfft);		
		for(int i=0; i<nfft; i++) { fprintf(pFile,"%3.5f \t", (output[k][i]));}
		fprintf(pFile, "\n");	
	}
	fprintf(pFile, "]; figure; surf(data); \n");
	fclose(pFile);
	free(output);


	//int	n = 8;
	//arrayMath arrayOp = arrayMath();
	//fftw_plan fft;

	// take fft of data


	//free(x);
	//free(absIn);
	//fftw_destroy_plan(fft);
	//fftw_cleanup();
	return 0;


}

