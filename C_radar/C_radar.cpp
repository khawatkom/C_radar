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

using namespace std;


int _tmain(int argc, _TCHAR* argv[])
{
	complex<double> *x;
	fftw_complex *in, *out;
	int N = (int) pow(2.0, 18.0);
	x = (complex<double>*) malloc(sizeof(complex<double>) * N);

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	in = reinterpret_cast<fftw_complex*>(x);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	fftw_plan p;
	int i;

	auto start = std::chrono::high_resolution_clock::now();
	/* prepare a cosine wave */
	for (i = 0; i < N; i++) {
		x[i] = complex<double>(cos(3*2*M_PI*i/N), sin(2*2*M_PI*i/N)); //cos(3 * 2*M_PI*i/N);		
	}
	auto elapsed = std::chrono::high_resolution_clock::now() - start;
	long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
	std::cout << "Data generation took: " << microseconds << " microseconds\n";
//	for (i = 0; i < N; i++)
//		printf("freq: %3d %+9.5f %+9.5f I\n", i, in[i][0], in[i][1]);

	p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);




	/* forward Fourier transform, save the result in 'out' */
	start = std::chrono::high_resolution_clock::now();
	for(int i=0; i<150; i++)
		fftw_execute(p);
	elapsed = std::chrono::high_resolution_clock::now() - start;
	microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
	printf("Elapsed time for %i point complex FFT is : ", N);
	std::cout <<microseconds << " microseconds\n";
	std::cout <<"One FFT is: " << microseconds/150 << " microseconds\n";
//	for (i = 0; i < N; i++)
//		printf("freq: %3d %+9.5f %+9.5f I\n", i, out[i][0], out[i][1]);
	fftw_destroy_plan(p);
	fftw_cleanup();
	return 0;


}

