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

	int	n = 4;
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
	for(int i=0; i<n; i++) x[i] = complex<double> (cos(2*M_PI*1*i/n), 0.0);
	for(int i=0; i<n; i++) cout << x[i] << "\n";
	cout << "End of actual data.\n";

	// shift data by k
	//int k = 1;
	//x = arrayOp.shiftArray(x, n, k);
	//for(int i=0; i<n; i++) cout << x[i] << "\n";
	//printf("End of %i shifted data.\n", k);
	
	// take fft of data
	in = reinterpret_cast<fftw_complex*>(x);
	fft = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(fft);
	for(int i=0; i<n; i++) {printf("freq %i: %f + %f*i \n", i, out[i][0], out[i][1]);}
	cout << "End of FFT of data.\n";



	// take absolute value of fftw_complex data type
	absIn = arrayOp.complexAbs(out, n);
	for(int i=0; i<n; i++) {cout << absIn[i] << "\n";}
	cout << "End of absolute value of FFT of data.\n";


	free(x);
	free(absIn);
	fftw_destroy_plan(fft);
	fftw_cleanup();

	//complex<double> *x;
	//fftw_complex *in, *out;
	//int n = (int) pow(2.0, 18.0);
	//x = (complex<double>*) malloc(sizeof(complex<double>) * n);

	//in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
	//in = reinterpret_cast<fftw_complex*>(x);
 //   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
	//fftw_plan p;
	//int i;

	//auto start = std::chrono::high_resolution_clock::now();
	///* prepare a cosine wave */
	//for (i = 0; i < n; i++) {
	//	x[i] = complex<double>(cos(3*2*m_pi*i/n), sin(2*2*m_pi*i/n)); //cos(3 * 2*m_pi*i/n);		
	//}
	//auto elapsed = std::chrono::high_resolution_clock::now() - start;
	//long long microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
	//std::cout << "data generation took: " << microseconds << " microseconds\n";
	//p = fftw_plan_dft_1d(n, in, out, fftw_forward, fftw_estimate);




	///* forward fourier transform, save the result in 'out' */
	//start = std::chrono::high_resolution_clock::now(	);
	//for(int i=0; i<150; i++)
	//	fftw_execute(p);
	//elapsed = std::chrono::high_resolution_clock::now() - start;
	//microseconds = std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count();
	//printf("elapsed time for %i point complex fft is : ", n);
	//std::cout <<microseconds << " microseconds\n";
	//std::cout <<"one fft is: " << microseconds/150 << " microseconds\n";

	//fftw_destroy_plan(p);
	//fftw_cleanup();
	return 0;


}

