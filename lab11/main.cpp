#include <iostream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

using namespace std;

constexpr double T = 1.0;
constexpr double tmax = 3*T;
constexpr double sigma = T/20;
const double omega = 2*M_PI / T;

inline double f0(double tt) {

	return sin(omega*tt) + sin(2*omega*tt) + sin(3*omega*tt);
} 

inline double delta() {

	return (rand() / (RAND_MAX + 1.0)) - (1./2.);
}

inline double gwage(double tt) {

	return (1. / (sigma * sqrt(2*M_PI))) * exp(-(tt*tt) / (2*sigma*sigma));
}

int main() {

	int k = 8, counter = 0; // dla innych przypadkow wystarczy zmienic parametr k
	const int N = pow(2, k);
	const double dt = tmax/N;

    double f[2*N], g1[2*N], g2[2*N];
    double t = 0.0;

    for(int i = 0; i<2*N; ++i) {

    	if(i % 2) { // nieparzyste

    		f[i] = 0.0;
    		g1[i] = g2[i] = 0.0;
    	}
    	else {

    		f[i] = f0(t) + delta();
    		g1[i] = g2[i] = gwage(t);
    	}

    	counter++;
    	if(!(counter % 2)) {

    		counter = 0;
    		t += dt;
    	}
    }

    FILE *fdelta = fopen("fdelta.dat", "w"), *fff = fopen("f.dat", "w");

    for(int i = 0; i<N; ++i)
    	fprintf(fdelta, "%f %f\n", i*dt, f[2*i]);

    for(int i = 0; i<N; ++i)
    	fprintf(fff, "%f %f\n", i*dt, f0(i*dt));

    fclose(fdelta);
    fclose(fff);

    gsl_fft_complex_radix2_forward(f, 1, N);
    gsl_fft_complex_radix2_forward(g1, 1, N);
    gsl_fft_complex_radix2_backward(g2, 1, N);

    double a1, a2, b1, b2;

    for(int i = 0; i<N; ++i) {

    	a1 = f[2*i];
    	b1 = f[2*i+1];
    	a2 = g1[2*i] + g2[2*i];
    	b2 = g1[2*i+1] + g2[2*i+1];
    	f[2*i] = a1*a2 - b1*b2;
    	f[2*i+1] = a1*b2 + b1*a2;
    }

    gsl_fft_complex_radix2_backward(f, 1, N);

    double max = 0.0;
    for(int i = 0; i<N; ++i) {

    	if(f[2*i] > max)
    		max = f[2*i];
    }

    FILE *output = fopen("output.dat", "w");

    for(int i = 0; i<N; ++i)
    	fprintf(output, "%f %f\n", i*dt,(f[2*i]) / (max / 2.5));

    fclose(output);

    return 0;
}