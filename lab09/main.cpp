#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <gsl/gsl_linalg.h>

constexpr double sigma = 4.;
constexpr double x0 = 2.;
constexpr double left = -3*sigma + x0;
constexpr double right = 3*sigma + x0;
constexpr double alfa = 0.5;
constexpr int N = 101;
constexpr int m = 4;
const double a0 = -x0*x0 / (2.*sigma*sigma);
const double a1 = x0 / (sigma*sigma);
const double a2 = -1. / (2.*sigma*sigma);

const double step = (right - left) / (N - 1);

inline double g(double x) {

	const double U = rand() / (RAND_MAX + 1.0);
	return (exp(a0 + a1*x + a2*x*x)*(1. + alfa*(U - 0.5)));
}

inline double f(double x) {

	return log(g(x));
}

inline double GG(double x, gsl_vector *r) {

	return (exp(gsl_vector_get(r, 0) + gsl_vector_get(r, 1)*x + gsl_vector_get(r, 2)*x*x + gsl_vector_get(r, 3)*x*x*x));
}

int main() {

    srand(time(NULL));
    double x[N];
    double fx[N];
    constexpr double approxStep = 0.1;
    double xtmp = left;


    gsl_vector *r = gsl_vector_calloc(m);
    gsl_matrix *G = gsl_matrix_calloc(m, m);

    for(int i = 0; i<N; ++i) {

    	x[i] = left + i*step;
    	fx[i] = f(x[i]);
    }

    for(int k = 0; k<m; ++k) {

    	double val = 0;

    	for(int j = 0; j<N; ++j)
    		val += fx[j]*pow(x[j], k);

    	gsl_vector_set(r, k, val);
    }

    for(int i = 0; i<m; ++i) {

    	for(int k = 0; k<m; ++k) {

    		double val = 0;
    		for(int j = 0; j<N; ++j)
    			val += pow(x[j], i+k);

    		gsl_matrix_set(G, i, k, val);
    	}
    }

    gsl_linalg_HH_svx(G, r);

   // std::cout << "wspolczynniki b\n";
   // for(int i = 0; i<m; ++i)
   // 	std::cout << gsl_vector_get(r, i) << std::endl;

   // std::cout << "wspolczynniki a\n";
   // std::cout << a0 << " " << a1 << " " << a2 << std::endl;
    
    // aproksymacja
    while(xtmp < right) {

    	std::cout << xtmp << " " << GG(xtmp, r) << '\n';
    	xtmp += approxStep;
    }
	
  //  for(int i=0; i<N; ++i)
//	    std::cout << x[i] << " " << exp(fx[i]) << '\n';

    gsl_matrix_free(G);
    gsl_vector_free(r);

    return 0;
}
