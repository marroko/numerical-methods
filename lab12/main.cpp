#include <iostream>
#include <cmath>

#include</opt/NR/numerical_recipes.c/nrutil.c>
#include</opt/NR/numerical_recipes.c/nrutil.h>
#include</opt/NR/numerical_recipes.c/gauleg.c>
#include</opt/NR/numerical_recipes.c/gammln.c>
#include</opt/NR/numerical_recipes.c/gaulag.c>
#include</opt/NR/numerical_recipes.c/gauher.c>

constexpr float left = 0.;
constexpr float right = 2.;
constexpr float alpha = 0.;

inline float c1(float x) { return x / (4*x*x + 1); }
inline float c2(float x, int k) { return pow(x, k); }
inline float c3(float x, int toPow) { return pow(sin(x), toPow); }
inline float c1a(float x, float a, float c) { return log(a*a*x*x + c*c) / (2*a*a); }
inline float c1Int() { return c1a(right, 2, 1) - c1a(left, 2, 1); }

inline int factorial(int n) { return (n == 1) ? 1 : (n * factorial(n-1)); }

void computeGauLag(int N, int k, const char *file) {

	FILE *f2 = fopen(file, "w");

	while (N < 21) {

		double sum = 0.;

		float *x = vector(1, N);
		float *w = vector(1, N);

    	gaulag(x, w, N, alpha);

    	for(int i = 1; i<=N; ++i)
    		sum += c2(x[i], k)*w[i];

    	fprintf(f2, "%d %f\n", N, sum - factorial(k));

    	//std::cout << sum << " " << factorial(k) << '\n';

    	free_vector(x, 1, N);
    	free_vector(w, 1, N);

    	N++;
    }
    fclose(f2);
}

int main() {

	int N = 2, k = 5;
	float sum = 0.;

	FILE *f1 = fopen("gauleg.dat", "w");
	FILE *f3 = fopen("gauher.dat", "w");

	/// 1 ///
	while (N < 21) {

		sum = 0.;

		float *x = vector(1, N);
		float *w = vector(1, N);

    	gauleg(left, right, x, w, N);

    	for(int i = 1; i<=N; ++i)
    		sum += c1(x[i])*w[i];

    	fprintf(f1, "%d %f\n", N, sum - c1Int());

    	free_vector(x, 1, N);
    	free_vector(w, 1, N);

    	N++;
    }

    /// 2 ///
    N = 2;
    computeGauLag(N, k, "gaulag5.dat");

    N = 2;
    k = 10;
    computeGauLag(N, k, "gaulag10.dat");

    /// 3 ///
    N = 2;
    float sum2 = 0.;
    float cdok = 0.1919832644;

    while (N < 16) {

		sum = 0.;
		sum2 = 0.;

		float *x = vector(1, N);
		float *w = vector(1, N);

    	gauher(x, w, N);

    	for(int i = 1; i<=N; ++i)
    		sum += c3(x[i], 2)*w[i];

    	for(int i = 1; i<=N; ++i)
    		sum2 += sum*c3(x[i], 4)*w[i];

    	fprintf(f3, "%d %f\n", N, sum2 - cdok);

    	free_vector(x, 1, N);
    	free_vector(w, 1, N);

    	N++;
    }

    fclose(f1);
    fclose(f3);

    return 0;
}