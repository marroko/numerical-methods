#include <iostream>
#include <cmath>
#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/gaussj.c"

// nalezy zmieniac tylko parametr n oraz funkcje w 
// liniach 61, 63-64, 96 - w celu rozpatrzenia roznych przypadkow

constexpr int n = 14; 
constexpr float step = 0.01;
constexpr float xLeft = -5.0;
constexpr float xRight = 5.0;
constexpr float h = (xRight - xLeft) / (n-1);

inline float f1(float x) { return (1. / (1. + x*x)); }
inline float f1derive(float x) { return ((f1(x + step) - f1(x - step)) / (2*step)); }
inline float f2(float x) { return cos(2.*x); }
inline float f2derive(float x) { return ((f2(x + step) - f2(x - step)) / (2*step)); }

float s(float x, float *cc, float *xw) {

    float sum = 0.0;
    float cube = 1. / (h*h*h);

    for(int i = 0; i<=n+1; ++i) {

        cube = 1. / (h*h*h);

        if(x >= xw[i-2] && x < xw[i-1])
            cube *= pow(x - xw[i-2], 3);
        else if(x >= xw[i-1] && x < xw[i])
            cube *= pow(h, 3) + 3*h*h*(x - xw[i-1]) + 3*h*pow(x - xw[i-1], 2) - 3*pow(x - xw[i-1], 3);
        else if(x >= xw[i] && x < xw[i+1])
            cube *= pow(h, 3) + 3*h*h*(xw[i+1] - x) + 3*h*pow(xw[i+1] - x, 2) - 3*pow(xw[i+1] - x, 3);
        else if(x >= xw[i+1] && x < xw[i+2])
            cube *= pow(xw[i+2] - x, 3);
        else 
            cube = 0.0;

        sum += cc[i]*cube;
    }

    return sum;
}

int main() {

    float *c = vector(0, n+1);
    float **b = matrix(1, n, 1, 1);
    float **A = matrix(1, n, 1, n);
    float *xw = vector(-2, n+3);
    float *yw = vector(1, n);
    float x = xLeft;

    for(int i = -2; i<=n+3; ++i)
        xw[i] = xLeft + h*(i-1);

    for(int i = 1; i<=n; ++i)
        yw[i] = f2(xw[i]);
    
    float alfa = f2derive(xLeft);
    float beta = f2derive(xRight);

    for(int i = 1; i<=n; ++i) {

        for(int j = 1; j<=n; ++j) {

    		if(i == j)
    			A[i][j] = 4.0;
    		else if(i == j+1 || j == i+1)
    			A[i][j] = 1.0;
    		else
    			A[i][j] = 0.0;
    	}
    }

    for(int i = 1; i<=n; ++i)
    	b[i][1] = yw[i];

    b[1][1] = yw[1] + (h / 3. * alfa);
    b[n][1] = yw[n] - (h / 3. * beta);

    A[1][2] = A[n][n-1] = 2.0;

    gaussj(A, n, b, 1);

    for(int i = 1; i<=n; ++i)
    	c[i] = b[i][1];

    c[0] = c[2] - alfa*h/3;
    c[n+1] = c[n-1] + beta*h/3;

    while(x < xRight) {

        std::cout << x << " " << f2(x) << " " << s(x, c, xw) << std::endl;
        x += step;
    }

    free_vector(xw, -2, n+3);
    free_vector(yw, 1, n);
    free_vector(c, 0, n+1);
    free_matrix(b, 1, n, 1, 1);
    free_matrix(A, 1, n, 1, n);

    return 0;
}