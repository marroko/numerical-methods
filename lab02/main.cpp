#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/ludcmp.c"
#include "/opt/NR/numerical_recipes.c/lubksb.c"

void printMatrix(float **A, int size, std::string name) {

	std::cout << "\t\t\t" << name << "\n";

	for(int i=1; i<=size; ++i) {

		for(int j=1; j<=size; ++j)
			std::cout << A[i][j] << " ";
		std::cout << "\n";
	}
	std::cout << "\n";
}

float ptrMatrix(float **A, int size) {

	float max = fabs(A[1][1]);

	for(int i=1; i<=size; ++i) {

		for(int j=1; j<=size; ++j) {

			if(max < A[i][j])
				max = A[i][j];
		}
	}
	return max;
}

int main() {

	constexpr int N = 4;
	int *indx = ivector(1, N);

	float **A, **revA, **LU, **L, **U, **multiArevA, d, detA = 1.0;
	float *a, *b, *c, *dd;

	A = matrix(1, N, 1, N);
	revA = matrix(1, N, 1, N);
	multiArevA = matrix(1, N, 1, N);
	LU = matrix(1, N, 1, N);
	L = matrix(1, N, 1, N);
	U = matrix(1, N, 1, N);

	for(int i=1; i<=N; ++i) {

		for(int j=1; j<=N; ++j) {

			A[i][j] = 1.0 / (i+j);
			LU[i][j] = A[i][j];
		}
	}

	ludcmp(LU, N, indx, &d);

	// obliczenie zawartosci macierzy L, U
	for(int i=1; i<=N; ++i) {

		for(int j=1; j<=N; ++j) {

			if(i == j) {

				L[i][j] = 1.0;
				U[i][j] = LU[i][j];
			}
			if(j > i) {

				L[i][j] = 0.0;
				U[i][j] = LU[i][j];
			}
			if(j < i) {

				L[i][j] = LU[i][j];
				U[i][j] = 0.0;
			}
		}
	}

	//obliczenie wyznacznika macierzy A
	for(int i=1; i<=N; ++i)
		detA *= U[i][i];


	//wyznaczanie macierzy odwrotnej 
	a = vector(1, N);
	b = vector(1, N);
	c = vector(1, N);
	dd = vector(1, N);

	for(int i=1; i<=N; ++i) {

		a[i] = 0.0;
		b[i] = 0.0;
		c[i] = 0.0;
		dd[i] = 0.0;
	}

	a[1] = 1.0;
	b[2] = 1.0;
	c[3] = 1.0;
	dd[4] = 1.0;

	lubksb(LU, N, indx, a);
	lubksb(LU, N, indx, b);
	lubksb(LU, N, indx, c);
	lubksb(LU, N, indx, dd);

	for(int j=1; j<=N; ++j)
		revA[j][1] = a[j];

	for(int j=1; j<=N; ++j)
		revA[j][2] = b[j];

	for(int j=1; j<=N; ++j)
		revA[j][3] = c[j];

	for(int j=1; j<=N; ++j)
		revA[j][4] = dd[j];
	
	//obliczanie iloczynu A oraz revA

	for(int i=1; i<=N; ++i) {

		for(int j=1; j<=N; ++j) {

			multiArevA[i][j] = 0.0;
			for(int k=1; k<=N; ++k)
				multiArevA[i][j] += A[i][k] * revA[k][j];
		}
	}

	//obliczanie wskaznikow uwarunkowania macierzy

	float ptrA = ptrMatrix(A, N);
	float ptrRevA = ptrMatrix(revA, N);
	float kappa = ptrA * ptrRevA;

	//wypisanie wynikow

	printMatrix(A, N, "Macierz A");
	printMatrix(LU, N, "Macierz LU");
	printMatrix(L, N, "Macierz L");
	printMatrix(U, N, "Macierz U");
	std::cout << "detA: " << detA << "\n\n";
	printMatrix(revA, N, "Macierz revA");
	printMatrix(multiArevA, N, "Macierz multiArevA");

	std::cout << "Wskazniki uwarunkowan: " << "\n"
			  << "ptrA:    " << ptrA << "\n"
			  << "ptrRevA: " << ptrRevA << "\n"
			  << "kappa:   " << kappa << "\n";


	free_matrix(A, 1, N, 1, N);
	free_matrix(revA, 1, N, 1, N);
	free_matrix(multiArevA, 1, N, 1, N);
	free_matrix(LU, 1, N, 1, N);
	free_matrix(U, 1, N, 1, N);
	free_matrix(L, 1, N, 1, N);
	free_ivector(indx, 1, N); 
	free_vector(a, 1, N); 
	free_vector(b, 1, N); 
	free_vector(c, 1, N); 
	free_vector(dd, 1, N); 
	
	return 0;
}
