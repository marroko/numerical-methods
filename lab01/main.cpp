#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/gaussj.c"

int main() {

	constexpr int N = 400; // rozmiar macierzy A: NxN
	constexpr float omega = 1.0;
	constexpr float vPoczatkowa = 0.0;
	constexpr float h = 0.1;

	float **A,**b;
//	Alokacja macierzy
	A = matrix(1, N, 1, N);
	b = matrix(1, N, 1, 1);
	
// 	Wypelnienie macierzy A i wektora b
	for(int i=1; i<=N; ++i) {

		for(int j=1; j<=N; ++j) {

			if(i == j)
				A[i][j] = 1.0;
			else if(i == j+2 && i>2)
				A[i][j] = 1.0;
			else if(j == i-1 && j>1)
				A[i][j] = omega*omega*h*h - 2;
		}
	}

	A[2][1] = -1.0;
	b[1][1] = 1.0;
	b[2][1] = vPoczatkowa * h;

	for(int i=3; i<=N; ++i)
		b[i][1] = 0.0;
	
//	Rozwiazanie ukladu rownan Ax=b - wywolanie procedury:
	gaussj(A,N,b,1);

//	Wypisanie rozwiazania, ktore procedura gaussj(A,N,b,1); zapisala w wektorze b.
	for(int i=1; i<=N; ++i) {
		std::cout << h*(i-1) << " "
	 			  << b[i][1] 
	 			  << "\n"; 
	}
    
//	Zwolnienie pamieci
	free_matrix(A, 1, N, 1, N);
	free_matrix(b, 1, N, 1, 1);    
	
	return 0;
}

