#include "../ObjectiveNR.h"

int main() {

	constexpr int N = 4;
	int *indx = ivector(1, N);

    Matrix A(N, N, "A"), revA(N, N, "revA"), LU(N, N, "LU"), L(N, N, "L"), U(N, N, "U");
    Matrix *multiArevA;

    float *a, *b, *c, *dd, d, detA = 1.0;

	for(int i=1; i<=N; ++i) {

		for(int j=1; j<=N; ++j) {

            A(i,j) = 1.0 / (i+j);
            LU(i,j) = A(i,j);
		}
	}

    ludcmp(LU(), N, indx, &d);

	// obliczenie zawartosci macierzy L, U
	for(int i=1; i<=N; ++i) {

		for(int j=1; j<=N; ++j) {

			if(i == j) {

                L(i,j) = 1.0;
                U(i,j) = LU(i,j);
			}
			if(j > i) {

                L(i,j) = 0.0;
                U(i,j) = LU(i,j);
			}
			if(j < i) {

                L(i,j) = LU(i,j);
                U(i,j) = 0.0;
			}
		}
	}

	//obliczenie wyznacznika macierzy A
	for(int i=1; i<=N; ++i)
        detA *= U(i,i);


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

    lubksb(LU(), N, indx, a);
    lubksb(LU(), N, indx, b);
    lubksb(LU(), N, indx, c);
    lubksb(LU(), N, indx, dd);

	for(int j=1; j<=N; ++j)
        revA(j,1) = a[j];

	for(int j=1; j<=N; ++j)
        revA(j,2) = b[j];

	for(int j=1; j<=N; ++j)
        revA(j,3) = c[j];

	for(int j=1; j<=N; ++j)
        revA(j,4) = dd[j];
	
	//obliczanie iloczynu A oraz revA

    multiArevA = A*revA;
    multiArevA->setName("multiArevA");

	//obliczanie wskaznikow uwarunkowania macierzy

    float ptrA = A.ptrMatrix();
    float ptrRevA = revA.ptrMatrix();
	float kappa = ptrA * ptrRevA;

	//wypisanie wynikow

    std::cout << A << LU << L << U;
	std::cout << "detA: " << detA << "\n\n";
    std::cout << revA << *multiArevA;

	std::cout << "Wskazniki uwarunkowan: " << "\n"
			  << "ptrA:    " << ptrA << "\n"
			  << "ptrRevA: " << ptrRevA << "\n"
			  << "kappa:   " << kappa << "\n";

    delete multiArevA;

	free_ivector(indx, 1, N); 
	free_vector(a, 1, N); 
	free_vector(b, 1, N); 
	free_vector(c, 1, N); 
	free_vector(dd, 1, N); 
	
	return 0;
}
