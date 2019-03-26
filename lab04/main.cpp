#include "../ObjectiveNR.h"

int main() {

    constexpr int N = 5;

    Matrix A(N, N, "A"), P(N, N, "przeksztalcenia P"), Y(N, N, "Y"), X(N, N, "X");
    FVector d(N, "d"), e(N, "e"), tmp(N, "tmp");
    float sum1 = 0.0;
    float sum2 = 0.0;
    float beta;

    for(int i = 1; i<=N; ++i) {

        for(int j=1; j<=N; ++j) {
            A(i,j) = sqrt(i + j);

            if(i == j)
                Y(i,j) = 1.0;
            else 
                Y(i,j) = 0.0;
        }
    }

    /* dla wygodnego odczytu elementow macierzy 
    nalezy zmienic 353 linijke w pliku 
    "ObjectiveNR.h" na:
    output << std::setw(12) << m(i,j) << " "; */

    P = A;

    // redukcji macierzy P (wlasciwie A) do postaci trójdiagonalnej

    tred2(P(), N, d(), e());

    std::cout << A << P;
    
    tqli(d(), e(), N, Y());

    // kolumny Y przechowuja kolejne wektory wlasne T
    std::cout << Y;

    // wartosci wlasne macierzy A
    std::cout << d;

    X = P*Y;

    // wektory wlasne xk macierzy A dla odpowiednich wartosci wlasnych
    std::cout << X;

    // obliczanie kolejnych wspolczynnikow beta
    for(int i=1; i <= N; ++i)
        tmp[i] = 0.0;

    for(int i = 1; i <= N; i++) {

        for(int k = 1; k <= N; k++) {

            for(int l = 1; l <= N; ++l)
                tmp[k] += A(k,l) * X(l,i);
        }

        for(int j = 1; j <= N; j++) {

            sum1 += X(j,i)*tmp[j];
            sum2 += X(j,i)*X(j,i);
        }

        beta = sum1 / sum2;
        sum1 = sum2 = 0.0;
        for(int a=1; a <= N; ++a)
            tmp[a] = 0.0;

        std::cout << "beta " << i << " = " << beta << std::endl;
    }

    // czy betak jest rowne obliczonym wczesniej wartosciom wlasnym?
	
	return 0;
}
