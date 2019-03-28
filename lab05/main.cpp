#include "../ObjectiveNR.h"

int main() {

    constexpr int N = 7, IT_MAX = 12;
    Matrix A(N, N, "A"), W(N, N, "W"), tmp(N, N, "tmp"), X(N, N, "X"), D(N, N, "D");
    FVector x(N, "x"), y(N, "y");
    float lambda;

    for(int i = 1; i<=N; ++i) {

        for(int j=1; j<=N; ++j)
            A(i,j) = 1.0 / sqrt(2.0 + abs(i - j));
    }

    W = A;

    for(int k = 1; k<=N; ++k) {

        for(int j = 1; j<=N; ++j)
            x[j] = 1.0; 
        
        for(int i = 1; i<=IT_MAX; ++i) {

            y = W*x;
            lambda = y*x / (x*x);
            //std::cout << i << " " << lambda << "\n";

            x = y * (1.0 / y.normVector());
        }
        //std::cout << k << " " << x;

        for(int j = 1; j<=N; ++j)
            X(j,k) = x[j];

        // iloczyn tensorowy
        for(int j = 1; j<=N; ++j) {
            
            for(int l = 1; l<=N; ++l)
                tmp(j,l) = x[j]*x[l];
        }
        W = W - tmp*lambda;  
    }

    std::cout << X;

    D = X.transpose() * A * X;

    std::cout << D;

    return 0;
}