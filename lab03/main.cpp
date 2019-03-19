#include "../ObjectiveNR.h"

#define max(X,Y) ((X)>(Y)? (X):(Y))
#define min(X,Y) ((X)<(Y)? (X):(Y))
#define abs(X) ((X)>0? (X):-(X))

constexpr int N = 400, m = 5;

FVector wstega(const Matrix &A, const FVector &x) {

    int jmin, jmax;
    FVector y(N, "y");

    for(int i=1; i<=N; i++) {

        jmin = max(1, 2-m);
        jmax = min(i+1+m, N);
        y[i] = 0.0;
        for(int j=jmin; j<=jmax; j++)
            y[i] += A(i,j) * x[j];
    }

    return y;
}

DVector wstega(const DMatrix &A, const DVector &x) {

    int jmin, jmax;
    DVector y(N, "y");

    for(int i=1; i<=N; i++) {

        jmin = max(1, 2-m);
        jmax = min(i+1+m, N);
        y[i] = 0.0;
        for(int j=jmin; j<=jmax; j++)
            y[i] += A(i,j) * x[j];
    }

    return y;
}

int main() {

    // zmienne do pomiaru czasu wykonania obliczen
    auto t1 = std::chrono::steady_clock::now();
    auto t2 = std::chrono::steady_clock::now();
    auto diff = std::chrono::duration<double, std::milli>(t2 - t1).count();

    Matrix A(N, N, "A"), x(N, 1, "x");
    DMatrix AA(N, N, "AA");
    FVector b(N, "b"), r(N, "r");
    DVector bb(N, "bb"), xx(N, "xx"), rr(N, "rr");
    float alfa;
    int k = 1;

	// r k - numer iteracji
	// r - wektor poprawek

    std::cout << "iteracja  normaR  alfa  normaX\n";
    float rNorm;
    double rrNorm;

    int i=0;

    // pomiar czasu dla 10 prob metoda Gaussa-Jordana

    while(i++<10) {

        for(int i=1; i<=N; ++i) {

            b[i] = i;
            bb[i] = i;
            x(i,1) = 0.0; // najpierw elementy wektora x rowne 0
            xx[i] = 0.0; // najpierw elementy wektora xx rowne 0

            for(int j=1; j<=N; ++j) {

                if(abs(i-j) <= m) {
                    A(i,j) = 1.0 / (1 + abs((i-j)));
                    AA(i,j) = 1.0 / (1 + abs((i-j)));
                }
                else {
                    A(i,j) = 0.0;
                    AA(i,j) = 0.0;
                }
            }
        }
        auto t1 = std::chrono::steady_clock::now();
        gaussj(A(), N, x(), 1);

        auto t2 = std::chrono::steady_clock::now();
        auto diff = std::chrono::duration<double, std::milli>(t2 - t1).count();

        std::cout << diff << '\n';
    }

    std::cout << "\n";
    k = 1;

    // obliczanie dla podwojnej precyzji
   do {

        rr = bb - wstega(AA, xx);

        alfa = (rr*rr) / (rr*wstega(AA, rr));
        xx = xx + rr*alfa;

        rrNorm = rr.normVector();

        std::cout << k << " " << rrNorm << " " << alfa << " " << xx.normVector() << '\n';

        k++;

   } while(rrNorm > 1e-6);
	
	return 0;
}
