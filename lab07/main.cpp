#include <iostream>
#include <cmath>

constexpr double leftX = -5.0;
constexpr double rightX = 5.0;
constexpr int N = 20;

inline double f(double x) { return std::exp(-(x*x)); }

inline double czeb(double m) { return 0.5*(rightX - leftX) * std::cos(M_PI * (2*m + 1) / (2*N + 2)) + (leftX + rightX); }

void interp(double x, const int N, double xNodes[], double yNodes[]) {

    double sum = 0.0;
    double multiply;

    for(int j=0; j<N+1; ++j) {

        multiply = 1.0;

        for(int k=0; k<N+1; ++k) {

            if(k != j)
                multiply *= (x - xNodes[k]) / (xNodes[j] - xNodes[k]);
        }

        multiply *= yNodes[j];
        sum += multiply;
    }

    std::cout << sum << "\n";
}

int main() {

    constexpr double step = 0.01; // co ile obliczamy wartosc wielomianu
    double currentPos = leftX;    // zaczynamy od -5.0

    constexpr double range = rightX - leftX;
    constexpr double littleRange = range / N;

    double nodesX[N+1];
    double nodesY[N+1];

    for(int i=0; i<N+1; ++i) {

        nodesX[i] = czeb(i); // w pierwszym zadaniu leftX + littleRange*i;
        nodesY[i] = f(nodesX[i]);
    }

    while(currentPos < rightX) {

        std::cout << currentPos << " ";
        interp(currentPos, N, nodesX, nodesY);
        currentPos += step;
    }

    return 0;
}