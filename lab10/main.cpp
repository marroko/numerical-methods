#include <iostream>
#include <cmath>

//gradient wskazuje w ktora strone funkcja rosnie najbardziej

constexpr double h = 0.1;
constexpr double delta = 1e-4;
constexpr double epsilon = 1e-3;
constexpr int maxIter = 1000;

inline double f(double x, double y) { return ((5./2.) * pow(x*x - y, 2) + pow(1 - x, 2)); }

inline double gradX(double x, double y) { return ((f(x + delta, y) - f(x - delta, y)) / (2*delta)); }
inline double gradY(double x, double y) { return ((f(x, y + delta) - f(x, y - delta)) / (2*delta)); }

inline double norm(double *rNext, double *r) { 

	return (sqrt(pow(rNext[0] - r[0], 2) + pow(rNext[1] - r[1], 2))); 
}

void next(double &x, double &y) {

	double xtmp = x;
	x -= h*gradX(xtmp, y);
	y -= h*gradY(xtmp, y);
}

int main() {

    double r[2] = { -0.75, 1.75 }; // r[0] --- x, r[1] --- y
    double rNext[2] = { -0.75, 1.75 };
    double normCond = 1.;
    int counter = 0;

    while(normCond > epsilon && counter < maxIter) {

    	r[0] = rNext[0];
    	r[1] = rNext[1];

    	next(rNext[0], rNext[1]);
    	normCond = norm(rNext, r);

    	counter++;

    	std::cout << rNext[0] << " " << rNext[1] << '\n';
    }
    // std::cout << counter << " " << r[0] << " " << r[1] << '\n';
    // std::cout << counter << " " << rNext[0] << " " << rNext[1] << '\n';

    return 0;
}