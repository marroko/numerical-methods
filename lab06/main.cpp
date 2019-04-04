#include <iostream>
#include <cmath>

double deltaX = 0.1;

inline double f(double x) {

	return (x-1.2)*(x-2.3)*((x-3.3)*(x-3.3));
}

inline double derive(double x, double delta) {

	return ((f(x+delta) - f(x-delta)) / (2*delta));
}

inline double nextStep(double x0, double x1, double (*ptr)(double)) {

	 return (x1 - (x1 - x0) / (ptr(x1) - ptr(x0)) * ptr(x1));
}

inline double u(double x) {

	return (f(x) / derive(x, deltaX));
}

void solution(double x0, double x1, double (*ptr)(double))  {

	double x = nextStep(x0, x1, ptr);
	double tmp = x;
	double epsilon = std::abs(x - x1);
    int i = 1;

    while(std::abs(x - x0) > 1e-6) {

        std::cout << i++ << "&" << x << "&" << epsilon << "&" << f(x) << " \\\\ \\hline\n"; // udogodnienie do LaTeX'a 
    	x0 = x1;
    	x1 = x;
    	x = nextStep(x0, x1, ptr);
    	epsilon = std::abs(x - tmp);
    	tmp = x;
    }
}

int main() {

   	double (*ptr)(double) = f;

   	// do rysowania funkcji
   	double x = 0.9;

   	while(x < 3.700001) {

        std::cout << x << " " << f(x) << " " << derive(x, deltaX);
   		x += deltaX;
   	}

    // niemodyfikowana 
    std::cout << "x = 1.2\n";
    solution(0.9, 1.0, ptr);
    std::cout << "\nx = 2.3\n";
    solution(1.7, 1.75, ptr);
    std::cout << "\nx = 3.3\n";
    solution(3.7, 3.65, ptr);

    // modyfikowana
    ptr = u;
    std::cout << "\nModyfikowana x = 3.3, delta = 0.1\n";
    solution(3.7, 3.65, ptr);

    // obliczenia dla deltaX = 0.001
    std::cout << "\nModyfikowana x = 3.3, delta = 0.001\n";
    deltaX = 0.001;
    solution(3.7, 3.65, ptr);

    return 0;
}
