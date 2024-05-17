#include <iostream>
#include "math.h"
#include <cmath>

double sum_two(double * delta, double a, double b) {
    double s = a + b;
    double bs = s - b;
    double as = s - a;
    *delta = (b - bs) + (a - as); 
    return s;
}

double two_mult(double * delta, double a, double b) {
    double p = a * b;
    *delta = fma(a, b, -p);
    return p;
}

int main() {
    double a = 2e-16;
    double b = 1e-16;
    double delta, delta2;
    double s = sum_two(&delta, a, b);
    std::cout << "Sum: " << s << " Delta: " << delta << std::endl;
    double p = two_mult(&delta2, a, b);
    std::cout << "Mult: " << p << " Delta: " << delta2 << std::endl;
    return 0;
}