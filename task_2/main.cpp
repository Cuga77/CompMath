#include <iostream>
#include <iomanip>
#include "math.h"
#include <cmath>

constexpr int TAYLOR_SERIES_TERMS = 11;

double cos_terms(double x, double n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double term = pow(-1, i) * pow(x, 2 * i) / tgamma(2 * i + 1);
        sum += term;
    }
    return sum;
}

long double double_double_cos(long double * delta, double x, double n) {
    long double cos_x = cosl(x);
    double cos_x_approx = cos_terms(x, n);
    *delta = cos_x - static_cast<long double>(cos_x_approx);
    return cos_x;
}

int main() {
    long double delta;
    double x = 1.5;
    long double cos_x = double_double_cos(&delta, x, TAYLOR_SERIES_TERMS);
    double double_cos_x = cos_terms(x, TAYLOR_SERIES_TERMS);
    std::cout << "Cos:\t\t   " << std::setprecision(50) << cos_x
    << "\n" << "Double double cos: " << std::setprecision(50) << double_cos_x 
    << "\n" << "Delta:\t\t   " << std::setprecision(50)<< delta << std::endl;
    return 0;
}