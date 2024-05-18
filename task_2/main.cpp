#include <iostream>
#include <iomanip>
#include "math.h"
#include <cmath>

constexpr int TAYLOR_SERIES_TERMS = 10;

double cos_terms(double x, double n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        double term = pow(-1, i) * pow(x, 2 * i) / tgamma(2 * i + 1);
        sum += term;
    }
    return sum;
}

double double_double_cos(double * delta, double x, double n) {
    double cos_x = cos(x);
    double cos_x_approx = cos_terms(x, n);
    double delta_cos = cos_x - cos_x_approx;
    *delta = delta_cos;
    return cos_x;
}

int main() {
    double delta;

    double x = 1.5;
    double cos_x = double_double_cos(&delta, x, TAYLOR_SERIES_TERMS);
    double double_cos_x = cos_terms(x, TAYLOR_SERIES_TERMS);
    std::cout << "Cos:\t\t   " << std::setprecision(35) << cos_x
    << "\n" << "Double double cos: " << std::setprecision(35) << double_cos_x 
    << "\n" << "Delta:\t\t   " << std::setprecision(35)<< delta << std::endl;
    return 0;
}