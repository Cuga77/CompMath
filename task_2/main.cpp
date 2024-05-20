#include <iostream>
#include <cmath>
#include <iomanip>

typedef long double ld;

constexpr ld eps = 1e-254;

ld cos_dd(ld x, ld eps) {
    // Уменьшение x до значения в пределах от -π до π
    x = fmod(x, 2 * M_PI);
    if (x < -M_PI) {
        x += 2 * M_PI;
    } else if (x > M_PI) {
        x -= 2 * M_PI;
    }

    ld term = 1;
    ld sum = term;
    int i = 0;

    while (std::abs(term) > eps) {
        term *= -x * x / ((2 * i + 1) * (2 * i + 2)); // -x^2 / (2n + 1)(2n + 2) 
        sum += term;
        i++;
    }

    return sum;
}

int main() {
    ld x, delta;
    std::cin >> x;
    std::cout << std::setprecision(50) << std::fixed;
    std::cout << "double_double cos x = "  << cos_dd(x, eps) << std::endl;
    std::cout << "\t      cos x = " << cos(x) << std::endl;
    delta = std::abs(cos_dd(x, eps) - cos(x));
    std::cout << "\t      delta = " << delta << std::endl;

    return 0;
}