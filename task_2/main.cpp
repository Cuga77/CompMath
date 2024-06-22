#include <iostream>
#include <cmath>
#include <cstdio>
#include <limits>

struct DoubleDouble {
    double hi;
    double lo;
    explicit DoubleDouble(double h = 0, double l = 0) : hi(h), lo(l) {}
};

DoubleDouble utilTwoSum(double a, double b) {
    double s = a + b;
    double t = s - a;
    double err = (a - (s - t)) + (b - t);
    return DoubleDouble(s, err);
}

DoubleDouble add(const DoubleDouble& a, const DoubleDouble& b) {
    DoubleDouble sum = utilTwoSum(a.hi, b.hi);
    double s = sum.hi;
    double e = sum.lo + a.lo + b.lo;
    sum = utilTwoSum(s, e);
    return sum;
}

DoubleDouble utilTwoProd(double a, double b) {
    double p = a * b;
    double err = std::fma(a, b, -p);
    return DoubleDouble(p, err);
}

DoubleDouble multiply(const DoubleDouble& a, const DoubleDouble& b) {
    DoubleDouble p1 = utilTwoProd(a.hi, b.hi);
    p1.lo += a.hi * b.lo + a.lo * b.hi;
    DoubleDouble prod = utilTwoSum(p1.hi, p1.lo);
    return prod;
}

DoubleDouble divide(const DoubleDouble& a, const DoubleDouble& b) {
    double q1 = a.hi / b.hi;
    DoubleDouble r = add(a, DoubleDouble(-q1 * b.hi, -q1 * b.lo));
    double q2 = r.hi / b.hi;
    DoubleDouble result = add(DoubleDouble(q1, 0), DoubleDouble(q2, 0));
    return result;
}

DoubleDouble negate(const DoubleDouble& a) {
    return DoubleDouble(-a.hi, -a.lo);
}

DoubleDouble sin(const DoubleDouble& x) {
    DoubleDouble term = x;
    DoubleDouble sum = term;
    DoubleDouble x_squared = multiply(x, x);
    int n = 1;

    while (true) {
        DoubleDouble denom1 = DoubleDouble(2 * n, 0.0);
        DoubleDouble denom2 = DoubleDouble(2 * n + 1, 0.0);
        DoubleDouble denom = multiply(denom1, denom2);

        term = multiply(term, negate(x_squared));
        term = divide(term, denom);

        DoubleDouble new_sum = add(sum, term);

        DoubleDouble relativeError = divide(term, new_sum);
        if (std::fabs(relativeError.hi + relativeError.lo) < std::numeric_limits<double>::epsilon()) {
            break;
        }

        sum = new_sum;
        ++n;
    }

    return sum;
}

int main() {
    double input;
    std::printf("Enter the value of x to compute sin(x): ");
    std::scanf("%lf", &input);

    DoubleDouble x(input, 0.0);
    DoubleDouble result = sin(x);
    
    std::printf("sin(%.40f) = %.40e + %.40e\nSummed result = %.40e\n", input, result.hi, result.lo, result.hi + result.lo);

    return 0;
}