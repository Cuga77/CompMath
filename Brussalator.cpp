#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>


void rk4(double x, double y, double alpha, double beta, double delta, double gamma, double h, std::vector<double>& x_values, std::vector<double>& y_values) {
    double k1_x = h * (alpha * x - beta * x * y);
    double k1_y = h * (delta * x * y - gamma * y);

    double k2_x = h * (alpha * (x + 0.5 * k1_x) - beta * (x + 0.5 * k1_x) * (y + 0.5 * k1_y));
    double k2_y = h * (delta * (x + 0.5 * k1_x) * (y + 0.5 * k1_y) - gamma * (y + 0.5 * k1_y));

    double k3_x = h * (alpha * (x + 0.5 * k2_x) - beta * (x + 0.5 * k2_x) * (y + 0.5 * k2_y));
    double k3_y = h * (delta * (x + 0.5 * k2_x) * (y + 0.5 * k2_y) - gamma * (y + 0.5 * k2_y));

    double k4_x = h * (alpha * (x + k3_x) - beta * (x + k3_x) * (y + k3_y));
    double k4_y = h * (delta * (x + k3_x) * (y + k3_y) - gamma * (y + k3_y));

    x += (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6;
    y += (k1_y + 2 * k2_y + 2 * k3_y + k4_y) / 6;

    x_values.push_back(x);
    y_values.push_back(y);
}

int main() {
    double alpha = 1.0, beta = 2.0, delta = 1.0, gamma = 0.3;
    double x = 1.0, y = 0.0;
    double h = 0.01;
    int steps = 1000;

    std::vector<double> x_values, y_values;

    for (int i = 0; i < steps; ++i) {
        rk4(x, y, alpha, beta, delta, gamma, h, x_values, y_values);
    }

    std::ofstream file("brusselator_data.txt");
    for (size_t i = 0; i < x_values.size(); ++i) {
        file << x_values[i] << " " << y_values[i] << std::endl;
    }
    file.close();

    std::cout << "Результаты сохранены в файл brusselator_data.txt" << std::endl;

    return 0;
}