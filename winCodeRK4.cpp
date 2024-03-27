#include <opencv2/opencv.hpp>
#include <vector>
#include <mutex>
#include <thread>
#include <iostream>


void visualize(const std::vector<double>& xs, const std::vector<double>& ys) {
    int rows = xs.size();
    int cols = ys.size();

    // Создаем изображение
    cv::Mat img(rows, cols, CV_8UC1);

    // Заполняем изображение значениями
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Преобразуем значение в диапазоне 0-255
            img.at<uchar>(i, j) = cv::saturate_cast<uchar>(xs[i] * ys[j] * 255);
        }
    }

    // Показываем изображение
    cv::imshow("Brusselator Visualization", img);
    cv::waitKey(0);
}

// Параметры брюсселятора
double a = 1.0;
double b = 3.0;

// Функция для первого уравнения брюсселятора
double f(double x, double y) {
    return a - (b + 1) * x + x * x * y;
}

// Функция для второго уравнения брюсселятора
double g(double x, double y) {
    return b * x - x * x * y;
}

void euler(double* x, double* y, double dt) {
    double dx = dt * f(*x, *y);
    double dy = dt * g(*x, *y);
    *x += dx;
    *y += dy;
}

void rungeKutta(double &x, double &y, double dt, double (*f)(double, double), double (*g)(double, double)) {
    double k1, k2, k3, k4;
    double l1, l2, l3, l4;

    k1 = dt * f(x, y);
    l1 = dt * g(x, y);

    k2 = dt * f(x + dt / 2, y + l1 / 2);
    l2 = dt * g(x + dt / 2, y + l1 / 2);

    k3 = dt * f(x + dt / 2, y + l2 / 2);
    l3 = dt * g(x + dt / 2, y + l2 / 2);

    k4 = dt * f(x + dt, y + l3);
    l4 = dt * g(x + dt, y + l3);

    x = x + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
    y = y + (1.0 / 6.0) * (l1 + 2 * l2 + 2 * l3 + l4);
}

int main() {
    double x = 0.0, y = 0.0;
    double dt = 0.01;

    std::vector<double> xs, ys;

    for (int i = 0; i < 10000; i++) {
        rungeKutta(x, y, dt, f, g);
        xs.push_back(x);
        ys.push_back(y);
    }

    visualize(xs, ys);

    return 0;
}