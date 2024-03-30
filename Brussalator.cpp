#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <matplotlibcpp.h> // Добавьте эту библиотеку для визуализации
#include <numpy/arrayobject.h> // Добавьте эту библиотеку для работы с массивами

namespace plt = matplotlibcpp;

// Параметры брюсселятора
const double a = 1.0;
const double b = 3.0;

// Функция для первого уравнения брюсселятора
double f(double x, double y) {
    return a - (b + 1) * x + x * x * y;
}

// Функция для второго уравнения брюсселятора
double g(double x, double y) {
    return b * x - x * x * y;
}

void rk4(double &x, double &y, double dt, double (*f)(double, double), double (*g)(double, double)) {
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
        rk4(x, y, dt, f, g);
        xs.push_back(x);
        ys.push_back(y);
    }

    // Преобразование векторов в numpy массивы для визуализации
    PyObject* x_array = PyArray_SimpleNewFromData(1, &xs.size(), NPY_DOUBLE, (void*)xs.data());
    PyObject* y_array = PyArray_SimpleNewFromData(1, &ys.size(), NPY_DOUBLE, (void*)ys.data());

    // Визуализация
    plt::plot(x_array, y_array);
    plt::show();

    // Освобождение памяти
    Py_DECREF(x_array);
    Py_DECREF(y_array);

    return 0;
}