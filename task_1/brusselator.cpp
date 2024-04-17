#include <iostream>
#include <cstring>
#include <cmath>
#include <fstream>

#include <SFML/Graphics.hpp>

// compile: g++ -o brusselator brusselator.cpp && .\a.exe || ./a.out

// Размеры окна
const int WINDOW_WIDTH = 200;
const int WINDOW_HEIGHT = 200;

// Параметры брюсселятора
double a = 2.9;
double b = 1.3;
double b1 = b + 1;
double h = 0.001;
int N = 1000;

//Инициализация начальных значений
double **x, **y, **vx, **vy;

// Функции для вычисления производных x и y
//b1 для ускорения работы
double X(double x, double y) {
    return a - b1 * x + x * x * y;
}

double Y(double x, double y) {
    return b * x - x * x * y;
}

void rk4(double** x, double** y, double** vx, double** vy, double dt) {
    double k1x, k2x, k3x, k4x;
    double k1y, k2y, k3y, k4y;
    double k1vx, k2vx, k3vx, k4vx;
    double k1vy, k2vy, k3vy, k4vy;

    for (int i = 0; i < WINDOW_WIDTH; i++) {
        for (int j = 0; j < WINDOW_HEIGHT; j++) {
            k1x = vx[i][j] * dt;
            k1y = vy[i][j] * dt;
            k1vx = X(x[i][j], y[i][j]) * dt;
            k1vy = Y(x[i][j], y[i][j]) * dt;

            k2x = (vx[i][j] + k1vx / 2) * dt;
            k2y = (vy[i][j] + k1vy / 2) * dt;
            k2vx = X(x[i][j] + k1x / 2, y[i][j] + k1y / 2) * dt;
            k2vy = Y(x[i][j] + k1x / 2, y[i][j] + k1y / 2) * dt;
            
            k3x = (vx[i][j] + k2vx / 2) * dt;
            k3y = (vy[i][j] + k2vy / 2) * dt;
            k3vx = X(x[i][j] + k2x / 2, y[i][j] + k2y / 2) * dt;
            k3vy = Y(x[i][j] + k2x / 2, y[i][j] + k2y / 2) * dt;

            k4x = (vx[i][j] + k3vx) * dt;
            k4y = (vy[i][j] + k3vy) * dt;
            k4vx = X(x[i][j] + k3x, y[i][j] + k3y) * dt;
            k4vy = Y(x[i][j] + k3x, y[i][j] + k3y) * dt;

            x[i][j] += (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
            y[i][j] += (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
            vx[i][j] += (k1vx + 2 * k2vx + 2 * k3vx + k4vx) / 6;
            vy[i][j] += (k1vy + 2 * k2vy + 2 * k3vy + k4vy) / 6;
        }
    }
}

int main() {

    //Параметры диффузии
    double D = 1.1;
    double dt = 0.001;
    double dx = 0.1;
    double dy = 0.1;     

    // Инициализация массивов для хранения значений x, y, vx, vy
    x = new double*[WINDOW_WIDTH];
    y = new double*[WINDOW_WIDTH];
    vx = new double*[WINDOW_WIDTH];
    vy = new double*[WINDOW_WIDTH];
    for (int i = 0; i < WINDOW_WIDTH; i++) {
        x[i] = new double[WINDOW_HEIGHT];
        y[i] = new double[WINDOW_HEIGHT];
        vx[i] = new double[WINDOW_HEIGHT];
        vy[i] = new double[WINDOW_HEIGHT];
        for (int j = 0; j < WINDOW_HEIGHT; j++) {
            x[i][j] = (double)rand() / RAND_MAX; // случайное число от 0 до 1
            y[i][j] = (double)rand() / RAND_MAX;
            vx[i][j] = (double)rand() / RAND_MAX;
            vy[i][j] = (double)rand() / RAND_MAX;
        }
    }

    std::ofstream file("output.txt");
    
    // Обновление системы в цикле
    for (int i = 0; i < WINDOW_WIDTH; i++) {
        for (int j = 0; j < WINDOW_HEIGHT; j++) {
            // Обновление значений с помощью метода Рунге-Кутты
            rk4(x, y, vx, vy, dt);
            // Запись текущих значений в файл
            file << x[i][j] << " " << y[i][j] << " " << vx[i][j] << " " << vy[i][j] << "\n";
        }
    }

    for (int i = 0; i < WINDOW_WIDTH; i++) {
        delete[] x[i];
        delete[] y[i];
        delete[] vx[i];
        delete[] vy[i];
    }
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;

    return 0;
}