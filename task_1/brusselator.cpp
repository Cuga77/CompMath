#include <iostream>
#include <cstring>
#include <cmath>
#include <fstream>

// compile: g++ -o brusselator brusselator.cpp 

// Размеры окна
const int WINDOW_WIDTH = 100;
const int WINDOW_HEIGHT = 100;

// Параметры брюсселятора
double a = 1.3;
double b = 3.5;
double h = 0.1;
double x00 = 0.1, y00 = 0.1, vx = 4.2, vy = 6.2;
int N = 1000;


double* xr, *xv, *yr, *yv;

double X(double x, double y) {
    return a - (b + 1) * x + x * x * y;
}

double Y(double x, double y) {
    return b * x - x * x * y;
}

void rk4() {
    double k[4], l[4], n[4], m[4];
    int i;

    xr[0] = x00;
    xv[0] = vx;
    yr[0] = y00;
    yv[0] = vy;
    for (i = 0; i < N; i++) {
        k[0] = xv[i] * h;
        l[0] = X(xr[i], yr[i]) * h;
        m[0] = yv[i] * h;
        n[0] = Y(xr[i], yr[i]) * h;
        k[1] = (xv[i] + 0.5 * l[0]) * h;
        l[1] = X(xr[i] + 0.5 * k[0], yr[i] + 0.5 * m[0]) * h;
        m[1] = (yr[i] + 0.5 * n[0]) * h;
        n[1] = Y(xr[i] + 0.5 * k[0], yr[i] + 0.5 * m[0]) * h;
        k[2] = (xv[i] + 0.5 * l[1]) * h;
        l[2] = X(xr[i] + 0.5 * k[1], yr[i] + 0.5 * m[1]) * h;
        m[2] = (yr[i] + 0.5 * n[1]) * h;
        n[2] = Y(xr[i] + 0.5 * k[1], yr[i] + 0.5 * m[1]) * h;
        k[3] = (xv[i] + l[2]) * h;
        l[3] = X(xr[i] + k[2], yr[i] + m[2]) * h;
        m[3] = (yr[i] + n[2]) * h;
        n[3] = Y(xr[i] + k[2], yr[i] + m[2]) * h;
        xr[i+1] = xr[i] + (1.0 / 6.0) * (k[0] + 2*k[1] + 2*k[2] + k[3]);
        xv[i+1] = xv[i] + (1.0 / 6.0) * (l[0] + 2*l[1] + 2*l[2] + l[3]);
        yr[i+1] = yr[i] + (1.0 / 6.0) * (m[0] + 2*m[1] + 2*m[2] + m[3]);
        yv[i+1] = yv[i] + (1.0 / 6.0) * (n[0] + 2*n[1] + 2*n[2] + n[3]);
    }
}

void updateConcentration(double* concentration, double D, double dt, double dx, double dy) {
    double* newConcentration = new double[WINDOW_WIDTH * WINDOW_HEIGHT];

    for (int x = 1; x < WINDOW_WIDTH - 1; ++x) {
        for (int y = 1; y < WINDOW_HEIGHT - 1; ++y) {
            int index = y * WINDOW_WIDTH + x;

            double laplacian = concentration[index - WINDOW_WIDTH] + concentration[index + WINDOW_WIDTH] + concentration[index - 1] + concentration[index + 1] - 4 * concentration[index];
            newConcentration[index] = concentration[index] + D * dt * laplacian / (dx * dy);
        }
    }

    memcpy(concentration, newConcentration, WINDOW_WIDTH * WINDOW_HEIGHT * sizeof(double));
    delete[] newConcentration;
}

int main() {

    //Параметры диффузии
    double D = 1.1;
    double dt = 0.6;
    double dx = 0.1;
    double dy = 0.1; 
    double intensity = 0.9;


    double* concentration = new double[WINDOW_WIDTH * WINDOW_HEIGHT];
    for (int i = 0; i < WINDOW_WIDTH * WINDOW_HEIGHT; ++i) {
        concentration[i] += 1.0;
    }

    xr = (double*)malloc(N * sizeof(double));
    xv = (double*)malloc(N * sizeof(double));
    yr = (double*)malloc(N * sizeof(double));
    yv = (double*)malloc(N * sizeof(double));


    for (int i = 0; i < N; i++) {
        // Обновление системы
        rk4();

        // Обновление массива концентрации на основе результатов runge()
        updateConcentration(concentration, D, dt, dx, dy);

        // Сохранение результатов в файл
        std::ofstream file;
        file.open("output.txt", std::ios_base::app);
        for (int j = 0; j < WINDOW_WIDTH * WINDOW_HEIGHT; ++j) {
            file << concentration[j] << " ";
        }
        file << "\n";
        file.close();
    }


    delete[] concentration;
    free(xr);
    free(xv);
    free(yr);
    free(yv);

    return 0;
}