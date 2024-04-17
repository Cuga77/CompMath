#include <iostream>
#include <cstring>
#include <cmath>
#include <fstream>

#include <SFML/Graphics.hpp>

// compile: g++ -O3 brusselator.cpp -lsfml-graphics -lsfml-window -lsfml-system -fopenmp && ./a.out


// Размеры окна
const int WINDOW_WIDTH = 100;
const int WINDOW_HEIGHT = 100;

// Параметры брюсселятора
double a = 1.0;
double b = 4.0;
double b1 = b + 1;
double h = 0.000001;
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

    #pragma omp parallel for
    for (int i = 0; i < WINDOW_WIDTH; i++) {
        for (int j = 0; j < WINDOW_HEIGHT; j++) {
            double half_k1vx = k1vx / 2;
            double half_k1vy = k1vy / 2;

            k2x = (vx[i][j] + half_k1vx) * dt;
            k2y = (vy[i][j] + half_k1vy) * dt;
            k2vx = X(x[i][j] + k1x / 2, y[i][j] + k1y / 2) * dt;
            k2vy = Y(x[i][j] + k1x / 2, y[i][j] + k1y / 2) * dt;

            double half_k2vx = k2vx / 2;
            double half_k2vy = k2vy / 2;

            k3x = (vx[i][j] + half_k2vx) * dt;
            k3y = (vy[i][j] + half_k2vy) * dt;
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

void compute_diffusion(double** arr, double D, double dt, double dx, double dy) {
    double** temp = new double*[WINDOW_WIDTH];
    for (int i = 0; i < WINDOW_WIDTH; i++) {
        temp[i] = new double[WINDOW_HEIGHT];
        memcpy(temp[i], arr[i], WINDOW_HEIGHT * sizeof(double));
    }

    #pragma omp parallel for
    for (int i = 1; i < WINDOW_WIDTH - 1; i++) {
        for (int j = 1; j < WINDOW_HEIGHT - 1; j++) {
            double diff_x = (temp[i - 1][j] - 2 * temp[i][j] + temp[i + 1][j]) / (dx * dx);
            double diff_y = (temp[i][j - 1] - 2 * temp[i][j] + temp[i][j + 1]) / (dy * dy);
            arr[i][j] += D * dt * (diff_x + diff_y);
        }
    }

    for (int i = 0; i < WINDOW_WIDTH; i++) {
        delete[] temp[i];
    }
    delete[] temp;
}

int main() {
    // Создание окна SFML
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Brusselator Visualization");

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
    
    while (window.isOpen()) {
        // Обработка событий
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // Обновление системы
        #pragma omp parallel for
        for (int i = 0; i < WINDOW_WIDTH; i++) {
            for (int j = 0; j < WINDOW_HEIGHT; j++) {
                rk4(x, y, vx, vy, dt);
                compute_diffusion(x, D, dt, dx, dy);
                compute_diffusion(y, D, dt, dx, dy);
                // compute_diffusion(vx, D, dt, dx, dy);
                // compute_diffusion(vy, D, dt, dx, dy);
            }
        }

        // Очистка окна
        window.clear();

        // Визуализация данных
        for (int i = 0; i < WINDOW_WIDTH; i++) {
            for (int j = 0; j < WINDOW_HEIGHT; j++) {
                sf::RectangleShape pixel(sf::Vector2f(1, 1));
                pixel.setPosition(i, j);
                pixel.setFillColor(sf::Color(x[i][j] * 255, y[i][j] * 255, 0));
                window.draw(pixel);
            }
        }
    // Отображение окна
    window.display();
    window.clear();
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