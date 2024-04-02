#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>

// Размеры окна
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 600;

// Параметры брюсселятора
double a = 1.0;
double b = 3.0;
double h = 0.01;
double x00 = 0.0, y00 = 0.0, vx = 0.0, vy = 0.0;
int N = 1000;

double* xr, *xv, *yr, *yv;

double X(double x, double y) {
    return a - (b + 1) * x + pow(x, 2) * y;
}

double Y(double x, double y) {
    return b * x - pow(x, 2) * y;
}

void runge() {
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

int main() {
    // Инициализация SFML
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Brusselator Diffusion with SFML");

    // Создание текстуры для визуализации
    sf::Texture texture;
    texture.create(WINDOW_WIDTH, WINDOW_HEIGHT);
    sf::Sprite sprite(texture);

    // Создание массива для хранения значений концентрации
    double* concentration = new double[WINDOW_WIDTH * WINDOW_HEIGHT];
    // Инициализация значений концентрации
    for (int i = 0; i < WINDOW_WIDTH * WINDOW_HEIGHT; ++i) {
        concentration[i] = 0.0;
    }

    // Выделение памяти для хранения результатов
    xr = (double*)malloc(N * sizeof(double));
    xv = (double*)malloc(N * sizeof(double));
    yr = (double*)malloc(N * sizeof(double));
    yv = (double*)malloc(N * sizeof(double));

    // Главный цикл
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        // Обновление состояния системы
        runge();

        // Обновление массива концентрации на основе результатов runge()
        for (int i = 0; i < N; ++i) {
            double intensity = xr[i]; // Используем xr как интенсивность
            for (int x = 0; x < WINDOW_WIDTH; ++x) {
                for (int y = 0; y < WINDOW_HEIGHT; ++y) {
                    int index = y * WINDOW_WIDTH + x;
                    concentration[index] = intensity; // Обновляем концентрацию
                }
            }
        }

        // Визуализация
        sf::Uint8* pixels = new sf::Uint8[WINDOW_WIDTH * WINDOW_HEIGHT * 4];
        for (int i = 0; i < WINDOW_WIDTH * WINDOW_HEIGHT; ++i) {
            double value = concentration[i];
            pixels[i * 4 + 0] = static_cast<sf::Uint8>(value * 255); // Red
            pixels[i * 4 + 1] = static_cast<sf::Uint8>(value * 255); // Green
            pixels[i * 4 + 2] = static_cast<sf::Uint8>(value * 255); // Blue
            pixels[i * 4 + 3] = 255; // Alpha
        }
        texture.update(pixels);
        delete[] pixels;

        window.clear();
        window.draw(sprite);
        window.display();
    }

    delete[] concentration;
    free(xr);
    free(xv);
    free(yr);
    free(yv);

    return 0;
}