#include <SFML/Graphics.hpp>
#include <iostream>
#include <cstring>
#include <cmath>

// compile: g++ -o brusselator brusselator.cpp -lsfml-graphics -lsfml-window -lsfml-system

// Размеры окна
u_int64_t WINDOW_WIDTH = 150;
u_int64_t WINDOW_HEIGHT = 150;

// Параметры брюсселятора
double a = 1.0;
double b = 4.0;
double b1 = b+1;
double h = 0.00001;
double x00 = 0.0, y00 = 0.0, vx = 2.0, vy = 1.5;
int N = 100;

// Параметры диффузии
double D = 0.9;
double dt = 0.1;
double dx = 0.1;
double dy = 0.1; 
double intensity = 0.9;

//Параметры изображения
double frame = 0.0;

double* xr, *xv, *yr, *yv;

double X(double x, double y) {
    return a - b1 * x + x * x * y;
}

double Y(double x, double y) {
    return b * x - x * x * y;
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
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Brusselator Diffusion with SFML");

    sf::Texture texture;
    texture.create(WINDOW_WIDTH, WINDOW_HEIGHT);
    sf::Sprite sprite(texture);


    double* concentration = new double[WINDOW_WIDTH * WINDOW_HEIGHT];
    for (int i = 0; i < WINDOW_WIDTH * WINDOW_HEIGHT; ++i) {
        concentration[i] += 1.0;
    }

    xr = (double*)malloc(N * sizeof(double));
    xv = (double*)malloc(N * sizeof(double));
    yr = (double*)malloc(N * sizeof(double));
    yv = (double*)malloc(N * sizeof(double));


    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        //Обновление системы
        runge();

        // Обновление массива концентрации на основе результатов runge()
        updateConcentration(concentration, D, dt, dx, dy);

        // Обновление массива концентрации на основе результатов runge()
        for (int x = 0; x < WINDOW_WIDTH; ++x) {
            for (int y = 0; y < WINDOW_HEIGHT; ++y) {
                int index = y * WINDOW_WIDTH + x;

                if (x == 0 || y == 0 || x == WINDOW_WIDTH - 1 || y == WINDOW_HEIGHT - 1) {
                    // Условия Дирихле: функция равна нулю на границе
                    concentration[index] = 0;
                } else {
                    // Обновляем концентрацию
                    concentration[index] = xr[index] + xv[index] + yr[index] + yv[index];
                }
            }
        }

        // Визуализация
        sf::Uint8* pixels = new sf::Uint8[WINDOW_WIDTH * WINDOW_HEIGHT * 4];
        for (int i = 0; i < WINDOW_WIDTH * WINDOW_HEIGHT; ++i) {
            double value = concentration[i];
            double normalizedValue = value / 255.0; 
            sf::Uint8 color = static_cast<sf::Uint8>(value * 255);
            // sf::Uint8 animatedColor = static_cast<sf::Uint8>((color + frame) % 256);
            sf::Uint8 animatedColor = static_cast<sf::Uint8>((color + sf::Uint8(128 * (1 + sin(frame / 50.0)))) % 256);   

            // Используем цветовую карту для отображения различных уровней концентрации
            sf::Color colorMap = sf::Color::Black;
            if (value < 0.5) {
                colorMap = sf::Color::Blue;    
            } else {
                colorMap = sf::Color::Red;
            }

            pixels[i * 4 + 0] = colorMap.r;
            pixels[i * 4 + 1] = colorMap.g;
            pixels[i * 4 + 2] = colorMap.b;
            pixels[i * 4 + 3] = 255; 
        }
        texture.update(pixels);
        delete[] pixels;

        frame++;

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