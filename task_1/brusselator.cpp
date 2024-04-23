#include <iostream>
#include <cmath>
#include <fstream>
#include <chrono>
#include <cstring>
#include <SFML/Graphics.hpp>

constexpr int WINDOW_WIDTH = 600;
constexpr int WINDOW_HEIGHT = 600;

constexpr double A = 0.6;
constexpr double B = 2.9;
constexpr double B1 = B + 1;

constexpr double Dd = 0.9;
constexpr double DT = 0.000192;
constexpr double DX = 0.015;
constexpr double DY = 0.028;

// compile: g++ -O3 brusselator.cpp -lsfml-graphics -lsfml-window -lsfml-system -fopenmp && ./a.out

double X(double x, double y) {
    return B - B1 * x + x * x * y;
}

double Y(double x, double y) {
    return B * x - x * x * y;
}

void rk4(double** x, double** y, double** vx, double** vy, double dt) {
    double k1x, k2x, k3x, k4x;
    double k1y, k2y, k3y, k4y;
    double k1vx, k2vx, k3vx, k4vx;
    double k1vy, k2vy, k3vy, k4vy;

    #pragma omp parallel for
    for (int i = 0; i < WINDOW_WIDTH; i++) {
        for (int j = 0; j < WINDOW_HEIGHT; j++) {
            k1x = vx[i][j] * dt;
            k1y = vy[i][j] * dt;
            k1vx = X(x[i][j], y[i][j]) * dt;
            k1vy = Y(x[i][j], y[i][j]) * dt;

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

void compute_diffusion(double** arr, double** temp, double D, double dt, double dx, double dy) {
    for (int i = 0; i < WINDOW_WIDTH; ++i) {
        memcpy(temp[i], arr[i], WINDOW_HEIGHT * sizeof(double));
    }

    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double Ddt = D * dt;

    #pragma omp parallel for
    for (int i = 1; i < WINDOW_WIDTH - 1; i++) {
        for (int j = 1; j < WINDOW_HEIGHT - 1; j++) {
            double diff_x = (temp[i - 1][j] - 2 * temp[i][j] + temp[i + 1][j]) / dx2;
            double diff_y = (temp[i][j - 1] - 2 * temp[i][j] + temp[i][j + 1]) / dy2;
            arr[i][j] += Ddt * (diff_x + diff_y);
        }
    }
}

int main() {
    double D = Dd;
    double dt = DT;
    double dx = DX;
    double dy = DY; 
    double **x, **y, **vx, **vy;

    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Brusselator");   
    x = new double*[WINDOW_WIDTH];
    y = new double*[WINDOW_WIDTH];
    vx = new double*[WINDOW_WIDTH];
    vy = new double*[WINDOW_WIDTH];
    double** temp = new double*[WINDOW_WIDTH];
    for (int i = 0; i < WINDOW_WIDTH; i++) {
        x[i] = new double[WINDOW_HEIGHT];
        y[i] = new double[WINDOW_HEIGHT];
        vx[i] = new double[WINDOW_HEIGHT];
        vy[i] = new double[WINDOW_HEIGHT];
        temp[i] = new double[WINDOW_HEIGHT];
        for (int j = 0; j < WINDOW_HEIGHT; j++) {
            x[i][j] = (double)rand() / RAND_MAX; // случайное число от 0 до 1
            y[i][j] = (double)rand() / RAND_MAX;
            vx[i][j] = (double)rand() / RAND_MAX;
            vy[i][j] = (double)rand() / RAND_MAX;
        }
    }

    sf::Texture texture;
    texture.create(WINDOW_WIDTH, WINDOW_HEIGHT);
    sf::Sprite sprite(texture);
    sf::Image image;
    image.create(WINDOW_WIDTH, WINDOW_HEIGHT);

    u_int16_t step_x = int(WINDOW_WIDTH / 3);
    u_int16_t step_y = int(WINDOW_HEIGHT / 3);
    
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        #pragma omp parallel for
        for (int i = 0; i < WINDOW_WIDTH; i += step_x) {
            for (int j = 0; j < WINDOW_HEIGHT; j += step_y) {
                rk4(x, y, vx, vy, dt);
                if (i % 15 == 0 && j % 25 == 0) {
                    compute_diffusion(x, x, D, dt, dx, dy);
                    compute_diffusion(y, y, D, dt, dx, dy);
                    compute_diffusion(vx, vx, D, dt, dx, dy);
                    compute_diffusion(vy, vy, D, dt, dx, dy); 
                }
            }
        }
        #pragma omp parallel for
        for (int i = 0; i < WINDOW_WIDTH; i++) {
            for (int j = 0; j < WINDOW_HEIGHT; j++) {
                sf::Color color(x[i][j]*255, y[i][j]*255, 90);
                image.setPixel(i, j, color);
            }
        }
        
        texture.update(image);
        window.clear();
        window.draw(sprite);
        window.display();
    }



    for (int i = 0; i < WINDOW_WIDTH; i++) {
        delete[] x[i];
        delete[] y[i];
        delete[] vx[i];
        delete[] vy[i];
        delete[] temp[i];
    }
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
    delete[] temp;

    return 0;
    
}