#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>
#include <SFML/Graphics.hpp>


// dx/dt = x^2y/2 + a - bx -x
// dy/dt = -x^2y/2 + bx

// with diffusion
// dx/dt = x^2y/2 + a - bx -x + D(d^2x/dx^2 + d^2x/dy^2)
// dy/dt = -x^2y/2 + bx + D(d^2y/dx^2 + d^2y/dy^2)


constexpr int WINDOW_WIDTH = 800;
constexpr int WINDOW_HEIGHT = 800;
constexpr int size = 2 * WINDOW_WIDTH * WINDOW_HEIGHT;

constexpr double A = 1.3;
constexpr double B = 3.5;
constexpr double B1 = (B + 1);

constexpr double H = 0.006;
constexpr double D = 0.00045;

// compile: g++ -O3 brusselator.cpp -lsfml-graphics -lsfml-window -lsfml-system -fopenmp && ./a.out


void rk4(double *X, double *t, double h, void (*f)(double *X, double *Xdot), double *k1, double *k2, double *k3, double *k4, double *Xtemp) {
    f(X, k1);
    for (int i = 0; i < size; i++) {
        k1[i] *= h;
        Xtemp[i] = X[i] + 0.5 * k1[i];
    }
    *t += 0.5 * h;
    f(Xtemp, k2);
    for (int i = 0; i < size; i++) {
        k2[i] *= h;
        Xtemp[i] = X[i] + 0.5 * k2[i];
    }
    f(Xtemp, k3);
    for (int i = 0; i < size; i++) {
        k3[i] *= h;
        Xtemp[i] = X[i] + k3[i];
    }
    *t += 0.5 * h;
    f(Xtemp, k4);
    for (int i = 0; i < size; i++) {
        k4[i] *= h;
    }
    for (int i = 0; i < size; i++) {
        X[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
    }
    *t += 0.5 * h;
}

int main() {
    double t = 0;
    double* x = new double[size];
    double* temp = new double[size];
    double* k1 = new double[size];
    double* k2 = new double[size];
    double* k3 = new double[size];
    double* k4 = new double[size];

    srand(time(NULL));
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Bruesselator"); 
    sf::Texture texture;
    texture.create(WINDOW_WIDTH, WINDOW_HEIGHT);
    sf::Sprite sprite(texture);
    sf::Image image;
    image.create(WINDOW_WIDTH, WINDOW_HEIGHT);  
    
    for (int i = 0; i < WINDOW_WIDTH; i++) {
        for (int j = 0; j < WINDOW_HEIGHT; j++) {
            x[(i * WINDOW_HEIGHT + j) * 2] = ((double)rand() / RAND_MAX) * 
                (1 + 0.1 * ((double)rand() / RAND_MAX - 0.5));
            x[(i * WINDOW_HEIGHT + j) * 2 + 1] = ((double)rand() / RAND_MAX) * 
                (1 + 0.1 * ((double)rand() / RAND_MAX - 0.5));
        }
    }
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        rk4(x, &t, H, [](double *X, double *Xdot) {
            for (int i = 0; i < size; i += 2) {
                double dx = (X[(i + 2) % size] - 2 * X[i] + X[(i - 2 + size) % size]) / (H * H);
                double dy = (X[(i + WINDOW_HEIGHT) % size] - 2 * X[i] + X[(i - WINDOW_HEIGHT + size) % size]) / (H * H); 
                double laplacian_x = dx + dy;

                double dx_y = (X[(i + 2) % size + 1] - 2 * X[i + 1] + X[(i - 2 + size) % size + 1]) / (H * H);
                double dy_y = (X[(i + WINDOW_HEIGHT) % size + 1] - 2 * X[i + 1] + X[(i - WINDOW_HEIGHT + size) % size + 1]) / (H * H); 
                double laplacian_y = dx_y + dy_y;
                
                // double laplacian_x = (X[(i + 2) % size] - 2 * X[i] + X[(i - 2 + size) % size] + X[(i + WINDOW_HEIGHT) % size] - 2 * X[i] + X[(i - WINDOW_HEIGHT + size) % size]) / (H * H); //laplacian = (f(x + h) - 2f(x) + f(x - h)) / h^2
                // double laplacian_y = (X[(i + 2) % size + 1] - 2 * X[i + 1] + X[(i - 2 + size) % size + 1] + X[(i + WINDOW_HEIGHT) % size + 1] - 2 * X[i + 1] + X[(i - WINDOW_HEIGHT + size) % size + 1]) / (H * H); //laplacian = (f(x + h) - 2f(x) + f(x - h)) / h^2

                Xdot[i] = X[i] * X[i] * X[i + 1] * 0.5 + A - B1 * X[i] + D * laplacian_x;
                Xdot[i + 1] = -X[i] * X[i] * X[i + 1] * 0.5 + B * X[i] + D * laplacian_y;
            }
        }, k1, k2, k3, k4, temp);
        for (int i = 0; i < WINDOW_WIDTH; i++) {
            for (int j = 0; j < WINDOW_HEIGHT; j++) {
                sf::Color color(255 * x[(i * WINDOW_HEIGHT + j) * 2 + 1],
                    255 * x[(i * WINDOW_HEIGHT + j) * 2], 0);
                image.setPixel(i, j, color);
            }
        }
        texture.update(image);
        sprite.setTexture(texture, true);
        window.clear();
        window.draw(sprite);
        window.display();
    }

    delete[] k1;
    delete[] k2;
    delete[] k3;
    delete[] k4;
    delete[] temp;
    delete[] x;
    return 0;
}