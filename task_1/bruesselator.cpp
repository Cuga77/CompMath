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


constexpr int WINDOW_WIDTH = 500;
constexpr int WINDOW_HEIGHT = 500;
constexpr int size = 2 * WINDOW_WIDTH * WINDOW_HEIGHT;

constexpr double A = 3.5;
constexpr double B = 1.3;
constexpr double B1 = (B + 1);

constexpr double Dd = (double)(1.5E-5);
constexpr double H = (double)(1E-3);

// compile: g++ -O3 bruesselator.cpp -lsfml-graphics -lsfml-window -lsfml-system -fopenmp && ./a.out

double X(double x, double y) {
    return A - B1*x + x*x*y;
}

double Y(double x, double y) {
    return B * x - x*x*y;
}

void rk4(double *X, double *t, double h, void (*f)(double *X, double *Xdot), double *k1, double *k2, double *k3, double *k4) {
    f(X, k1);
    for (int i = 0; i < size; i++) {
        k1[i] *= h;
        k2[i] = X[i] + k1[i] * 0.5;
    }
    f(k2, k2);
    for (int i = 0; i < size; i++) {
        k2[i] *= h;
        k3[i] = X[i] + k2[i] * 0.5;
    }
    f(k3, k3);
    for (int i = 0; i < size; i++) {
        k3[i] *= h;
        k4[i] = X[i] + k3[i];
    }
    f(k4, k4);
    for (int i = 0; i < size; i++) {
        k4[i] *= h;
        X[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
    }
    *t += h;
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
            x[(i * WINDOW_HEIGHT + j) * 2] = (double)rand() / RAND_MAX;
            x[(i * WINDOW_HEIGHT + j) * 2 + 1] = (double)rand() / RAND_MAX;
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
                Xdot[i] = X[i] * X[i + 1] / 2 + A - B1 * X[i] - X[i];
                Xdot[i + 1] = -X[i] * X[i + 1] / 2 + B * X[i];
            }
        }, k1, k2, k3, k4);
        for (int i = 0; i < WINDOW_WIDTH; i++) {
            for (int j = 0; j < WINDOW_HEIGHT; j++) {
                sf::Color color = sf::Color(255 * x[(i * WINDOW_HEIGHT + j) * 2], 255 * x[(i * WINDOW_HEIGHT + j) * 2 + 1], 0);
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