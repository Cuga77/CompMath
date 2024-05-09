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


constexpr int WINDOW_WIDTH = 600;
constexpr int WINDOW_HEIGHT = 600;

constexpr double A = 4.6;
constexpr double B = 1.2;
constexpr double B1 = (B + 1);

constexpr double Dd = (double)(1.5E-5);
constexpr double H = (double)(1.5E-4);

// compile: g++ -O3 bruesselator.cpp -lsfml-graphics -lsfml-window -lsfml-system -fopenmp && ./a.out

double X(double x, double y) {
    return A - B1*x + x*x*y;
}

double Y(double x, double y) {
    return B * x - x*x*y;
}


// //todo: Доделать, чтобы работало без segfault
// void rk4(double *X, double *t, double h, void (*f)(double *X, double *Xdot), double *k1, double *k2, double *k3, double *k4, double *Xtemp) {
//     int size = 2 * WINDOW_WIDTH * WINDOW_HEIGHT;

//     f(X, k1);
//     for (int i = 0; i < size; i++) {
//         Xtemp[i] = X[i] + 0.5 * h * k1[i];
//     }

//     f(Xtemp, k2);
//     for (int i = 0; i < size; i++) {
//         Xtemp[i] = X[i] + 0.5 * h * k2[i];
//     }

//     f(Xtemp, k3);
//     for (int i = 0; i < size; i++) {
//         Xtemp[i] = X[i] + h * k3[i];
//     }

//     f(Xtemp, k4);
//     for (int i = 0; i < size; i++) {
//         X[i] = X[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
//     }
// }

void rk4 (double x[], double *t, double h) {
    double k1x, k2x, k3x, k4x;
    double k1y, k2y, k3y, k4y;
    double temp_x, temp_y; 
    k1x = h * X(x[0], x[1]);
    k1y = h * Y(x[0], x[1]);
    k2x = h * X(x[0] + 0.5 * k1x, x[1] + 0.5 * k1y);
    k2y = h * Y(x[0] + 0.5 * k1x, x[1] + 0.5 * k1y);
    k3x = h * X(x[0] + 0.5 * k2x, x[1] + 0.5 * k2y);
    k3y = h * Y(x[0] + 0.5 * k2x, x[1] + 0.5 * k2y);
    k4x = h * X(x[0] + k3x, x[1] + k3y);
    k4y = h * Y(x[0] + k3x, x[1] + k3y);
    temp_x = x[0] + (1.0 / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x);
    temp_y = x[1] + (1.0 / 6.0) * (k1y + 2 * k2y + 2 * k3y + k4y);
    x[0] = temp_x;
    x[1] = temp_y;
    *t += h;
}


int main() {
    double t = 0;
    double* x = new double[WINDOW_WIDTH * WINDOW_HEIGHT * 2];
    double* temp = new double[WINDOW_WIDTH * WINDOW_HEIGHT * 2];
    double* k1 = new double[WINDOW_WIDTH * WINDOW_HEIGHT * 2];
    double* k2 = new double[WINDOW_WIDTH * WINDOW_HEIGHT * 2];
    double* k3 = new double[WINDOW_WIDTH * WINDOW_HEIGHT * 2];
    double* k4 = new double[WINDOW_WIDTH * WINDOW_HEIGHT * 2];

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
        
        #pragma parallel for
        for (int i = 0; i < WINDOW_WIDTH; i++) {
            for (int j = 0; j < WINDOW_HEIGHT; j++) {
                // x[(i * WINDOW_HEIGHT + j) * 2] += (double)rand() / RAND_MAX - 0.0001;
                // x[(i * WINDOW_HEIGHT + j) * 2 + 1] += (double)rand() / RAND_MAX - 0.0001;
                rk4(&x[(i * WINDOW_HEIGHT + j) * 2], &t, H);
                // rk4(&x[(i * WINDOW_HEIGHT + j) * 2], &t, H, [](double *X, double *Xdot) {
                //     Xdot[0] = X[0] * X[0] * X[1] / 2 + A - B1 * X[0] - X[0];
                //     Xdot[1] = -X[0] * X[0] * X[1] / 2 + B * X[0];
                // }, k1, k2, k3, k4, temp);
                

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