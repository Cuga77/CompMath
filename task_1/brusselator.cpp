#include <iostream>
#include <cmath>
#include <fstream>
#include <cstring>
#include <SFML/Graphics.hpp>

constexpr int WINDOW_WIDTH = 800;
constexpr int WINDOW_HEIGHT = 800;
constexpr int size = 2 * WINDOW_WIDTH * WINDOW_HEIGHT;

constexpr double A = 1.3;
constexpr double B = 7.1;
constexpr double B1 = (B + 1);

constexpr double H = 0.004;
constexpr double Dx = 0.00007;
constexpr double Dy = 0.00002;

// compile: g++ -O3 brusselator.cpp -lsfml-graphics -lsfml-window -lsfml-system  && ./a.out


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
            x[(j * WINDOW_WIDTH + i) * 2] = ((double)rand() / RAND_MAX) * 
                (1 + 0.1 * ((double)rand() / RAND_MAX - 0.1));
            x[(j * WINDOW_WIDTH + i) * 2 + 1] = ((double)rand() / RAND_MAX) * 
                (1 + 0.1 * ((double)rand() / RAND_MAX - 0.1));
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
                
                //границы не включены
                // double y_plus_h = (i + WINDOW_HEIGHT < size) ? X[i + WINDOW_HEIGHT] : 0;
                // double y_minus_h = (i - WINDOW_HEIGHT >= 0) ? X[i - WINDOW_HEIGHT] : 0;
                // double dx2_y = (y_plus_h - 2 * X[i] + y_minus_h) / (H * H);

                // double laplacian_x = dx2_x + dx2_y;

                // double x_plus_h_y = X[i + 2 + 1];
                // double x_minus_h_y = X[i - 2 + 1];
                // double dx2_x_y = (x_plus_h_y - 2 * X[i + 1] + x_minus_h_y) / (H * H);

                // double y_plus_h_y = (i + WINDOW_HEIGHT < size) ? X[i + WINDOW_HEIGHT + 1] : 0;
                // double y_minus_h_y = (i - WINDOW_HEIGHT >= 0) ? X[i - WINDOW_HEIGHT + 1] : 0;
                // double dx2_y_y = (y_plus_h_y - 2 * X[i + 1] + y_minus_h_y) / (H * H);


                double x_plus_h = X[(i + 2) % size];
                double x_minus_h = X[(i - 2 + size) % size];
                double dx2_x = (x_plus_h - 2 * X[i] + x_minus_h) / (H * H);

                double y_plus_h = X[(i + 2 * WINDOW_WIDTH) % size];
                double y_minus_h = X[(i - 2 * WINDOW_WIDTH + size) % size];
                double dx2_y = (y_plus_h - 2 * X[i] + y_minus_h) / (H * H);

                double laplacian_x = dx2_x + dx2_y;

                double x_plus_h_y = X[(i + 2) % size + 1]; 
                double x_minus_h_y = X[(i - 2 + size) % size + 1];
                double dx2_x_y = (x_plus_h_y - 2 * X[i + 1] + x_minus_h_y) / (H * H);

                double y_plus_h_y = X[(i + 2 * WINDOW_WIDTH) % size + 1];
                double y_minus_h_y = X[(i - 2 * WINDOW_WIDTH + size) % size + 1]; 
                double dx2_y_y = (y_plus_h_y - 2 * X[i + 1] + y_minus_h_y) / (H * H);

                double laplacian_y = dx2_x_y + dx2_y_y;

                Xdot[i] = X[i] * X[i] * X[i + 1] * 0.5 + A - B1 * X[i] + Dx * laplacian_x;
                Xdot[i + 1] = -X[i] * X[i] * X[i + 1] * 0.5 + B * X[i] + Dy * laplacian_y;
            }
        }, k1, k2, k3, k4, temp);
        for (int i = 0; i < WINDOW_WIDTH; i++) {
            for (int j = 0; j < WINDOW_HEIGHT; j++) {
                sf::Color color(255 * x[(j * WINDOW_WIDTH + i) * 2 + 1],
                    255 * x[(j * WINDOW_WIDTH + i) * 2], 150);
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