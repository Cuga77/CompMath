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


constexpr int WINDOW_WIDTH = 1000;
constexpr int WINDOW_HEIGHT = 1000;

constexpr double A = 4.6;
constexpr double B = 1.2;
constexpr double B1 = (B + 1);

constexpr double Dd = (double)(1.5E-10);
constexpr double H = (double)(1.5E-10);

// compile: g++ -O3 bruesselator.cpp -lsfml-graphics -lsfml-window -lsfml-system -fopenmp && ./a.out

double X(double x, double y) {
    return A - B1*x + x*x*y;
}

double Y(double x, double y) {
    return B * x - x*x*y;
}
//todo: Переделать
// //диффузия ПОСЛЕ
// void diffusion(double** x, double** y, double D, double H, double dx, double dy) {
//     double** temp = new double*[WINDOW_WIDTH];
//     #pragma omp parallel for
//     for (int i = 0; i < WINDOW_WIDTH; i++) {
//         temp[i] = new double[WINDOW_HEIGHT];
//         for (int j = 0; j < WINDOW_HEIGHT; j++) {
//             int i_minus = (i > 0) ? i - 1 : WINDOW_WIDTH - 1;
//             int i_plus = (i < WINDOW_WIDTH - 1) ? i + 1 : 0;
//             int j_minus = (j > 0) ? j - 1 : WINDOW_HEIGHT - 1;
//             int j_plus = (j < WINDOW_HEIGHT - 1) ? j + 1 : 0;
//             double dx2_x = (x[i_minus][j] - 2 * x[i][j] + x[i_plus][j]) / (dx * dx);
//             double dy2_x = (x[i][j_minus] - 2 * x[i][j] + x[i][j_plus]) / (dy * dy);
//             temp[i][j] = x[i][j] + dt * X(x[i][j], y[i][j], dx2_x, dy2_x);
//             double dx2_y = (y[i_minus][j] - 2 * y[i][j] + y[i_plus][j]) / (dx * dx);
//             double dy2_y = (y[i][j_minus] - 2 * y[i][j] + y[i][j_plus]) / (dy * dy);
//             temp[i][j] = y[i][j] + dt * Y(x[i][j], y[i][j], dx2_y, dy2_y);
//             x[i][j] = temp[i][j];
//             y[i][j] = temp[i][j];
//         }
//         free(temp[i]);
//     }
//     free(temp);
// }

//todo: Переделать
//Диффузия по 4 соседним точкам
// void diffusion(double** x, double** y, double D, double dt, double dx, double dy) {
//     double** temp = new double*[WINDOW_WIDTH];
//     #pragma omp parallel for
//     for (int i = 0; i < WINDOW_WIDTH; i++) {
//         temp[i] = new double[WINDOW_HEIGHT];
//         for (int j = 0; j < WINDOW_HEIGHT; j++) {
//             int i_minus = (i > 0) ? i - 1 : i;
//             int i_plus = (i < WINDOW_WIDTH - 1) ? i + 1 : i;
//             int j_minus = (j > 0) ? j - 1 : j;
//             int j_plus = (j < WINDOW_HEIGHT - 1) ? j + 1 : j;
//             double dx2_x = (x[i_minus][j] - 2 * x[i][j] + x[i_plus][j]) / (dx * dx);
//             double dy2_x = (x[i][j_minus] - 2 * x[i][j] + x[i][j_plus]) / (dy * dy);
//             temp[i][j] = x[i][j] + dt * (X(x[i][j], y[i][j], dx2_x, dy2_x) + D * (dx2_x + dy2_x));
//             // double dx2_x = (x[i][j] - 2 * x[i][j] + x[i][j]) / (dx * dx);
//             // double dy2_x = (x[i][j] - 2 * x[i][j] + x[i][j]) / (dy * dy);
//             // temp[i][j] = x[i][j] + dt * (X(x[i][j], y[i][j], dx2_x, dy2_x) + D * (dx2_x + dy2_x));
//             double dx2_y = (y[i_minus][j] - 2 * y[i][j] + y[i_plus][j]) / (dx * dx);
//             double dy2_y = (y[i][j_minus] - 2 * y[i][j] + y[i][j_plus]) / (dy * dy);
//             temp[i][j] = y[i][j] + dt * (Y(x[i][j], y[i][j], dx2_y, dy2_y) + D * (dx2_y + dy2_y));
//             // double dx2_y = (y[i][j] - 2 * y[i][j] + y[i][j]) / (dx * dx);
//             // double dy2_y = (y[i][j] - 2 * y[i][j] + y[i][j]) / (dy * dy);
//             // temp[i][j] = y[i][j] + dt * (Y(x[i][j], y[i][j], dx2_y, dy2_y) + D * (dx2_y + dy2_y));
//             x[i][j] = temp[i][j];
//             y[i][j] = temp[i][j];
//         }
//         free(temp[i]);
//     }
//     free(temp);
// }
//todo: Переделать
//Дисперсия с дискретным оператором Лапласа
// void diffusion_with_Laplace_discrete(double x[], double H, ) {
//     double** temp_x = new double*[WINDOW_WIDTH];
//     double** temp_y = new double*[WINDOW_WIDTH];
//     double dx2 = dx * dx;
//     double dy2 = dy * dy;
//     #pragma omp parallel for
//     for (int i = 0; i < WINDOW_WIDTH; i++) {
//         temp_x[i] = new double[WINDOW_HEIGHT];
//         temp_y[i] = new double[WINDOW_HEIGHT];
//         for (int j = 0; j < WINDOW_HEIGHT; j++) {
//             int i_minus = (i > 0) ? i - 1 : 1;
//             int i_plus = (i < WINDOW_WIDTH - 1) ? i + 1 : 1;
//             int j_minus = (j > 0) ? j - 1 : 1;
//             int j_plus = (j < WINDOW_HEIGHT - 1) ? j + 1 : 1;
//             double laplacian_x = (x[i_minus][j] + x[i_plus][j] + x[i][j_minus] + x[i][j_plus] - 4 * x[i][j]) / dx2;
//             double laplacian_y = (y[i_minus][j] + y[i_plus][j] + y[i][j_minus] + y[i][j_plus] - 4 * y[i][j]) / dy2;
//             temp_x[i][j] = x[i][j] + H * D * laplacian_x;
//             temp_y[i][j] = y[i][j] + H * D * laplacian_y;
//         }
//     }
//     for (int i = 0; i < WINDOW_WIDTH; i++) {
//         for (int j = 0; j < WINDOW_HEIGHT; j++) {
//             x[i][j] = temp_x[i][j];
//             y[i][j] = temp_y[i][j];
//         }
//         delete[] temp_x[i];
//         delete[] temp_y[i];
//     }
//     delete[] temp_x;
//     delete[] temp_y;
// }


//// todo: RK4  не должна знать о показаниях x и y
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

    temp_x = x[0] + (1.0 / 9.0) * (2 * k1x + 3 * k3x + 4 * k4x);
    temp_y = x[1] + (1.0 / 9.0) * (2 * k1y + 3 * k3y + 4 * k4y);

    x[0] = temp_x;
    x[1] = temp_y;

    *t += h;
}

//todo: Переделать
//диффузия которая была ДО
// void diffusion(double** arr, double** temp, double D, double H, double dx, double dy) {
//     for (int i = 0; i < WINDOW_WIDTH; ++i) {
//         memcpy(temp[i], arr[i], WINDOW_HEIGHT * sizeof(double));
//     }
//     double dx2 = dx * dx;
//     double dy2 = dy * dy;
//     double DH = D * H;
//     #pragma omp parallel for
//     for (int i = 0; i < WINDOW_WIDTH; i++) {
//         for (int j = 0; j < WINDOW_HEIGHT; j++) {
//             double diff_x = (temp[(i > 0 ? i - 1 : WINDOW_WIDTH - 1)][j] - 2 * temp[i][j] + temp[(i < WINDOW_WIDTH - 1 ? i + 1 : 0)][j]) / dx2;
//             double diff_y = (temp[i][(j > 0 ? j - 1 : WINDOW_HEIGHT - 1)] - 2 * temp[i][j] + temp[i][(j < WINDOW_HEIGHT - 1 ? j + 1 : 0)]) / dy2;
//             arr[i][j] += DH * (diff_x + diff_y);
//         }
//     }
// }


int main() {
    double* x = new double[WINDOW_WIDTH * WINDOW_HEIGHT * 2];
    double t = 0.0;

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
        for (int i = 0; i < WINDOW_WIDTH; i++) {
            for (int j = 0; j < WINDOW_HEIGHT; j++) {
                x[(i * WINDOW_HEIGHT + j) * 2] += (double)rand() / RAND_MAX - 0.5;
                x[(i * WINDOW_HEIGHT + j) * 2 + 1] += (double)rand() / RAND_MAX - 0.5;
                rk4(&x[(i * WINDOW_HEIGHT + j) * 2], &t, H);

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
    delete[] x;
    return 0;
}