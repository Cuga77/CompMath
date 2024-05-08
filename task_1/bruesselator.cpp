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


constexpr int WINDOW_WIDTH = 400;
constexpr int WINDOW_HEIGHT = 400;

constexpr double A = 1.0;
constexpr double B = 3.5;
constexpr double B1 = (B + 1);

constexpr double Dd = (double)(1.1E-4);
constexpr double DT = 0.0001;
constexpr double DX = 0.0015;
constexpr double DY = 0.0028;

// compile: g++ -O3 bruesselator.cpp -lsfml-graphics -lsfml-window -lsfml-system -fopenmp && ./a.out

double X(double x, double y, double dx2, double dy2) {
    return A - B1*x + x*x*y + Dd * (dx2 + dy2);
}

double Y(double x, double y, double dx2, double dy2) {
    return B * x - x*x*y + Dd * (dx2 + dy2);
}

//диффузия ПОСЛЕ

void diffusion(double** x, double** y, double** new_x, double** new_y, double D, double dt, double dx, double dy) {
    #pragma omp parallel for
    for (int i = 0; i < WINDOW_WIDTH; i++) {
        for (int j = 0; j < WINDOW_HEIGHT; j++) {
            int i_minus = (i > 0) ? i - 1 : WINDOW_WIDTH - 1;
            int i_plus = (i < WINDOW_WIDTH - 1) ? i + 1 : 0;
            int j_minus = (j > 0) ? j - 1 : WINDOW_HEIGHT - 1;
            int j_plus = (j < WINDOW_HEIGHT - 1) ? j + 1 : 0;

            double dx2_x = (x[i_minus][j] - 2 * x[i][j] + x[i_plus][j]) / (dx * dx);
            double dy2_x = (x[i][j_minus] - 2 * x[i][j] + x[i][j_plus]) / (dy * dy);
            new_x[i][j] = x[i][j] + dt * X(x[i][j], y[i][j], dx2_x, dy2_x);

            double dx2_y = (y[i_minus][j] - 2 * y[i][j] + y[i_plus][j]) / (dx * dx);
            double dy2_y = (y[i][j_minus] - 2 * y[i][j] + y[i][j_plus]) / (dy * dy);
            new_y[i][j] = y[i][j] + dt * Y(x[i][j], y[i][j], dx2_y, dy2_y);
        }
    }
}

void rk4 (double** x, double** y, double** vx, double** vy, double dt) {
    double** k1x = new double*[WINDOW_WIDTH];
    double** k1y = new double*[WINDOW_WIDTH];
    double** k1vx = new double*[WINDOW_WIDTH];
    double** k1vy = new double*[WINDOW_WIDTH];
    double** k2x = new double*[WINDOW_WIDTH];
    double** k2y = new double*[WINDOW_WIDTH];
    double** k2vx = new double*[WINDOW_WIDTH];
    double** k2vy = new double*[WINDOW_WIDTH];
    double** k3x = new double*[WINDOW_WIDTH];
    double** k3y = new double*[WINDOW_WIDTH];
    double** k3vx = new double*[WINDOW_WIDTH];
    double** k3vy = new double*[WINDOW_WIDTH];
    double** k4x = new double*[WINDOW_WIDTH];
    double** k4y = new double*[WINDOW_WIDTH];
    double** k4vx = new double*[WINDOW_WIDTH];
    double** k4vy = new double*[WINDOW_WIDTH];
    for (int i = 0; i < WINDOW_WIDTH; i++) {
        k1x[i] = new double[WINDOW_HEIGHT];
        k1y[i] = new double[WINDOW_HEIGHT];
        k1vx[i] = new double[WINDOW_HEIGHT];
        k1vy[i] = new double[WINDOW_HEIGHT];
        k2x[i] = new double[WINDOW_HEIGHT];
        k2y[i] = new double[WINDOW_HEIGHT];
        k2vx[i] = new double[WINDOW_HEIGHT];
        k2vy[i] = new double[WINDOW_HEIGHT];
        k3x[i] = new double[WINDOW_HEIGHT];
        k3y[i] = new double[WINDOW_HEIGHT];
        k3vx[i] = new double[WINDOW_HEIGHT];
        k3vy[i] = new double[WINDOW_HEIGHT];
        k4x[i] = new double[WINDOW_HEIGHT];
        k4y[i] = new double[WINDOW_HEIGHT];
        k4vx[i] = new double[WINDOW_HEIGHT];
        k4vy[i] = new double[WINDOW_HEIGHT];
    }

    #pragma omp parallel for
    for (int i = 0; i < WINDOW_WIDTH; i++) {
        for (int j = 0; j < WINDOW_HEIGHT; j++) {
            k1x[i][j] = dt * X(x[i][j], y[i][j], DX, DY);
            k1y[i][j] = dt * Y(x[i][j], y[i][j], DX, DY);
            k1vx[i][j] = dt * X(vx[i][j], vy[i][j], DX, DY);
            k1vy[i][j] = dt * Y(vx[i][j], vy[i][j], DX, DY);
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < WINDOW_WIDTH; i++) {
        for (int j = 0; j < WINDOW_HEIGHT; j++) {
            k2x[i][j] = dt * X(x[i][j] + k1x[i][j] / 2, y[i][j] + k1y[i][j] / 2, DX, DY);
            k2y[i][j] = dt * Y(x[i][j] + k1x[i][j] / 2, y[i][j] + k1y[i][j] / 2, DX, DY);
            k2vx[i][j] = dt * X(vx[i][j] + k1vx[i][j] / 2, vy[i][j] + k1vy[i][j] / 2, DX, DY);
            k2vy[i][j] = dt * Y(vx[i][j] + k1vx[i][j] / 2, vy[i][j] + k1vy[i][j] / 2, DX, DY);
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < WINDOW_WIDTH; i++) {
        for (int j = 0; j < WINDOW_HEIGHT; j++) {
            k3x[i][j] = dt * X(x[i][j] + k2x[i][j] / 2, y[i][j] + k2y[i][j] / 2, DX, DY);
            k3y[i][j] = dt * Y(x[i][j] + k2x[i][j] / 2, y[i][j] + k2y[i][j] / 2, DX, DY);
            k3vx[i][j] = dt * X(vx[i][j] + k2vx[i][j] / 2, vy[i][j] + k2vy[i][j] / 2, DX, DY);
            k3vy[i][j] = dt * Y(vx[i][j] + k2vx[i][j] / 2, vy[i][j] + k2vy[i][j] / 2, DX, DY);
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < WINDOW_WIDTH; i++) {
        for (int j = 0; j < WINDOW_HEIGHT; j++) {
            k4x[i][j] = dt * X(x[i][j] + k3x[i][j], y[i][j] + k3y[i][j], DX, DY);
            k4y[i][j] = dt * Y(x[i][j] + k3x[i][j], y[i][j] + k3y[i][j], DX, DY);
            k4vx[i][j] = dt * X(vx[i][j] + k3vx[i][j], vy[i][j] + k3vy[i][j], DX, DY);
            k4vy[i][j] = dt * Y(vx[i][j] + k3vx[i][j], vy[i][j] + k3vy[i][j], DX, DY);
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < WINDOW_WIDTH; i++) {
        for (int j = 0; j < WINDOW_HEIGHT; j++) {
            x[i][j] += (k1x[i][j] + 2 * k2x[i][j] + 2 * k3x[i][j] + k4x[i][j]) / 6;
            y[i][j] += (k1y[i][j] + 2 * k2y[i][j] + 2 * k3y[i][j] + k4y[i][j]) / 6;
            vx[i][j] += (k1vx[i][j] + 2 * k2vx[i][j] + 2 * k3vx[i][j] + k4vx[i][j]) / 6;
            vy[i][j] += (k1vy[i][j] + 2 * k2vy[i][j] + 2 * k3vy[i][j] + k4vy[i][j]) / 6;
        }
    }

    diffusion(x, y, x, y, Dd, dt, DX, DY);
    // diffusion(vy, vx, vx, vy, Dd, dt, DX, DY);

    for (int i = 0; i < WINDOW_WIDTH; i++) {
        delete[] k1x[i];
        delete[] k1y[i];
        delete[] k1vx[i];
        delete[] k1vy[i];
        delete[] k2x[i];
        delete[] k2y[i];
        delete[] k2vx[i];
        delete[] k2vy[i];
        delete[] k3x[i];
        delete[] k3y[i];
        delete[] k3vx[i];
        delete[] k3vy[i];
        delete[] k4x[i];
        delete[] k4y[i];
        delete[] k4vx[i];
        delete[] k4vy[i];
    }

    delete[] k1x;
    delete[] k1y;
    delete[] k1vx;
    delete[] k1vy;
    delete[] k2x;
    delete[] k2y;
    delete[] k2vx;
    delete[] k2vy;
    delete[] k3x;
    delete[] k3y;
    delete[] k3vx;
    delete[] k3vy;
    delete[] k4x;
    delete[] k4y;
    delete[] k4vx;
    delete[] k4vy;
}


//диффузия которая была ДО

// void diffusion(double** arr, double** temp, double D, double dt, double dx, double dy) {
//     for (int i = 0; i < WINDOW_WIDTH; ++i) {
//         memcpy(temp[i], arr[i], WINDOW_HEIGHT * sizeof(double));
//     }

//     double dx2 = dx * dx;
//     double dy2 = dy * dy;
//     double Ddt = D * dt;

//     #pragma omp parallel for
//     for (int i = 0; i < WINDOW_WIDTH; i++) {
//         for (int j = 0; j < WINDOW_HEIGHT; j++) {
//             double diff_x = (temp[(i > 0 ? i - 1 : WINDOW_WIDTH - 1)][j] - 2 * temp[i][j] + temp[(i < WINDOW_WIDTH - 1 ? i + 1 : 0)][j]) / dx2;
//             double diff_y = (temp[i][(j > 0 ? j - 1 : WINDOW_HEIGHT - 1)] - 2 * temp[i][j] + temp[i][(j < WINDOW_HEIGHT - 1 ? j + 1 : 0)]) / dy2;
//             arr[i][j] += Ddt * (diff_x + diff_y);
//         }
//     }
// }

int main() {
    double D = Dd;
    double dt = DT;
    double t = 0;
    double dx = DX;
    double dy = DY; 
    double **x, **y, **vx, **vy;

    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Bruesselator");   
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
            x[i][j] = (double)rand() / RAND_MAX; // rand num 0 to 1
            y[i][j] = (double)rand() / RAND_MAX;
            vx[i][j] = (double)rand() / RAND_MAX;
            vy[i][j] = (double)rand() / RAND_MAX;
        }
    }

    std::ofstream outfile("output.txt");

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

        //Для графиков
        // #pragma omp parallel for
        // for (int i = 0; i < WINDOW_HEIGHT; i += 1) {
        //     for (int j = 0; j < WINDOW_WIDTH; j += 1) {
        //     rk4(x, y, vx, vy, dt);
        //     outfile << x[i][j] << " " << vy[i][j] << " " << t << "\n";
        //     t += dt;
        //     }
        // }

        #pragma omp parallel for
        for (int i = 0; i < WINDOW_WIDTH; i+= step_x) {
            for (int j = 0; j < WINDOW_HEIGHT; j+= step_y) {
                rk4(x, y, vx, vy, dt);
            }
        }
        #pragma omp parallel for
        for (int i = 0; i < WINDOW_WIDTH; i++) {
            for (int j = 0; j < WINDOW_HEIGHT; j++) {
                sf::Color color(std::log(x[i][j])*255, std::log(y[i][j])*255, 0);
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

    outfile.close();    

    return 0;
    
}