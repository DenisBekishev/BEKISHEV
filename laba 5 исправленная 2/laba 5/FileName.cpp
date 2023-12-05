#include<iostream>
#include<cmath>

const double eps = 0.00001;

using namespace std;

double function(double x) {
    return pow(x * x * x - 1, -0.5);
}

double cube_simpson(double x, double y) {
    return pow(x * x / (1 + y * y), 2);
}

void trapezoidal(double a, double b, int n) {
    double h = (b - a) / n;
    double integral = 0.5 * (function(a) + function(b));

    for (int i = 1; i < n; i++) {
        integral += function(a + i * h);
    }

    integral *= h;
    cout << "Integral by Trapezoidal rule: " << integral << endl;
}

void simpson(double a, double b, int n) {
    double h = (b - a) / n;
    double integral = function(a) + function(b);

    for (int i = 1; i < n; i++) {
        if (i % 2 != 0)
            integral += 4 * function(a + i * h);
        else
            integral += 2 * function(a + i * h);
    }

    integral *= h / 3;
    cout << "Integral by Simpson's rule: " << integral << endl;
}

void cube_simpson(double a, double b, double c, double d, int m, int n) {
    double hx = (b - a) / (2 * n);
    double hy = (d - c) / m;
    double integral = 0;

    for (int i = 0; i < 2 * n; i += 2) {
        for (int j = 0; j < m; j++) {
            integral += cube_simpson(a + i * hx, c + j * hy);
            integral += 4 * cube_simpson(a + (i + 1) * hx, c + j * hy);
            integral += cube_simpson(a + (i + 2) * hx, c + j * hy);
            integral += 4 * cube_simpson(a + i * hx, c + (j + 1) * hy);
            integral += 16 * cube_simpson(a + (i + 1) * hx, c + (j + 1) * hy);
            integral += 4 * cube_simpson(a + (i + 2) * hx, c + (j + 1) * hy);
            integral += cube_simpson(a + i * hx, c + (j + 2) * hy);
            integral += 4 * cube_simpson(a + (i + 1) * hx, c + (j + 2) * hy);
            integral += cube_simpson(a + (i + 2) * hx, c + (j + 2) * hy);
        }
    }

    integral *= (hx * hy / 9);
    cout << "Integral by Cube Simpson's rule: " << integral << endl;
}

int main() {
    double a = 1.3;
    double b = 2.621;

    int trapezoidal_n = 100; // Number of intervals for Trapezoidal rule
    int simpson_n = 100;     // Number of intervals for Simpson's rule
    int cube_simpson_m = 2;  // Number of intervals along y-axis for Cube Simpson's rule
    int cube_simpson_n = 4;  // Number of intervals along x-axis for Cube Simpson's rule

    trapezoidal(a, b, trapezoidal_n);
    cout << "--------------------------------------------------------------------" << endl;
    simpson(a, b, simpson_n);
    cout << "--------------------------------------------------------------------" << endl;
    double c = 1;
    double d = 2;
    cube_simpson(a, b, c, d, cube_simpson_m, cube_simpson_n);

    return 0;
}
