#include <iostream>
#include <cmath>

using namespace std;


const int n = 2;


void F(double x[n], double result[n]) {
    result[0] = x[0] - x[1] - 6 * log(x[0]) - 1;
    result[1] = x[0] - 3 * x[1] - 6 * log(x[0]) - 2;
}


void J(double x[n], double Jacobian[n][n]) {
    Jacobian[0][0] = 1 - 6 / x[0];
    Jacobian[0][1] = -1;
    Jacobian[1][0] = 1 - 6 / x[0];
    Jacobian[1][1] = -3;
}

int main() {
    
    double x[n] = { 0.5, 0.2 };

    
    double epsilon1 = 1e-9;
    double epsilon2 = 1e-9;

    
    int maxIterations = 100;

    
    for (int k = 0; k < maxIterations; ++k) {
        
        double residual[n];
        F(x, residual);

        
        double Jacobian[n][n];
        J(x, Jacobian);

        
        double deltaX[n] = {
            -(residual[0] * Jacobian[1][1] - residual[1] * Jacobian[0][1]) /
            (Jacobian[0][0] * Jacobian[1][1] - Jacobian[1][0] * Jacobian[0][1]),

            -(residual[1] * Jacobian[0][0] - residual[0] * Jacobian[1][0]) /
            (Jacobian[0][0] * Jacobian[1][1] - Jacobian[1][0] * Jacobian[0][1])
        };

        
        x[0] += deltaX[0];
        x[1] += deltaX[1];

        
        double delta1 = sqrt(residual[0] * residual[0] + residual[1] * residual[1]);
        double delta2 = sqrt(deltaX[0] * deltaX[0] + deltaX[1] * deltaX[1]);
        cout << "Iteration " << k + 1 << ": delta1 = " << delta1 << ", delta2 = " << delta2 << endl;

        
        if (delta1 < epsilon1 && delta2 < epsilon2) {
            cout << "Converged to solution: x = {" << x[0] << ", " << x[1] << "}" << endl;
            break;
        }

        
        if (k >= maxIterations - 1) {
            cout << "Maximum number of iterations reached without convergence." << endl;
            break;
        }
    }

    return 0;
}
