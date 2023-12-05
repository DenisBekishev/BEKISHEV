#include <iostream>
#include <cmath>
using namespace std;

const int MAX_SIZE = 100;

void gaussElimination(double A[MAX_SIZE][MAX_SIZE], double B[MAX_SIZE], int n, double x[MAX_SIZE]) {
    for (int i = 0; i < n; i++) {
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        for (int k = i; k < n; k++) {
            swap(A[maxRow][k], A[i][k]);
        }
        swap(B[maxRow], B[i]);

        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            B[k] -= factor * B[i];
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        x[i] = B[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }

    
    double residual[MAX_SIZE];
    double maxResidual = 0.0;
    for (int i = 0; i < n; i++) {
        residual[i] = B[i];
        for (int j = 0; j < n; j++) {
            residual[i] -= A[i][j] * x[j];
        }
        if (abs(residual[i]) > maxResidual) {
            maxResidual = abs(residual[i]);
        }
    }


    cout << "Solution vector:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    cout << "Residual vector:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "residual[" << i << "] = " << residual[i] << endl;
    }


}

int main() {
    int n;
    cout << "Enter the size of the system: ";
    cin >> n;

    double A[MAX_SIZE][MAX_SIZE];
    double B[MAX_SIZE];
    double x[MAX_SIZE];

    cout << "Enter the augmented matrix:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
        cin >> B[i];
    }

    
    double A0[MAX_SIZE][MAX_SIZE];
    double b0[MAX_SIZE];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A0[i][j] = A[i][j];
        }
        b0[i] = B[i];
    }

    gaussElimination(A, B, n, x);

    
    double b_new[MAX_SIZE];
    for (int i = 0; i < n; i++) {
        b_new[i] = 0.0;
        for (int j = 0; j < n; j++) {
            b_new[i] += A0[i][j] * x[j];
        }
    }

    
    double x_new[MAX_SIZE];
    gaussElimination(A0, b_new, n, x_new); 

    
    double delta = 0.0;
    for (int i = 0; i < n; i++) {
        double error = abs(x[i] - x_new[i]) / abs(x[i]);
        if (error > delta) {
            delta = error;
        }
    }

    cout << "Relative error using your method: " << delta << endl;

    return 0;
}
