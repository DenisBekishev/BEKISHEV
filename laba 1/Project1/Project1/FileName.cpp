#include <iostream>
#include <iomanip>

using namespace std;

void matrix(double** matrix_full, double* matrix_answers, int n) {
	cout << endl << "�������: " << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << setw(10) << matrix_full[i][j];
		}
		cout << setw(10) << matrix_answers[i];
		cout << endl;
	}
	cout << endl;
}

double* ve�torNevyazki(double** matrix_full, double* matrix_answers, double* answers, int n) {
	double* tempSum = new double[n] {};
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			tempSum[i] += matrix_full[i][j] * answers[j];
		}
	}
	for (int i = 0; i < n; i++) {
		tempSum[i] = tempSum[i] - matrix_answers[i];
	}
	return tempSum;
}

double* find_x(double** matrix_full, double* matrix_answers, int n) {
	for (int i = 0; i < n; i++) {
		double max_number = abs(matrix_full[i][i]);
		int max_column = i;
		for (int j = i + 1; j < n; j++) {
			if (abs(matrix_full[j][i]) > max_number) {
				max_number = abs(matrix_full[j][i]);
				max_column = j;
			}
		}
		if (max_number == 0) {
			cout << "Bad type of matrix all main elements == 0";
			return 0;
		}
		swap(matrix_full[i], matrix_full[max_column]);
		swap(matrix_answers[i], matrix_answers[max_column]);
		matrix(matrix_full, matrix_answers, n);
		double divider = matrix_full[i][i];
		for (int k = i; k < n; k++) {
			matrix_full[i][k] /= divider;
		}
		matrix_answers[i] /= divider;
		double* save_string = new double[n + 1];
		for (int j = i; j < n; j++) {
			save_string[j] = matrix_full[i][j];
		}
		save_string[n] = matrix_answers[i];
		for (int j = i + 1; j < n; j++) {
			matrix_answers[j] -= save_string[n] * matrix_full[j][i];
			for (int k = n - 1; k >= i; k--) {
				matrix_full[j][k] -= save_string[k] * matrix_full[j][i];
			}
		}
		delete[] save_string;
	}
	matrix(matrix_full, matrix_answers, n);

	double* answers = new double[n];
	for (int i = n - 1; i >= 0; i--) {
		if (i + 1 == n) {
			answers[i] = matrix_answers[i] / matrix_full[i][i];
		}
		else {
			double sum = 0;
			for (int j = n - 1; j > i; j--) {
				sum += matrix_full[i][j] * answers[j];
			}
			answers[i] = (matrix_answers[i] - sum) / matrix_full[i][i];
		}
	}
	return answers;
}

double pogreshnost(double* answers, double* new_answers, int n) {
	double* result = new double[n];
	for (int i = 0; i < n; i++) {
		result[i] = new_answers[i] - answers[i];
	}
	double max_up = result[0];
	for (int i = 1; i < n; i++) {
		if (max_up < result[i]) {
			max_up = result[i];
		}
	}
	double max_down = answers[0];
	for (int i = 1; i < n; i++) {
		if (max_down < answers[i]) {
			max_down = answers[i];
		}
	}
	return max_up / max_down;
}

int main() {
	setlocale(LC_CTYPE, "Russian");
	cout << "������� ������ �������" << endl << "n = ";
	int n;
	cin >> n;
	cout << endl;
	double** matrix_full = new double* [n];
	for (int i = 0; i < n; i++) {
		matrix_full[i] = new double[n];
	}
	for (int i = 0; i < n; i++) {
		cout << "������� " << i + 1 << " ������" << endl;
		for (int j = 0; j < n; j++) {
			cin >> matrix_full[i][j];
		}
		cout << endl;
	}
	double** matrix_full_copy = new double* [n];
	for (int i = 0; i < n; i++) {
		matrix_full_copy[i] = new double[n];
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			matrix_full_copy[i][j] = matrix_full[i][j];
		}
	}
	cout << "������� ������ �" << endl;
	double* matrix_answers = new double[n];
	for (int i = 0; i < n; i++) {
		cin >> matrix_answers[i];
	}
	double* matrix_answers_copy = new double[n];
	for (int i = 0; i < n; i++) {
		matrix_answers_copy[i] = matrix_answers[i];
	}
	matrix(matrix_full, matrix_answers, n);

	double* answers = find_x(matrix_full, matrix_answers, n);
	if (answers == 0) return 0;
	for (int i = 0; i < n; i++) {
		cout << "X" << i + 1 << " = " << answers[i] << " ";
	}
	cout << endl;
	cout << "������ �������: ";
	double* vectorNevyaz = ve�torNevyazki(matrix_full_copy, matrix_answers_copy, answers, n);
	for (int i = 0; i < n; i++) {
		cout << vectorNevyaz[i] << " ";
	}
	double norma = vectorNevyaz[0];
	for (int i = 0; i < n; i++) {
		if (vectorNevyaz[i] > norma) {
			norma = vectorNevyaz[i];
		}
	}
	cout << endl << "����� �������: " << norma << endl;

	double* new_answers = find_x(matrix_full_copy, answers, n);
	cout << endl << "����� X:" << endl;
	for (int i = 0; i < n; i++) {
		cout << "X" << i + 1 << " = " << new_answers[i] << " ";
	}
	cout << endl << "������������� �����������: " << pogreshnost(answers, new_answers, n);

	
	delete[] answers;
	delete[] matrix_answers;
	delete[] matrix_answers_copy;
	for (int i = 0; i < n; i++) {
		delete[] matrix_full[i];
	}
	delete[] matrix_full;
	for (int i = 0; i < n; i++) {
		delete[] matrix_full_copy[i];
	}
	delete[] matrix_full_copy;
}