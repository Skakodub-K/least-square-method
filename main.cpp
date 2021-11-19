#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include "matrix_algebra.h"
#include <windows.h>

using namespace std;
int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	MatrixAlgebra mnk;
	int n, k;

	while (true) {
		cout << "! k <= n !" << endl;
		cout << "¬ведите число экспериментальных точек(n) :" << endl;
		cin >> n;
		cout << "¬ведите степень полинома(k) : " << endl;
		cin >> k;
		if (k <= n)
			break;
	}

	mnk.pow_k = k + 1;
	mnk.siz = n;

	vector<double> x(n, 0);
	vector<double> y(n, 1);
	vector<vector<double>> A(k + 1, vector <double>(k + 1, 0));
	vector<double> z(k + 1, 0);

	mnk.Set(A, x, y);
	z = mnk.System_solution(A);

	//reverse(z.begin(), z.end()); чтобы сперва шел коэфициент со старшей степенью
	cout << "Z: " << endl;
	for (size_t i = 0; i < k + 1; i++) {
		cout << z[i] << " ";
	}
	printf("\n");
	return 0;
}
