#include "matrix_algebra.h"

void MatrixAlgebra::Set(vector<vector<double>> A_, vector<double>x_, vector<double>y_)
{
	A = A_;		x = x_;		 y = y_;
}
double Deter(vector<vector<double>> Matrix) {

	double deter = 1.0;
	for (size_t k = 0; k + 1 < Matrix.size(); k++) {
		for (size_t i = k + 1; i + 1 < Matrix[k].size(); i++) {
			Matrix[i][k] = Matrix[i][k] / Matrix[k][k];
			for (size_t j = k + 1; j < Matrix.size(); j++)
				Matrix[i][j] = Matrix[i][j] - Matrix[i][k] * Matrix[k][j];
		}
	}
	for (size_t i = 0; i < Matrix.size(); i++) {
		for (size_t j = 0; j < Matrix[i].size(); j++) {
			if (i == j)
				deter *= Matrix[i][j];
		}
	}
	return deter;
}
double MatrixAlgebra::Sum(vector<double>x_, vector<double>y_, int k_) {

	double sum = 0.0;
	if (k_ == 0)
		return siz;
	for (int i = 1; i < siz; i++) {
		sum += pow(x_[i], k_) * y_[i];
	}
	return sum;
}
vector<vector<double>>MatrixAlgebra::Mult_A(vector<vector<double>>mA, vector<vector<double>>mB) {

	vector<vector<double>>new_A(pow_k, vector <double>(pow_k, 0));

	for (int i = 0; i < pow_k; i++) {
		for (int j = 0; j < pow_k; j++) {
			new_A[i][j] = 0.0;
			for (int h = 0; h < pow_k; h++)
				new_A[i][j] += mA[i][h] * mB[h][j];
		}
	}
	return new_A;
}
vector<vector<double>>MatrixAlgebra::Get_A() {

	srand(time(NULL));
	cout << "enter x" << endl;

	for (int i = 0; i < x.size(); i++)
		cin >> x[i];

	for (int i = 0; i < pow_k; i++) {
		for (int j = 0; j < pow_k; j++) {
			A[i][j] = Sum(x, y, i + j);
		}
	}
	return A;
}
vector<double>MatrixAlgebra::System_solution(vector<vector<double>>) {

	srand(time(NULL));
	vector<double> b(pow_k, 0);

	A = Get_A();
	c.resize(pow_k);

	cout << "enter y" << endl;
	for (int i = 0; i < siz; i++)
		cin >> y[i];

	for (int i = 1; i < pow_k; i++)
		b[i] = Sum(x, y, i);

	for (int i = 0; i < siz; i++)
		b[0] += y[i];

	for (int i = 0; i < pow_k; i++) {
		cout << b[i] << " ";
	}
	cout << endl;

	A = Rev_A(A);
	c = Ab(A, b);

	return c;
}
vector<double>MatrixAlgebra::Ab(vector<vector<double>> A_, vector<double>b_) {

	vector<double> c_;
	int nrow = A_.size();
	int ncol = A_[0].size();

	if (b_.size() != ncol)
		cout << "Error" << endl;
	else {
		c_.resize(nrow);

		for (int i = 0; i < nrow; i++) {
			c_[i] = 0;
			for (int j = 0; j < ncol; j++) {
				c_[i] += A_[i][j] * b_[j];
			}
		}
	}
	return c_;
}
vector<vector<double>>MatrixAlgebra::Rev_A(vector<vector<double>>_A) {

	double   cm = 0, s = 0, t, temp, alpha = 0, beta = 0, gamma = 0, zeta, eps, EPS = 1e-8, SQR;
	short  steps = 0, MAX_STEPS = 40;
	temp = Deter(_A);

	if (temp == 0)
		cout << " Error !Determinant = 0!" << endl;
	else {
		cout << "_A" << endl;

		for (int i = 0; i < pow_k; i++) {
			for (int j = 0; j < pow_k; j++) {
				cout << _A[i][j] << " ";
			}
			cout << endl;
		}

		vector<vector<double>> M(pow_k, vector <double>(pow_k, 0));
		vector<vector<double>> RevM(pow_k, vector <double>(pow_k, 0));
		vector<vector<double>> V(pow_k, vector <double>(pow_k, 0));
		vector<vector<double>> Vt(pow_k, vector <double>(pow_k, 0));
		vector<vector<double>> R(pow_k, vector <double>(pow_k, 0));
		vector<vector<double>> Rt(pow_k, vector <double>(pow_k, 0));
		vector<double> D(pow_k, 0);
		vector<vector<double>> Ut(pow_k, vector <double>(pow_k, 0));
		vector<vector<double>> Mult(pow_k, vector <double>(pow_k, 0));
		vector<vector<double>> S(pow_k, vector <double>(pow_k, 0));
		vector<vector<double>> S1(pow_k, vector <double>(pow_k, 0));

		for (int i = 0; i < pow_k; i++) {
			D[i] = 0;
			for (int j = 0; j < pow_k; j++) {

				if (i == j)
					V[i][j] = 1;
				else
					V[i][j] = 0;

				RevM[i][j] = 0;
				R[i][j] = 0;
			}
		}
		while (steps <= MAX_STEPS) {
			for (int j = 1; j < pow_k; j++) {
				for (int k = 0; k < j; k++) {

					alpha = 0;
					beta = 0;
					gamma = 0;

					for (int i = 0; i < pow_k; i++) {
						alpha += pow(_A[i][k], 2);
						beta += pow(_A[i][j], 2);
						gamma += _A[i][k] * _A[i][j];
					}

					if (gamma == 0)
						t = 0;
					else {

						zeta = (beta - alpha) / (2 * gamma);
						t = sign(zeta) / (abs(zeta) + sqrt(1 + pow(zeta, 2)));
					}

					cm = 1 / sqrt(1 + t * t);
					s = t * cm;

					for (int i = 0; i < pow_k; i++) {

						t = _A[i][k];
						_A[i][k] = cm * t - s * _A[i][j];
						_A[i][j] = s * t + cm * _A[i][j];
						t = V[i][k];
						V[i][k] = cm * t - s * V[i][j];
						V[i][j] = s * t + cm * V[i][j];
					}
				}
			}
			for (int i = 0; i < pow_k; i++) {
				for (int m = 0; m < pow_k; m++) {

					Ut[m][i] = _A[i][m];
					Vt[m][i] = V[i][m];
				}
			}

			S = Mult_A(_A, Vt);
			Mult = Mult_A(Ut, _A);
			eps = 0;

			for (int i = 0; i < pow_k; i++) {
				for (int m = 0; m < pow_k; m++) {
					if (i != m) {
						eps += pow(Mult[i][m], 2);
					}
				}
			}
			steps++;

			if (eps < EPS) {
				break;
			}
		}
		if (steps > MAX_STEPS) {
			cout << "error" << "\n";
		}
		else {
			for (int j = 0; j < pow_k; j++) {
				SQR = 0;
				for (int i = 0; i < pow_k; i++) {

					SQR += pow(_A[i][j], 2);
				}
				D[j] = sqrt(SQR);
			}

			for (int j = 0; j < pow_k; j++)
				D[j] = 1.0 / D[j];

			for (int i = 0; i < pow_k; i++) {
				for (int j = 0; j < pow_k; j++) {

					_A[i][j] = _A[i][j] * D[j];
					Vt[i][j] = V[j][i];
				}
			}

			S1 = Mult_A(V, Vt);

			for (int k = 0; k < pow_k; k++) {
				for (int m = 0; m < pow_k; m++) {

					R[k][m] = 0;
					R[k][m] += _A[k][m] * D[m];
				}
			}

			for (int i = 0; i < pow_k; i++)
				for (int j = 0; j < pow_k; j++)
					Rt[i][j] = R[j][i];

			RevM = Mult_A(V, Rt);
		}

		_A = RevM;
		cout << "Rev_A" << endl;
		for (int i = 0; i < pow_k; i++) {
			for (int j = 0; j < pow_k; j++) {
				cout << _A[i][j] << " ";
			}
			cout << endl;
		}

		return _A;
	}
	exit(EXIT_FAILURE);
	return _A;
}
double MatrixAlgebra::sign(double G) {
	if (G > 0)
		G = 1;
	if (G < 0)
		G = -1;
	if (G == 0)
		G = 0;
	return G;
}