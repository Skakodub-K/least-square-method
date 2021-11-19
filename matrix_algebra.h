#include <vector>
#include <random>
#include <time.h>
#include <cmath>
#include <iostream>
using namespace std;

class MatrixAlgebra
{
public:
	int pow_k = 0;
	int siz = 0;
	vector<double>System_solution(vector<vector<double>>);
	void Set(vector<vector<double>>, vector<double>, vector<double>);
private:
	vector<double> x;
	vector<double> y;
	vector<double> b;
	vector<double> c;
	double sign(double G);
	vector<vector<double>> A;
	vector<vector<double>>Get_A();
	vector<double>Ab(vector<vector<double>>, vector<double>);
	vector<vector<double>>Rev_A(vector<vector<double>>);
	vector<vector<double>>Mult_A(vector<vector<double>>, vector<vector<double>>);
	double Sum(vector<double>, vector<double>, int);
};
