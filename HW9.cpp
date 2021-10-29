#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <valarray>
using namespace std;

void Print(double a, double b, double c, double d, double P, double Q, double i)
{
    cout << endl
         << "step " << i + 1 << endl;
    cout << "a = " << a << endl
         << "b = " << b << endl
         << "c = " << c << endl
         << "d = " << d << endl;
    cout << "Q[" << i + 1 << "] = " << Q << endl
         << "P[" << i + 1 << "] = " << P << endl;
}

void Progonka(vector<valarray<double>> A)
{
    size_t size = A.size();
    vector<double> P(size + 1), Q(size + 1), y(size);
    //обойду неудобные случаи вне цикла: начало и конец
    // начало
    P[0] = 0;
    Q[0] = 0;
    cout << "P[0] = 0\nQ[0] = 0\n";
    double a = 0, b = A[0][0], c = A[0][1], d = A[0][size];
    P[1] = -c / (b + a * P[0]);
    Q[1] = (d - a * Q[0]) / (b + a * P[0]);
    Print(a, b, c, d, P[1], Q[1], 0);
    //обход циклом P[2] -> P[size-1] и аналогично для Q
    for (size_t i = 1; i < size - 1; ++i)
    {
        a = A[i][i - 1];
        b = A[i][i];
        c = A[i][i + 1];
        d = A[i][size];
        P[i + 1] = -c / (b + a * P[i]);
        Q[i + 1] = (d - a * Q[i]) / (b + a * P[i]);
        Print(a, b, c, d, P[i + 1], Q[i + 1], i);
    }
    //конец матрицы
    a = A[size - 1][size - 2];
    b = A[size - 1][size - 1];
    c = 0;
    d = A[size - 1][size];
    P[size] = -c / (b + a * P[size - 1]);
    Q[size] = (d - a * Q[size - 1]) / (b + a * P[size - 1]);
    Print(a, b, c, d, P[size], Q[size], size - 1);
    cout << endl
         << "results:" << endl;
    //обход для получения иксов
    y[size - 1] = Q[size];
    for (int i = size - 2; i >= 0; --i)
        y[i] = Q[i + 1] + P[i + 1] * y[i + 1];
    for (size_t x = 0; x <= 1; x+=0.1)
        cout << "y(" << x << ") = " << y[x*10] << endl;
}

double zero(double x) { return 0; }

void FiniteDifferenceMethod(double a, double b, unsigned n,
    double (*K)(double), double (*L)(double), double (*M)(double), double (*F)(double),
    double R, double S, double T, double V, double W, double Z) {
    double h = (b - a) / (n - 1);
    vector<double> x;
    for (unsigned i = 0; i < n; ++i) {
        x.push_back(a + i * h);
    }
    vector<valarray<double>> matrix;
    matrix.resize(n+1);
    for (auto& line : matrix) {
        line.resize(n + 1);
        line.apply(zero);
    }
    matrix[0][0] = -R / h + S;
    matrix[0][1] = R / h;
    matrix[0][n] = T;
    for (unsigned i = 1; i < n; ++i) {
        double x_i = x[i];
        matrix[i][i - 1] = K(x_i) / pow(h, 2) - L(x_i) / (2 * h);
        matrix[i][i] = -2 * K(x_i) + M(x_i);
        matrix[i][i + 1] = K(x_i) / pow(h, 2) + L(x_i) / (2 * h);
        matrix[i][n] = F(x_i);
    }
    matrix[n][n - 2] = V / h;
    matrix[n][n-1] = -(V / h + W);
    matrix[n][n] = -Z;
    Progonka(matrix);
}

double K(double x) { return x*x - 1; }
double L(double x) { return x-3; }
double M(double x) { return -1; }
double F(double x) { return 0; }

int main()
{
    double a = 0, b = 1, step = 0.1, R = 1, S = 0, T = 0, V = 1, W = 1, Z = -0.75;
    ;
    FiniteDifferenceMethod(a, b, (b - a) / step,
        K, L, M, F,
        R, S, T, V, W, Z);
    
}