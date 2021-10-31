#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <valarray>
using namespace std;

void Progonka(vector<vector<double>> A, double begin, double end, double h)
{
    size_t size = A.size();
    vector<double> P(size + 1), Q(size + 1), y(size);
    //обойду неудобные случаи вне цикла: начало и конец
    // начало
    P[0] = 0;
    Q[0] = 0;
    double a = 0, b = A[0][0], c = A[0][1], d = A[0][size];
    P[1] = -c / (b + a * P[0]);
    Q[1] = (d - a * Q[0]) / (b + a * P[0]);
    //обход циклом P[2] -> P[size-1] и аналогично для Q
    for (size_t i = 1; i < size - 1; ++i)
    {
        a = A[i][i - 1];
        b = A[i][i];
        c = A[i][i + 1];
        d = A[i][size];
        P[i + 1] = -c / (b + a * P[i]);
        Q[i + 1] = (d - a * Q[i]) / (b + a * P[i]);
    }
    //конец матрицы
    a = A[size - 1][size - 2];
    b = A[size - 1][size - 1];
    c = 0;
    d = A[size - 1][size];
    P[size] = -c / (b + a * P[size - 1]);
    Q[size] = (d - a * Q[size - 1]) / (b + a * P[size - 1]);
    cout << endl
         << "results:" << endl;
    y[size - 1] = Q[size];
    for (int i = size - 2; i >= 0; --i)
        y[i] = Q[i + 1] + P[i + 1] * y[i + 1];
    for (double i = 0; begin + i * h <= end; ++i){
        cout << "y(" << begin + i*h << ") = " << y[i] << endl;
    }
}

double zero(double x) { return 0; }

void FiniteDifferenceMethod(double a, double b, double h,
    double (*K)(double), double (*L)(double), double (*M)(double), double (*F)(double),
    double R, double S, double T, double V, double W, double Z) {
    size_t n = (b - a) / h + 1;
    vector<double> x;
    for (unsigned i = 0; i < n; ++i) {
        x.push_back(a + i * h);
    }
    vector<vector<double>> matrix;
    matrix.resize(n);
    for (auto& line : matrix) line.assign(n+1, 0);
    matrix[0][0] = -R / h + S;
    matrix[0][1] = R / h;
    matrix[0][n] = T;
    for (unsigned i = 1; i + 1 < n; ++i) {
        matrix[i][i - 1] = K(x[i]) / pow(h, 2) - L(x[i]) / (2 * h);
        matrix[i][i] = -2 * K(x[i]) / pow(h, 2) + M(x[i]);
        matrix[i][i + 1] = K(x[i]) / pow(h, 2) + L(x[i]) / (2 * h);
        matrix[i][n] = F(x[i]);
    }
    matrix[n-1][n - 2] = V / h;
    matrix[n-1][n-1] = - V / h - W;
    matrix[n-1][n] = -Z;
    for(auto line : matrix){
        for(auto el : line){
            cout << el << " ";
        }
        cout << endl;
    }
    Progonka(matrix, a, b, h);
}

double K(double x) { return 2*x*x; }
double L(double x) { return x; }
double M(double x) { return 1; }
double F(double x) { return 2*sqrt(x); }

int main()
{
    double a = 1, b = 10, step = 2.25, R = 0, S = 1, T = 2, V = 0, W = 1, Z = 2*sqrt(10);
    ;
    FiniteDifferenceMethod(a, b, step,
        K, L, M, F,
        R, S, T,
        V, W, Z);
    
}