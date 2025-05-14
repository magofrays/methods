#include <iostream>
#include <cmath>
#include <iomanip>
#include <nlohmann/json.hpp>

using namespace std;

const double C = 5; // C = sqrt(2)*omega_0*t
const double eps = 1e-12;
const int max_iter = 100000;

double f(double lambda)
{
    return std::sqrt(lambda * (lambda - 1)) +
           std::log(std::sqrt(lambda) + std::sqrt(lambda - 1)) - C;
}

double df(double lambda)
{
    return 0.5 * ((2 * lambda - 1) / std::sqrt((lambda - 1) * lambda) + 1 / (std::sqrt(lambda - 1) * std::sqrt(lambda)));
}

double newton(double lambda0)
{
    double lambda = lambda0;
    int i = 0;
    for (; i < max_iter; ++i)
    {
        double delta = f(lambda) / df(lambda);
        lambda -= delta;

        if (abs(delta) < eps)
            break;
    }
    std::cout << "newton iters: " << i << "\n";
    return lambda;
}

double chord(double a, double b)
{
    if (f(a) * f(b) >= 0)
    {
        cerr << "Неверный интервал для метода хорд!" << endl;
        return NAN;
    }
    int i = 0;
    double c;
    for (; i < max_iter; ++i)
    {
        c = a - f(a) * (b - a) / (f(b) - f(a));

        if (abs(f(c)) < eps)
            break;

        f(c) * f(a) < 0 ? b = c : a = c;
    }
    std::cout << "chord iters: " << i << "\n";
    return c;
}

int main()
{
    double lambda_newton = newton(5);
    cout << "Метод Ньютона: λ = "
         << fixed << setprecision(25) << lambda_newton
         << " (ошибка: " << fabs(f(lambda_newton)) << ")" << endl;

    double lambda_chord = chord(1, 10);
    if (!isnan(lambda_chord))
    {
        cout << "Метод хорд:     λ = "
             << fixed << setprecision(25) << lambda_chord
             << " (ошибка: " << fabs(f(lambda_chord)) << ")" << endl;
    }

    return 0;
}