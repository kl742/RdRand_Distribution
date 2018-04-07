#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <conio.h>
#include <time.h>
#include <math.h>
#include <immintrin.h>
#include <malloc.h>
#include <locale.h>
static double x[257], y[256];
const double x1e = 7.69711747013104972;
const double Se = 3.9496598225815571993e-3; 
const unsigned __int64 max = 18446744073709551615;
const double x1n = 3.6541528853610088;
const double Sn = 4.92867323399e-3; 

double UniformDistribution(double a, double b)
{
	unsigned __int64 r; 
	_rdrand64_step(&r);
	return a + r*((b - a) / max);
}


double ExponentialDistSmirnov(double lambda)
{
	return (-1 / lambda)*log(UniformDistribution(0, 1));
}

void ReactExponential() {
	y[0] = exp(-x1e);
	x[256] = 0;
	x[0] = Se / y[0];
	for (unsigned i = 1; i <= 255; ++i)
	{
		x[i] = -log(y[i - 1]);
		y[i] = y[i - 1] + Se/x[i];
	}
}

double ExpZ() {
	int i = 0, idReact;
	unsigned int r;
	double e;
	for (;;) {
		_rdrand32_step(&r);
		idReact = r & 255;
		e = UniformDistribution(0, x[idReact]); 
		if (e < x[idReact + 1]) // Если попали в любой прямоугольник кроме 1го
			return e;
		if (idReact == 0) // Если попали в 1й прямоугольник
			return x1e + ExpZ(); // рекурсия, в случае попадания в 1й прямоугольник
		if (UniformDistribution(y[idReact - 1], y[idReact]) < exp(-e))
			return e;
	}
	return 0;
}

double ExponentialZiggurat(double lambda) {
	return ExpZ() / lambda;
}

double NormalCentralTheory(double m, double sigma)
{
	double S = 0, F;
	for (int i = 0;i < 50; i = i + 1)
	{
		S = S + UniformDistribution(0, 1);
	}
	F = (S - 25) / 2.041241452319;
	return m + F*sigma;
}




double NormalBoxMuller(double m, double sigma)
{
	double s, x, y, z;
	do {
		x = UniformDistribution(-1, 1);
		y = UniformDistribution(-1, 1);
		s = x*x + y*y;
	} while (s > 1 || s <= 0);
	z = x*pow((-2 * log(s) / s), 0.5);
	return m + z*sigma;
}

void ReactNormal() {
	y[0] = exp(-.5 * x1n * x1n);
	x[0] = Sn / y[0];
	x[256] = 0;
	for (unsigned i = 1; i <= 255; ++i)
	{
		x[i] = sqrt(-2 * log(y[i - 1]));
		y[i] = y[i - 1] + Sn / x[i];
	}
}

double SNZig() {
	int idReact, r;
	for (;;) {
		_rdrand32_step(&r);
		idReact = r & 255;
		double xn = UniformDistribution(0, x[idReact]);
		if (xn < x[idReact + 1])
			return ((signed)r > 0) ? xn : -xn;
		if (idReact == 0)
		{
			static double z = -1;
			double yn;
			if (z > 0)
			{
				xn = ExponentialDistSmirnov(x1n);
				z -= 0.5 * xn * xn;
			}
			if (z <= 0)
			{
				do {
					xn = ExponentialDistSmirnov(x1n);
					yn = ExponentialDistSmirnov(1);
					z = yn - 0.5 * xn * xn;
				} while (z <= 0);
			}
			xn += x1n;
			return ((signed)r > 0) ? xn : -xn;
		}
		if (UniformDistribution(y[idReact - 1], y[idReact]) < exp(-.5 * xn * xn))
			return ((signed)r > 0) ? xn : -xn;
	}
	return 0; 
}

double NormalZiggurat(double m, double sigma) {
	return m + SNZig() * sigma;
}

double CauchyDistSmirnov(double x0, double gamma)
{
	return x0 + gamma*tan(M_PI*(UniformDistribution(0,1) - 0.5));
}

double CauchyСircle(double x0, double gamma) {
	int i = 0;
	double x, y;
	do {
		x = UniformDistribution(-1, 1);
		y = UniformDistribution(-1, 1);
		i = i + 1;
	} while (x * x + y * y > 1.0 || y == 0.0);
	return x0 + gamma * x / y;
	
}

double CauchyNormal(double x0, double gamma) {
	return x0 + gamma * NormalZiggurat(0,1) / NormalZiggurat(0,1);
}

double LaplaceDist(double alpha, double beta)
{
	if (UniformDistribution(-1, 1) > 0) return beta + ExponentialDistSmirnov(alpha);
	else return beta - ExponentialDistSmirnov(alpha);
}

double LogNormal(double m, double sigma)
{
	return exp(NormalZiggurat(m, sigma));
}

double RayleighN(double sigma) 
{
	double N1 = NormalZiggurat(0, sigma);
	double N2 = NormalZiggurat(0, sigma);
	return sqrt(N1*N1 + N2*N2);
}

double RayleighE(double sigma) 
{
	return sigma * sqrt(ExponentialDistSmirnov(0.5));
}

double LogisticDistSmirnov(double m, double s)
{
	double U = UniformDistribution(0,1);
	return m + s*log(U/(1-U));
}