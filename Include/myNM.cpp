/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/

#include "myNM.h"

double myFunc(const double _x) {

	return  _x * _x * _x;;

}
void printVec(double* _vec, int _row)
{
	for (int i = 0; i < _row; i++)
		printf("Vector[%d] = %.1f \n", i, _vec[i]);
	printf("\n");
}

void FileoutVec1_csv(const char* _filename, double* _vec1, int _row) {
	FILE* _fp;
	_fp = fopen(_filename, "wt");
	if (_fp == NULL) {
		printf("Fail");
		return;
		}
	
	for (int i = 0; i < _row; i++) {
		fprintf(_fp, "%f\n", _vec1[i]);
	}
	fclose(_fp);
}

void FileoutVec2_csv(const char* _filename, double* _vec1,double* _vec2, int _row) {
	FILE* _fp;
	_fp = fopen(_filename, "wt");
	if (_fp == NULL) {
		printf("Fail");
		return;
	}

	for (int i = 0; i < _row; i++) {
		fprintf(_fp, "%f %f\n", _vec1[i],_vec2[i]);
	}
	fclose(_fp);
}

void gradient1D(double x[], double y[], double dydx[], int m) {
	//
	// dydx[0]    = 3-Point Forward Difference
	// dydx[1~19] = 2-Point central difference
	// dydx[20]   = 3-Point Backwared difference
	if (sizeof(x) != sizeof(y)) {
		printf("WARNING!: Length of X and Y are not equal\n");
		return;
	}
	double h = x[1] - x[0];     //time

	dydx[0] = (-3.0 * y[0] + 4.0 * y[1] - y[2]) / (2 * h);     // 3-Point Forward Difference
	for (int k = 1; k < m - 1; k++) {                        // 2-Point central Difference
		dydx[k] = ((y[k + 1] - y[k - 1]) / (2 * h));           // 3-Point Backwared difference
	}

	dydx[m - 1] = (y[m - 3] - 4.0 * y[m - 2] + 3.0 * y[m - 1]) / (2*h);

	//return dydx
	//return output value through PrintVec
}


void gradientFunc(double myFunc(const double x), double x[], double dydx[], int m) {
	
	double* y;

	y = (double*)malloc(sizeof(double) * m);
	
  for (int i =0; i<m; i++){
     y[i]= myFunc(x[i]); 	
   }

   // get gradient dydx
     gradient1D(x,y,dydx,m);
	free(y);
}
//====================================================================Integration===================================================
double IntegrateRect(double _x[], double _y[], int _m) {
	int N = _m - 1;
	double I = 0;
	for (int i = 0; i < N; i++)
		I += _y[i] * (_x[i + 1] - _x[i]);

	return I;
}

double trapz(double _x[], double _y[], int _m) {
	int N = _m - 1;

	double Tra = 0;
	for (int i = 0; i < N; i++) {
		Tra = Tra + (0.5) * ((_y[i] + _y[i + 1]) * (_x[i + 1] - _x[i]));
	}
	return Tra;
}

double Simpson13(double myFunc(const double _x), double _a, double _b, int _m) {  // (function, initial value, final value, iteration number
	
	register Count i;

	
	double h = (_b - _a) / _m;

	double result = myFunc(_a) + myFunc(_b) + (4 * myFunc(_b - h));

	LOOP(i,_m-2,1,2){
		double xi = _a + i * h;
		result = result + (4 * myFunc(xi) + 2 * myFunc(xi + h));

	}
	return result * (h / 3);

}

double Simpson38(double myFunc(const double x), double _a, double _b, int _m) {
	register Count i;
	double h      = (_b - _a) / _m;
	double result = myFunc(_a) + myFunc(_b);

	LOOP(i, _m - 1, 1, 3) {
		double xi = _a + i * h;
		result = result + 3*(myFunc(xi) + myFunc(xi + h));

	}

	LOOP(k, _m - 2, 3, 3) {
		double xi = _a + k * h;
		result = result + (2 * myFunc(xi));
	}
	
	return result * ((3*h) / 8 );
}

//========================================================1st order ODE VIP =======================================================

void odeEU(double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h) {
	register Count i;
	double N        = (_tf - _t0) / _h;
	double t[101]  = { 0.0, };

	t[0]   = _t0;
	_y[0]  = 0.0;

	LOOP(i, N, 0, 1) {
		t[i + 1] = t[i] + _h;
		_y[i + 1] = _y[i] + myFunc(t[i], _y[i]) * _h;
		/*printf("%f\n", _y[i]);*/
	}

}

void odeEM(double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h) {
	register Count i;
	double N       = (_tf - _t0) / _h;
	double t[101]  = { 0.0, };
	double Slope1  = 0.0;
	double Slope2  = 0.0;
	double yEU     = 0.0;
	t[0]  = _t0;
	_y[0] = 0.0;

	LOOP(i, N, 0, 1) {
		t[i + 1]  = t[i] + _h;
		Slope1    = myFunc(t[i], _y[i]);
		yEU       = _y[i] + Slope1 * _h;
		Slope2    = myFunc(t[i], yEU);
		_y[i + 1] = _y[i] + (Slope1 + Slope2) * _h / 2;
	}
}

void odeRK2(double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h) {
	register Count i;
	double N = (_tf - _t0) / _h;
	double t[101] = { 0,0, };
	double K1 = 0.0;
	double K2 = 0.0;
	double alpha = 1.0;
	double C2 = 1 / (2 * alpha);
	double C1 = 1 - C2;
	t[0]  = _t0;
	_y[0] = 0.0;

	LOOP(i, N, 0, 1) {
		t[i + 1] = t[i] + _h;
		K1 = myFunc(t[i], _y[i]);
		K2 = myFunc(t[i] + _h, _y[i] + _h * K1);

		_y[i + 1] = _y[i] + (C1 * K1 + C2 * K2) * _h;
	}

}

void odeRK3(double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h) {
	register Count i;
	double N = (_tf - _t0) / _h;
	double t[101] = { 0,0, };

	double alpha2 = 1 / 2;
	double alpha3 = 1.0;

	double K1 = 0.0;
	double K2 = 0.0;
	double K3 = 0.0;
	t[0] = _t0;
	_y[0] = 0.0;

	LOOP(i, N, 0, 1) {
		t[i + 1] = t[i] + _h;
		K1 = myFunc(t[i], _y[i]);
		K2 = myFunc(t[i] + _h/2, _y[i] + _h * K1/2);
		K3 = myFunc(t[i] + _h, _y[i] - K1 * _h + 2 * K2 * _h);
		_y[i + 1] = _y[i] + _h / 6 * (K1 + 4 * K2 + K3);
		
	}

}


void ode(double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h, int method,int _row) {
	switch (method) {
	case 1:
		printf("Euler Method\n");
		odeEU(myFunc, _y, _t0, _tf, _h);
		FileoutVec1_csv("Euler_method.csv", _y, _row);
		break;
	case 2:
		printf("Euler Modified Method\n");
		odeEM(myFunc, _y, _t0, _tf, _h);
		FileoutVec1_csv("Euler_Modified_method.csv", _y, _row);
		break;
	case 3:
		printf("Runge Kutta 2nd\n");
		odeRK2(myFunc, _y, _t0, _tf, _h);
		FileoutVec1_csv("Runge_Kutta_2nd_method.csv", _y, _row);
		break;
	case 4:
		printf("Runge Kutta 2nd\n");
		odeRK3(myFunc, _y, _t0, _tf, _h);
		FileoutVec1_csv("Runge_Kutta_3rd_method.csv", _y, _row);
		break;
	}
}
//==================================================================2nd ODE IVP====================================================================//

void sys2RK2 (void myFunc(const double _t, const double _Y[], double _dYdt[]), double _y1[], double _y2[], double _t0, double _tf, double _h, double _y1_init, double _y2_init) {

	register Count i;



	double N             = (_tf - _t0) / _h;
	double t[101]        = { 0.0, };
	double K1[2]         = { 0.0, };   //[K1_y, k1_z]
	double K2[2]         = { 0.0, };    //[K2_y, K2_z]
	double Yin_k1[2]     = { 0.0, };    //[y   , y_dot] k1
	double Yin_k2[2]     = { 0.0, };   //[y   , y_dot] k2
	double K1_y = 0.0;
	double K1_z = 0.0;
	double k2_y = 0.0;
	double K2_z = 0.0;
	double C2            = 0.5;
	double C1            = 0.5;
	t[0]                 = _t0;
	_y1[0]               = _y1_init;
	_y2[0]               = _y2_init;

	LOOP(i, N, 0, 1) {
		t[i + 1] = t[i] + _h;
		Yin_k1[0] = _y1[i];       // myfunc 2nd order diff equ(t, y[], dydt[](slope->K)-> input y, y_dot
		Yin_k1[1] = _y2[i];
		myFunc(t[i], Yin_k1,K1);           //K1이 dydt의 매개변수로 들어간다 
		double K1_y = K1[0];              //K1_y
		double K1_z = K1[1];              //K1_z
		printf("%f\n", K1_y);
		Yin_k2[0] = _y1[i] + K1_y * _h;
		Yin_k2[1] = _y2[i] + K1_z * _h;

		K2[0] = K1[0] + K1_y * _h;
		K2[1] = K1[1] + K1_z * _h;
		myFunc(t[i] + _h, Yin_k2, K2);

		double K2_y = K2[0];               //K2_y
		double K2_z = K2[1];               //K2_z
                 
		_y1[i + 1] = _y1[i] + 0.5 * (K1_y + K2_y) * _h;
		_y2[i + 1] = _y2[i] + 0.5 * (K1_z + K2_z) * _h;

	}

	
}