/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/

#include "myNM.h"

//double myFunc(const double _x) {
//
//	return  _x * _x * _x;;
//
//}
void printVec1(double* _vec, int _row)
{
	for (int i = 0; i < _row; i++)
		printf("Vector[%d] = %.10f \n", i, _vec[i]);
	printf("\n");
}

void printVec2(double* _vec1,double* _vec2, int _row)
{
	for (int i = 0; i < _row; i++)
		printf("Vector1[%d] = %.10f Vector2[%d]=%.10f \n", i, _vec1[i],i,_vec2[i]);
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
		fprintf(_fp, "%f, %f\n", _vec1[i],_vec2[i]);
	}
	fclose(_fp);
}

void FileoutVec3_csv(const char* _filename, double* _vec1, double* _vec2, double* _vec3, int _row) {
	FILE* _fp;
	_fp = fopen(_filename, "wt");
	if (_fp == NULL) {
		printf("Fail");
		return;
	}

	for (int i = 0; i < _row; i++) {
		fprintf(_fp, "%f, %f, %f\n", _vec1[i], _vec2[i], _vec3[i]);
	}
	fclose(_fp);
}

double factorial(double _x)
{
	if (_x <= 1)
		return 1;
	else
		return _x * factorial(_x - 1);
}


//////////////////////////////////Function sinTaylor : approximate sine function Taylor series)////////////////////////////////
double sinTaylor(double _x) // input unit [rad]///angle   
{
	/////////////////////////Function Parameter//////////////////////////

	int     N_max = 20;
	double  epsilon = 1e-5;


	int     N = 0;
	double  S_N = 0.0;
	double  S_N_prev;
	double  rel_chg;
	do
	{

		N = N + 1;

		S_N_prev = S_N;
		S_N = 0.0;
		for (int k = 0; k < N; k++)
		{
			S_N = S_N + (pow(-1, k) * pow(_x, 2 * k + 1) / factorial(2 * k + 1));

		}
		rel_chg = fabs((S_N - S_N_prev) / S_N_prev);
	} while (N < N_max && rel_chg >= epsilon);

	return S_N;

}

double sindTaylor(double _x) /// input unit [deg]
{

	double S_N_deg = sinTaylor(_x * D2R);

	return S_N_deg;

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

double DifThreeFor(double _t[], double _X[], int _index) {
	double k = 240.0;
	double h = (_t[_index + 1] - _t[_index])/10;
	double dTdx = (-3 * _X[_index] + 4 * _X[_index+1] - _X[_index+2]) / (2 * h); // error h^2
	return -1*k*dTdx;
}

double DifThreeBack(double _t[], double _X[], int _index) {
	double k = 240.0;
	double h= (_t[_index] - _t[_index-1]) / 10;
	double dTdx = (_X[_index - 2] - 4 * _X[_index - 1] + 3 * _X[_index]) / (2 * h); //error h^2
	printf("%f,%f,%f\n", k, h, dTdx);
	return -1 * k * dTdx;
}


//====================================================================Integration===================================================
double IntegrateRect(double _x[], double _y[], int _m) {
	int N = _m - 1;
	double I = 0;
	for (int i = 0; i < N-1; i++)
		I += _y[i] * (_x[i + 1] - _x[i]);

	return I;
}

double Integralmid(double _x[], double _y[], int _m) {
	register Count i;

	int N = _m - 1;
	double I = 0.0;

	LOOP(i, N - 1, 0, 1) {
		I = I + (_x[i + 1] - _x[i]) * (_y[i + 1] + _y[i]) * 0.5;
	}
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

double trapz_func(double myFunc(const double _x), double _a, double _b, int _m) {
	
	register Count i;
	double h = (_b - _a) / _m ;
	double Tra = 0.0;

		LOOP(i, _m-1, 0, 1) {
			double xi = _a + i * h;
			Tra = Tra + 0.5 * ((myFunc(xi) + myFunc(xi + h)) * h);
			
			//printf("%f, %f\n", xi,Tra);
		}
		return Tra;
}

double Simpson13(double _x[], double _a, double _b, int _m) {  // (function, initial value, final value, iteration number
	int N = _m - 1;	
	double h =(_b - _a) / N;
	register Count i;
	double Simpson13 = 0.0;

	Simpson13 = _x[0] + _x[N] + _x[N-1];
	
	LOOP(i, _m - 2, 1, 2) {
		Simpson13 = Simpson13 + (4 * _x[i] + 2 * _x[i+1]);

	}
	

	return Simpson13 * (h / 3);

}

double Simpson13func(double myFunc(const double _x), double _a, double _b, int _m) {  // (function, initial value, final value, iteration number
	
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

void odeEU(double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h,double _t[]) {
	register Count i;
	double N        = ((_tf - _t0) / _h) +1.0;
	//double t[101]  = { 0.0, };

	_t[0]   = _t0;
	_y[0]  = 0.0;

	LOOP(i, N-1, 0, 1) {
		_t[i + 1] = _t[i] + _h;
		_y[i + 1] = _y[i] + myFunc(_t[i], _y[i]) * _h;
		/*printf("%f\n", _y[i]);*/
	}

}

void odeEM(double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h,double _t[]) {
	register Count i;
	double N       = ((_tf - _t0) / _h) +1.0;
	//double t[101]  = { 0.0, };
	double Slope1  = 0.0;
	double Slope2  = 0.0;
	double yEU     = 0.0;
	_t[0]  = _t0;
	_y[0] = 0.0;

	LOOP(i, N-1, 0, 1) {
		_t[i + 1]  = _t[i] + _h;
		Slope1    = myFunc(_t[i], _y[i]);
		yEU       = _y[i] + Slope1 * _h;
		Slope2    = myFunc(_t[i], yEU);
		_y[i + 1] = _y[i] + (Slope1 + Slope2) * _h / 2;
		_t[i + 1] = _t[i] + _h;
	}
}

void odeRK2(double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h,double _t[]) {
	register Count i;
	double N = (_tf - _t0) / _h+1.0;
	//double t[101] = { 0,0, };
	double K1 = 0.0;
	double K2 = 0.0;
	double alpha = 1.0;
	double C2 = 1 / (2 * alpha);
	double C1 = 1 - C2;
	_t[0]  = _t0;
	_y[0] = 0.0;

	LOOP(i, N-1, 0, 1) {
		_t[i + 1] = _t[i] + _h;
		K1 = myFunc(_t[i], _y[i]);
		K2 = myFunc(_t[i] + _h, _y[i] + _h * K1);

		_y[i + 1] = _y[i] + (C1 * K1 + C2 * K2) * _h;
		_t[i + 1] = _t[i] + _h;
	}
}

void odeRK3(double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h,double _t[]) {
	register Count i;
	double N = (_tf - _t0) / _h+1.0;
	//double t[101] = { 0,0, };

	double alpha2 = 0.5;
	double alpha3 = 1.0;

	double K1 = 0.0;
	double K2 = 0.0;
	double K3 = 0.0;
	_t[0] = _t0;
	_y[0] = 0.0;

	LOOP(i, N-1, 0, 1) {
		_t[i + 1] = _t[i] + _h;
		K1 = myFunc(_t[i], _y[i]);
		K2 = myFunc(_t[i] + _h/2.0, _y[i] + _h * K1/2.0);
		K3 = myFunc(_t[i] + _h, _y[i] - K1 * _h + 2.0 * K2 * _h);
		_y[i + 1] = _y[i] + _h / 6.0 * (K1 + 4.0 * K2 + K3);
		_t[i + 1] = _t[i] + _h;
	}

}

void odeRK4(double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h,double _t[]) {
	register Count i;
	double N = (_tf - _t0) / _h+1;
	//double t[101] = { 0.0, };

	double alpha2 = 0.5;
	double alpha3 = 0.5;

	double K1 = 0.0;
	double K2 = 0.0;
	double K3 = 0.0;
	double K4 = 0.0;

	_t[0] = _t0;
	_y[0] = 0.0;

	LOOP(i, N - 1, 0, 1) {
		K1 = myFunc(_t[i], _y[i]);
		K2 = myFunc(_t[i] + 0.5 * _h, _y[i] + _h * 0.5 * K1);
		K3 = myFunc(_t[i] + 0.5 * _h, _y[i] + _h * 0.5 * K2);
		K4 = myFunc(_t[i] + _h, _y[i] + _h * K3);
		_y[i + 1] = _y[i] + _h / 6.0 * (K1 + 2 * K2 + 2 * K3 + K4);
		_t[i+1] = _t[i] + _h;
	}

}

void odePC(double myFunc(const double _t, const double _y), double _y[], double _ypc[], double _t0, double _tf, double _h, double _y0, double _t[]) {
	
	register Count i=0;
	int N_max     = 100;
	double tol    = pow(10, -5);
	double N      = (_tf - _t0) / _h + 1;
	double K1     = 0.0;
	double K2     = 0.0;
	double K3     = 0.0;
	double sub    = 0.0;
	double yEU    = 0.0;
	_y[0]         = _y0;
	_ypc[0]       = _y0;

	int Max_count = 0.0;
	
	do
	{
		Max_count++;

		K1           = myFunc(_t[i], _y[i]);
		yEU          = _y[i] + K1 * _h;
		K2           = myFunc(_t[i] + _h, yEU);
		_y[i + 1]    = _y[i] + 0.5 * (K1 + K2)*_h; //Euler Modified Method

		_ypc[0]      = _y[1];
		K3           = myFunc(_t[i] + _h, _y[i + 1]);
		_ypc[i + 1] = _y[i] + 0.5 * (K1 + K3) * _h;

		sub = fabs(_ypc[i + 1] - _y[i + 1] / _y[i + 1]);

		_t[i + 1] = _t[i] + _h;
		i++;


	} while (Max_count < N_max && sub >= tol);
}


void ode(double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h, int method, int _row, double _t[]) {
	
	switch (method) {
	case 1:
		printf("Euler Method\n");
		odeEU(myFunc, _y, _t0, _tf, _h,_t);
		FileoutVec1_csv("Euler_method.csv", _y, _row);
		break;
	case 2:
		printf("Euler Modified Method\n");
		odeEM(myFunc, _y, _t0, _tf, _h, _t);
		FileoutVec1_csv("Euler_Modified_method.csv", _y, _row);
		break;
	case 3:
		printf("Runge Kutta 2nd\n");
		odeRK2(myFunc, _y, _t0, _tf, _h, _t);
		FileoutVec1_csv("Runge_Kutta_2nd_method.csv", _y, _row);
		break;
	case 4:
		printf("Runge Kutta 3nd\n");
		odeRK3(myFunc, _y, _t0, _tf, _h, _t);
		printVec1(_y, _row);
		FileoutVec1_csv("Runge_Kutta_3rd_method.csv", _y, _row);
		break;
	case 5:
		printf("Runge Kutta 4th\n");
		odeRK4(myFunc, _y, _t0, _tf, _h, _t);
		printVec1(_y, _row);
		FileoutVec1_csv("Runge_Kutta_4th_method.csv", _y, _row);
	
	}
		
}
//==================================================================2nd ODE IVP====================================================================//

void sys2RK2 (void myFunc(const double _t, const double _Y[], double _dYdt[]), double _y1[], double _y2[], double _t0, double _tf, double _h, double _y1_init, double _y2_init, double _t[]) {

	register Count i;



	double N             = (_tf - _t0) / _h +1;
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
	_t[0]                 = _t0;
	_y1[0]               = _y1_init;
	_y2[0]               = _y2_init;

	LOOP(i, N-1, 0, 1) {
		_t[i + 1] = _t[i] + _h;
		Yin_k1[0] = _y1[i];       // myfunc 2nd order diff equ(t, y[], dydt[](slope->K)-> input y, y_dot
		Yin_k1[1] = _y2[i];
		myFunc(_t[i], Yin_k1,K1);           //K1이 dydt의 매개변수로 들어간다 
		double K1_y = K1[0];              //K1_y
		double K1_z = K1[1];              //K1_z
		
		//Slope 2
		Yin_k2[0] = _y1[i] + K1_y * _h;
		Yin_k2[1] = _y2[i] + K1_z * _h;

		myFunc(_t[i] + _h, Yin_k2, K2);

		double K2_y = K2[0];               //K2_y
		double K2_z = K2[1];               //K2_z
                 
		_y1[i + 1] = _y1[i] + 0.5 * (K1_y + K2_y) * _h;
		_y2[i + 1] = _y2[i] + 0.5 * (K1_z + K2_z) * _h;

	}

	
}

void sys2RK4(void myFunc(const double _t, const double _Y[], double _dYdt[]), double _y1[], double _y2[], double _t0, double _tf, double _h, double _y1_init, double _y2_init,double _t[]) {


	register Count i;

	double N                  = (_tf - _t0) / _h +1;
	//double _t[101]             = { 0.0, };
	double K1[2]              = { 0.0, };     //[K1_y, k1_z]
	double K2[2]              = { 0.0, };     //[K2_y, K2_z]
	double K3[2]              = { 0.0, };     //[k3_y, k3_z]
	double K4[2]              = { 0.0, };     //[k4_y, k4_z]
	double Yin_k1[2]          = { 0.0, }; //[y   , y_dot] k1
	double Yin_k2[2]          = { 0.0, }; //[y   , y_dot] k2
	double Yin_k3[2]          = { 0.0, }; //[y   , y_dot] k3
	double Yin_k4[2]          = { 0.0, }; //[y   , y_dot] k4
	double K1_y               = 0.0;
	double K1_z               = 0.0;
	double K2_y               = 0.0;
	double K2_z               = 0.0;
	double K3_y               = 0.0;
	double K3_z               = 0.0;
	double K4_y               = 0.0;
	double K4_z               = 0.0;
	double C2                 = 0.5;
	double C1                 = 0.5;
	double alpha              = 0.5;
	_t[0] = _t0;
	_y1[0] = _y1_init;
	_y2[0] = _y2_init;

	LOOP(i, N - 1, 0, 1) {
		_t[i + 1]             = _t[i] + _h;
		
		//Slope 1
		Yin_k1[0]            = _y1[i];
		Yin_k1[1]            = _y2[i];
		myFunc(_t[i], Yin_k1, K1);      // K1 ->dydt
		K1_y                 = K1[0];  // K1_y 
		K1_z                 = K1[1];  // K1_z

		// Slope 2
		Yin_k2[0]            = _y1[i] + alpha * K1_y * _h;
		Yin_k2[1]            = _y2[i] + alpha * K1_z * _h;

		myFunc(_t[i]+0.5*_h, Yin_k2, K2);

		K2_y                 = K2[0];
		K2_z                 = K2[1];

		//Slope 3

		Yin_k3[0]            = _y1[i] + alpha * K2_y * _h;
		Yin_k3[1]            = _y2[i] + alpha * K2_z * _h;

		myFunc(_t[i] + 0.5 * _h, Yin_k3, K3);
		K3_y                 = K3[0];
		K3_z                 = K3[1];

		//Slope 4

		Yin_k4[0] = _y1[i] + K3_y * _h;
		Yin_k4[1] = _y2[i] + K3_z * _h;

		myFunc(_t[i] + _h, Yin_k4, K4);
		K4_y = K4[0];
		K4_z = K4[1];

		_y1[i + 1] = _y1[i] + (K1_y + 2 * K2_y + 2 * K3_y + K4_y) * _h / 6.0;
		_y2[i + 1] = _y2[i] + (K1_z + 2 * K2_z + 2 * K3_z + K4_z) * _h / 6.0;

	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
//Matrix time = createMat(N, 1);
//for (int i = 0; i < N; i++)
//{
//
//	time.at[i][0] = t0 + i * h;
//
//}
//for (int i = 0; i < N; i++)
//{
//	ta[i] = t0 + i * h;
//
//}
//printf("/*---------------------------*/\n");
//printf("/ 2nd Order ODE : Pendulum /\n");
//printf("/*---------------------------*/\n");
//printf("ans for (a)\n");
//
//for (int i = 0; i < N; i++)
//	printf("t= %f\ttheta= %f\t w= %f\n", t0 + i * h, theta[i], w[i]);
//printf("\n\n");
//
//// Copy and paste the output in MATLAB and PLOT
//for (int i = 0; i < N; i++)
//	printf("%f\t%f\t%f\n", t0 + i * h, theta[i], w[i]);
//Matrix omega = createMat(N, 1);
//for (int i = 0; i < N; i++)
//{
//	omega.at[i][0] = w[i];
//}
//
//Matrix acceleration = gradient(time, omega);
//printMat(acceleration, "ans for b");
//// Problem 3 (b) Ends
//// Problem 3 (c) starts 
//Matrix thet = arr2Mat(theta, N, 1);
//
//double Q3_I = integral2(thet, t0, tf, N);
//
//
//printf("\nanswer of c) theta(t=0to2) = %f\n\n", Q3_I;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////



///* Problem 1*/
////void main(void)
////{
////
////
////   double V[] = { 0, 2.5, 6.2,8.2, 9.6, 12.5, 13.1, 14.1, 16.6, 18.3 };// [V]
////   double I[] ={0, 2, 5.0,7.0, 8.0, 10, 11,  12, 14 , 15};//[mA]
////   double n = sizeof(V) / sizeof(double);
////
////   Matrix V_m = arr2Mat(V, n, 1);
////   Matrix I_m = arr2Mat(I, n, 1);
////   Matrix res = linearFit(I_m , V_m); // kiloohm
////   printMat(res, "rest");
////   Matrix Vhat = createMat(n, 1);
////   /*Finding Residuals rk = yk - (a1xk-a0)*/
////   Matrix resi = createMat(n, 1);
////   double Error = 0;
////   initMat(resi, 0);
////   residual(I_m, V_m, resi, n, res);
////   TotalError(I_m, V_m, n, res, Error, resi);
////   printMat(resi, "residual");
////   initMat(Vhat, 0);
////   double len = 2 * n;
////   for (int i = 0; i < n; i++)
////   {
////      Vhat.at[i][0] = res.at[0][0]* I_m.at[i][0] + res.at[1][0];
////   }
////   printMat(Vhat, "result");
////
////
////   return;
////
////
////}
///* Problem 1 Ends */
//
///* Problem 2 */
////void main(void)
////{
////   double u[] = { 4.72,12.49,20.03,28.33,37.47,41.43,48.38,55.06,66.77,59.16,54.45,47.21,42.75,32.71,25.43,8.15 };
////   double v[] = { 7.18,7.3,7.37,7.42,7.47,7.5,7.53,7.55,7.58,7.56,7.55,7.53,7.51,7.47,7.44,7.28 };
////   /*u = Ae^BV = ln(u) = Bv + lnA*/
////   double n = sizeof(u) / sizeof(double);
////   Matrix Y = createMat(n, 1);
////   Matrix V = arr2Mat(v, n, 1);
////for (int i = 0; i < n; i++)
////{
////   Y.at[i][0] = log(u[i]);
////}
////Matrix coef = linearFit(V, Y);
////double B = coef.at[0][0];
////double A = exp(coef.at[1][0]);
////printMat(coef, "coefficient");
////Matrix uhat = createMat(n, 1);
////
////for (int i = 0; i < n ; i++)
////{
////   uhat.at[i][0] = A * exp(B * V.at[i][0]);
////}
////printMat(uhat, "res");
////}
//
///* Problem 3 */
////double diffforward_3(Matrix _x, Matrix _y);
////double diffbackward_3(Matrix _x, Matrix _y);
////void main(void)
////{
////   double x[] = { 0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0 };
////   double T[] = { 473.0,446.3,422.6,401.2,382.0,364.3,348.0,332.7,318.1,304.0,290.1 };
////   double len = sizeof(x) / sizeof(double);
////   Matrix xmat = arr2Mat(x, len, 1);
////   Matrix tmat = arr2Mat(T, len, 1);
////   for (int i = 0; i < len; i++)
////   {
////      xmat.at[i][0] /= 100;// m- transformation
////   }
////   double res1 = diffforward_3(xmat, tmat);
////   double res2 = diffbackward_3(xmat, tmat);
////   double k = 240;
////   res1 *= (-k);
////   res2 *= (-k);
////   printf("1 : %f 2 :%f", res1, res2);
////}
///* Proble 3 Ends */
//
///* Problem 4 */
////double w =10.0*pow(10.0,3);
////double L = 1.0;
////double E = 200.0 * pow(10.0, 9);
////double I = 2.1 * pow(10.0, -4);
////double a = 0.0;
////double b = 1.0;
////double N = 100.0;
////double h = (b - a) / 100.0;
////double myfunc(const double x);
////void main(void)
////{
////   
////   double interes = 0.0;
////   interes = integral(myfunc, a, b, N);
////   printf("%f", interes);
////
////}
////double myfunc(const double x)
////{
////   double num = pow((-w * (pow(x, 2) / 2)), 2);
////   double den = 2.0 * E * I;
////   
////   return num/den;
////}
///* Problem 4 Ends 02:43*/
//
///* Problem 5 02:43*/
////void main(void)
////{
////   
////}
//
///*ODE IVP Exercise*/
//#define PI 3.14159265368979323846264338327950288412
//
//// Problem1: Single equation of 1st order ODE
//double odeFunc_rc(const double t, const double v);
//
//// Problem2: Single equation of 2nd order ODE
//
//
//// Single Equation : odeEM, odeEU
//
//
///*-------------------------------------------------------------------------------\
//@ Numerical Methods by Young-Keun Kim - Handong Global University
//
//Author           : Jan Park, YKKIM
//Created          : 2021-06-03
//Modified         : 2021-06-03  by YKKIM
//Language/ver     : C++ in MSVS2017
//
//Description      : [Tutorial]ODE_IVP_student.c
//-------------------------------------------------------------------------------*/
//
//
////#include "../../include/myNM.h"
//#include <stdio.h>
//#include <math.h>
//
//#define ODE_EU 0
//#define ODE_EM 1
//#define ODE_RK2 2
//#define ODE_RK4 3
//
////  PI is defined in  myNM.h
//#define PI 3.14159265368979323846264338327950288412
//
//// Problem1: Single equation of 1st order ODE
//double odeFunc_rc(const double t, const double v);
//// Problem2: Single equation of 2nd order ODE
//void odeFunc_mck(const double t, const double Y[], double dYdt[]);
//// Single Equation : odeEM, odeEU
//void odeEU(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
//void odeEM(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
//void odeRK2(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
//void odeRK4(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
//// 2nd order Equations : sys2RK2, sys2RK4
//void sys2rk2(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);
//void sys2RK4(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);
//void odeFunc_something(const double t, const double Y[], double dYdt[]);
//
//




/////////////////////////////////////////////////////////Example Pendulum////////////////////////
//double a = 0.0;
//double b = 9.0;
//double h = 0.1;
//double N = (b - a) / h;
//double init_theta = PI / 2;
//double init_theta_dot = 0.2;
//double t[91] = { 0.0, };
//double y[91] = { 0.0, };
//double y_dot[91] = { 0.0, };
//
//sys2RK4(Pendulum, y, y_dot, a, b, h, init_theta, init_theta_dot, t);
//FileoutVec3_csv("Pendulum.csv", y, y_dot, t, N);
//printVec2(y, y_dot, N);