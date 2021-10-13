#pragma once
/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [YOUR NAME]
Created          : 26-03-2018
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H
#define     _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
# define PI                          double(3.141592653589793238462)
# define CROSSLINE                   printf("========================================================================\n")
# define LOOP(Count, End, a,b)       for(int Count = a; Count<End; Count=Count+b )
typedef unsigned int     Count;
//================================================= Differentiation===============================================
// Return the dy/dx results for the input data. (truncation error: O(h^2)
void gradient1D(double x[], double y[], double dydx[], int m);	
void gradientFunc(double func(const double x), double x[], double dydx[], int m);
//=================================================Integration===================================================
// Integration using rectangular method for discrete data inputs

double IntegrateRect(double _x[], double _y[], int _m);
double trapz(double _x[], double _y[], int _m);
double Simpson13(double myFunc(const double _x), double _a, double _b, int _m);
double Simpson38(double myFunc(const double x), double _a, double _b, int _m);

//============================================1st order ODE-VIP=============================================
void odeEU  (double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h);
void odeEM  (double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h);
void odeRK2 (double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h);
void odeRK3 (double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h);
void ode    (double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h,int method,int row);
void sys2RK2(void func(const double _t, const double _Y[], double _dYdt[]), double _y1[], double _y2[], double _t0, double _tf, double _h, double _y1_init, double _y2_init);
enum ode {EU=1,EM,RK2,RK3};



//void gradientFunc(double func(const double x), double x[], double dydx[], int m);
double myFunc(const double _x);
void printVec(double* _vec, int _row);
void FileoutVec1_csv(const char* _filename, double* _vec1, int _row);
void FileoutVec2_csv(const char* _filename, double* _vec1, double* _vec2, int _row);
#endif