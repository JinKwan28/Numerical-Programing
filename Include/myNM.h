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
# define D2R                         (PI/180)
# define CROSSLINE                   printf("========================================================================\n")
# define LOOP(Count, End, a,b)       for(int Count = a; Count<End; Count=Count+b )
# define MINUS                       double (-1)
typedef unsigned int     Count;
//==============================================Taylor series======================================================//
extern double  factorial(double _x);
extern double  sinTaylor(double _x);
extern double  sindTaylor(double _x);

//================================================= Differentiation===============================================
// Return the dy/dx results for the input data. (truncation error: O(h^2)
extern void gradient1D(double x[], double y[], double dydx[], int m);	
extern void gradientFunc(double func(const double x), double x[], double dydx[], int m);
extern double DifThreeFor(double _t[], double _X[], int _index);
extern double DifThreeBack(double _t[], double _X[], int _index);


//=================================================================================================================//
//=================================================Integration=====================================================//
// Integration using rectangular method for discrete data inputs

extern double IntegrateRect(double _x[], double _y[], int _m);
extern double Integralmid(double _x[], double _y[], int _m);
extern double trapz(double _x[], double _y[], int _m);
extern double trapz_func(double myFunc(const double _x), double _a, double _b, int _m);
extern double Simpson13(double _x[], double _a, double _b, int _m);
extern double Simpson13func(double myFunc(const double _x), double _a, double _b, int _m);
extern double Simpson38(double myFunc(const double x), double _a, double _b, int _m);


//=================================================================================================================//
//============================================1st order ODE-VIP=============================================
extern void odeEU  (double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h, double _t[]);
extern void odeEM  (double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h, double _t[]);
extern void odeRK2 (double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h, double _t[]);
extern void odeRK3 (double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h, double _t[]);
extern void odeRK4 (double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h, double _t[]);
extern void odePC  (double myFunc(const double _t, const double _y), double _y[], double _ypc[], double _t0, double _tf, double _h, double _y0,double _t[]);
extern void ode    (double myFunc(const double _t, const double _y), double _y[], double _t0, double _tf, double _h,int method,int row,double _t[]);
extern void sys2RK2(void func(const double _t, const double _Y[], double _dYdt[]), double _y1[], double _y2[], double _t0, double _tf, double _h, double _y1_init, double _y2_init,double _t[]);
extern void sys2RK4(void func(const double _t, const double _Y[], double _dYdt[]), double _y1[], double _y2[], double _t0, double _tf, double _h, double _y1_init, double _y2_init,double _t[]);
enum ode {EU=1,EM,RK2,RK3,RK4};


extern void printVec1(double* _vec, int _row);
extern void printVec2(double* _vec1, double* _vec2, int _row);
extern void FileoutVec1_csv(const char* _filename, double* _vec1, int _row);
extern void FileoutVec2_csv(const char* _filename, double* _vec1, double* _vec2, int _row);
extern void FileoutVec3_csv(const char* _filename, double* _vec1, double* _vec2, double*_vec3, int _row);
#endif