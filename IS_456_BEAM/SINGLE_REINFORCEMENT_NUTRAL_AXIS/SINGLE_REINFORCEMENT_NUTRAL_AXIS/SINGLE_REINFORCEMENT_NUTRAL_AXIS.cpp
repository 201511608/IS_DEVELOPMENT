// SINGLE_REINFORCEMENT_NUTRAL_AXIS.cpp : Defines the entry podouble for the console application.
//

#include "stdafx.h"
#include<stdio.h>
#include <iostream>
#include <vector>

using namespace std;

class BEAM
{

public:
double nutral_axis_single_rect(double fst, double Ast, double fck, double b) 
{
		return (0.87*fst*Ast) / (0.36*fck*b);
	}

double Moment_single_rect_EQsteel(double fst,double fck, double Ast,double d , double b)
{
	return 0.87*fst*Ast*d*(1-((Ast*fst)/(b*d*fck)));
}
double Moment_single_rect_EQconc( double fck,double xu, double Ast, double d, double b)
{
	double a = 0.36*(xu/d)*(1 - 0.42*(xu/d))*b*d*d*fck;
		return a;
}
double Moment_limit_single_rect(double xu_d_lim,double b, double d,double fck)
{
	return 0.36*xu_d_lim*(1-0.42*xu_d_lim)*b*d*d*fck;
}
double Es()
{
	return 200000;
}
double e_s_max(double fst)
{
	return (0.87*fst) / Es() + 0.002;

}
double Xu_d_limit(double ec, double es)
{
	return ec / (ec + es);
}
void Prdouble_values(vector<vector<double>> matrix2d)
{
	double d = 400;
	double fck = 20;
	double fst = 415;
	for (double b = 100; b <= 500; b = b + 100)
	{

		for (double Ast = 100; Ast <= 1000; Ast = Ast + 100)
		{
			matrix2d.push_back({ Ast,
								   b,
				(nutral_axis_single_rect(fst, Ast, fck, b)),
				Moment_single_rect_EQsteel(fst,fck, Ast, d,b),
				Moment_single_rect_EQconc(fck,nutral_axis_single_rect(fst, Ast, fck, b), Ast,  d, b),
				Moment_limit_single_rect(Xu_d_limit(0.0035,e_s_max(415)),b, d,fck),
				(Moment_single_rect_EQsteel(fst,fck, Ast, d,b) < Moment_limit_single_rect(Xu_d_limit(0.0035,e_s_max(415)),b, d,fck))*1.0 });
		}
	}

	for (double i = 0; i <= matrix2d.size() - 1; i++)
	{
		cout << endl << " Breadth:" << matrix2d[i][1]
			<< " Ast: " << matrix2d[i][0]
			<< " Nutral Axis:" << matrix2d[i][2]
			<< " Moment_from_steelEQ:" << matrix2d[i][3]
			<< " Moment_from_concEQ:" << matrix2d[i][4]
			<< " Limiting_Moment:" << matrix2d[i][5]
	     	<< " Mu<MuLimt_check :" << matrix2d[i][6];
	}
}

};



int main()
{
	vector<vector<double>> v2d;


	vector<double>v1;
	vector<double>::iterator it;
	it = v1.begin();
	cout << &it ;
	double k = 0;

	vector<vector<vector<double>>> matrix3d;// (3, vector<vector<double>>(6, vector<double>(10, 0)));
	//vector<vector<double>> matirx2d(100*100, vector<double>(3, 0));
	vector<vector<double>> matrix2d;


	

	BEAM beam1;
	double fst = 415, Ast = 1000, fck = 20, b = 300;
	beam1.nutral_axis_single_rect(fst, Ast, fck, b); 
	beam1.Prdouble_values(matrix2d);


	getchar();
    return 0;
}

