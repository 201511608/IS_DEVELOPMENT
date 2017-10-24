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
	// ultimate moment of resistance of tension steel  
	// [IS 456 ANNEXURE G]
	return 0.87*fst*Ast*d*(1-((Ast*fst)/(b*d*fck)));
}
double Moment_single_rect_EQsteel_in_Ast_P(double fst, double fck, double d, double b,double Ast_P)
{
	return 0.87*fst*Ast_P*(1 - ((1.005*fst*Ast_P) / fck))*b*d*d;  //Ast_P  0 to 1 do not multiply 100
}
double Moment_single_rect_EQconc( double fck,double xu, double Ast, double d, double b)
{
	double a = 0.36*(xu/d)*(1 - 0.42*(xu/d))*b*d*d*fck;
		return a;
}
double Moment_limit_single_rect(double xu_d_lim,double b, double d,double fck)
{
	// Concrete fail Limits due to exceeding in Strain
	return 0.36*xu_d_lim*(1-0.42*xu_d_lim)*b*d*d*fck;
}
double Es()
{
	return 200000;
}
double E_c(double fck)
{
	return 5000 * sqrt(fck);
}
double f_cr(double fck)
{
	return 0.75 * sqrt(fck);
}
double z_leverarm(double d, double x_by_d)
{
	// levearm from compression in concrete and tensile in steel
	return d*(1-0.416*x_by_d);
}
double e_s_max(double fst)
{
	return (0.87*fst) / Es() + 0.002;

}
double xu_by_d(double fst, double fck, double Ast, double d, double b)
{
	// xu_by_d from balanced C=T condition
	return (0.87*fst*Ast) / (0.36*fck*b*d);
}
double xu_by_d(double Mu, double fck, double d, double b)
{
	// x by d for given Mu
	return 1.2-sqrt((1.2)*(1.2) - (6.68*Mu) / (fck*b*d*d));
}
double Xu_by_d_limit(double ec, double es)
{
	return ec / (ec + es);
}
double Tensile_force_noyield(double fst, double Ast)  // not 0.87 * fst * Ast
{
	return fst*Ast;
}
double Tensile_force_yield(double fst, double Ast)  // not 0.87 * fst * Ast
{
	return 0.87*fst*Ast;
}
double Compression_force(double fck,double b,double x)
{
	return 0.36*fck*b*x;
}
double Ast_P(double Ast,double b,double d)
{
	return (Ast/(b*d)); // 0 to 1 NOT in terms of 100
}
double Ast_req_P(double fck, double fst, double xu_by_d)
{
	// for given Mu b d we get Ast
	// pno 50
	return 41.3*(fck / fst)*(xu_by_d); // 0 to 1 NOT in terms of 100
}
double Ast_req_P_limit(double fck, double fst, double xu_by_d_limit)
{
	// for given Mu b d we get Ast
	// pno 50
	return 41.3*(fck/fst)*(xu_by_d_limit); // 0 to 1 NOT in terms of 100
}
double Ast_req(double Mu, double fst, double z_leverarm)
{
	// Ast required from b d Mu and d in must be grater than min requre 
	// d > min depth d= Min_depth_from_Mu function
	return Mu/(0.87*fst*z_leverarm);
}
double Ast_req(double Mu, double fst, double b, double d,double fck)
{
	// Ast required from b d Mu and d in must be grater than min requre 
	// d > min depth d= Min_depth_from_Mu function
	// 45 pc varghis
	double q = 0.50-sqrt(0.25- (1.15*Mu / (fck*b*d*d)) );
	return q*b*d*fck/fst;
}
double Ast_max(double b, double d)
{
	// d > min depth d= Min_depth_from_Mu function
	// Clause [26.5.1]
	return 0.04*b*d;
}
double Ast_max_p()
{
	// Clause [26.5.1]
	return 0.04; // In terms of 0 to 1
}
double Ast_min(double fst, double b, double d)
{
	// d > min depth d= Min_depth_from_Mu function
	//[26.5.1]
	return (0.85/fst)*b*d; 
}
double Ast_min_p(double fst)
{
	// d > min depth d= Min_depth_from_Mu function
	//[26.5.1]
	return (0.85 / fst) ; //in terms of 0 to 1 
}
double Min_depth_from_Mu(double Xu_by_d_limit,double Mu,double fck,double b) {
	// Minimum depth from Mu
	//pno 52 pc vergish
// For any Input of Mu  min d  output wil be given 
	double K = 0.36*(Xu_by_d_limit*(1 - 0.42*Xu_by_d_limit));
	return  sqrt((Mu) / (K*fck*b)); 
}
double es_SteelStrain_compression(double ec_concreteStrain, double nutral_axis, double effectivecover)
{
	// Due to compression in steel
	return (ec_concreteStrain * (nutral_axis-effectivecover)) / nutral_axis;
}
double es_SteelStrain_tensile(double ec_concreteStrain, double d, double nutral_axis)
{
	// Due to tension in steel
	return (ec_concreteStrain * (d - nutral_axis)) / nutral_axis;
}
double ec_ConcreteStrain(double es_SteelStrain_tensile, double d, double nutral_axis)
{
	// Due to compression in Concrete
	return (es_SteelStrain_tensile * nutral_axis )/(nutral_axis - nutral_axis);
}




//Print
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
				Moment_single_rect_EQconc(0,nutral_axis_single_rect(fst, Ast, fck, b), Ast,  d, b),
				Moment_limit_single_rect(Xu_by_d_limit(0.0035,e_s_max(415)),b, d,fck),
				(Moment_single_rect_EQsteel(fst,fck, Ast, d,b) < Moment_limit_single_rect(Xu_by_d_limit(0.0035,e_s_max(415)),b, d,fck))*1.0 });
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
	     	<< " Mu<MuLimt_check [0 concrete failed 1 steel yield ]:" << matrix2d[i][6];
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
	vector<vector<double>> matrix2d;


	

	BEAM beam1;
	double fst = 415, Ast = 1000, fck = 20, b = 300;
	beam1.nutral_axis_single_rect(fst, Ast, fck, b); 
	beam1.Prdouble_values(matrix2d);


	getchar();
    return 0;
}

