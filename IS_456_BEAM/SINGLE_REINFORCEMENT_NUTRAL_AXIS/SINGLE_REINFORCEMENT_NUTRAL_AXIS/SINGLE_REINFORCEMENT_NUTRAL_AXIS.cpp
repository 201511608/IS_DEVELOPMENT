// SINGLE_REINFORCEMENT_NUTRAL_AXIS.cpp : Defines the entry podouble for the console application.
//

#include "stdafx.h"
#include<stdio.h>
#include <iostream>
#include <vector>

using namespace std;

class BEAM
{
	// Dimension   [ N  MM ]
public:
	double nutral_axis_single_rect(double fst, double Ast, double fck, double b)
	{
		return (0.87*fst*Ast) / (0.36*fck*b);
	}
	double Moment_single_rect_EQsteel(double fst, double fck, double Ast, double d, double b)
	{
		// ultimate moment of resistance of tension steel  
		// [IS 456 ANNEXURE G]
		return 0.87*fst*Ast*d*(1 - ((Ast*fst) / (b*d*fck)));
	}
	double Moment_single_rect_EQsteel_in_Ast_P(double fst, double fck, double d, double b, double Ast_P)
	{
		return 0.87*fst*Ast_P*(1 - ((1.005*fst*Ast_P) / fck))*b*d*d;  //Ast_P  0 to 1 do not multiply 100
	}
	double Moment_single_rect_EQconc(double fck, double xu, double Ast, double d, double b)
	{
		double a = 0.36*(xu / d)*(1 - 0.42*(xu / d))*b*d*d*fck;
		return a;
	}
	double Moment_limit_single_rect(double xu_d_lim, double b, double d, double fck)
	{
		// Concrete fail Limits due to exceeding in Strain
		return 0.36*xu_d_lim*(1 - 0.42*xu_d_lim)*b*d*d*fck;
	}
	double Center_of_Concrete_Compression(double xu)
	{
		// Centroid of Concrete compression stress  // distance from top fibre to centoidal concrete compressive stress 
		return 0.42*xu;
	}
	double Effective_cover(double clearcover, double dia)
	{
		// Define top and botton cover differently by calling this function !! with different diameter
		return clearcover + (dia / 2);
	}
	double modular_ratio_m(double E_s, double E_c)
	{
		return (E_s / E_c);
	}
	double E_s()
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
		return d*(1 - 0.416*x_by_d);
	}
	double e_s_max(double fst)
	{
		return (0.87*fst) / E_s() + 0.002;

	}
	double xu_by_d(double fst, double fck, double Ast, double d, double b)
	{
		// xu_by_d from balanced C=T condition
		return (0.87*fst*Ast) / (0.36*fck*b*d);
	}
	double xu_by_d(double Mu, double fck, double d, double b)
	{
		// x by d for given Mu
		return 1.2 - sqrt((1.2)*(1.2) - (6.68*Mu) / (fck*b*d*d));
	}
	double Xu_by_d_limit(double ec, double es)
	{
		// es =e_s_max() ??
		return ec / (ec + es);
	}
	double Tensile_force_Steel_noyield(double fst, double Ast)  // not 0.87 * fst * Ast
	{
		return fst*Ast;
	}
	double Tensile_force_steel_yield(double fst, double Ast)  // not 0.87 * fst * Ast
	{
		return 0.87*fst*Ast;
	}
	double Compression_force_concrete(double fck, double b, double x)
	{
		return 0.36*fck*b*x;
	}
	double Ast_P(double Ast, double b, double d)
	{
		return (Ast / (b*d)); // 0 to 1 NOT in terms of 100
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
		return 41.3*(fck / fst)*(xu_by_d_limit); // 0 to 1 NOT in terms of 100
	}
	double Ast_req_P_Balanced(double fck, double fst, double xu_by_d_limit)
	{
		// for given Mu b d we get Ast
		// pno 50
		return 41.3*(fck / fst)*(xu_by_d_limit); // 0 to 1 NOT in terms of 100
	}
	double Ast_req(double Mu, double fst, double z_leverarm)
	{
		// Ast required from b d Mu and d in must be grater than min requre 
		// d > min depth d= Min_depth_from_Mu function
		return Mu / (0.87*fst*z_leverarm);
	}
	double Ast_req(double Mu, double fst, double b, double d, double fck)
	{
		// Ast required from b d Mu and d in must be grater than min requre 
		// d > min depth d= Min_depth_from_Mu function
		// 45 pc varghis
		double q = 0.50 - sqrt(0.25 - (1.15*Mu / (fck*b*d*d)));
		return q*b*d*fck / fst;
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
		return (0.85 / fst)*b*d;
	}
	double Ast_min_p(double fst)
	{
		// d > min depth d= Min_depth_from_Mu function
		//[26.5.1]
		return (0.85 / fst); //in terms of 0 to 1 
	}
	double Min_depth_from_Mu(double Xu_by_d_limit, double Mu, double fck, double b) {
		// Minimum depth from Mu
		//pno 52 pc vergish
	// For any Input of Mu  min d  output wil be given 
		double K = 0.36*(Xu_by_d_limit*(1 - 0.42*Xu_by_d_limit));
		return  sqrt((Mu) / (K*fck*b));
	}
	double es_SteelStrain_compression(double ec_concreteStrain, double nutral_axis, double effectivecover)
	{
		// Due to compression in steel
		// effectivecover=clearcover + effectivecover + d/2 ;
		return (ec_concreteStrain * (nutral_axis - effectivecover)) / nutral_axis;
	}
	double es_SteelStrain_tensile(double ec_concreteStrain, double d, double nutral_axis)
	{
		// Due to tension in steel
		return (ec_concreteStrain * (d - nutral_axis)) / nutral_axis;
	}
	double ec_ConcreteStrain(double es_SteelStrain_tensile, double d, double nutral_axis)
	{
		// Due to compression in Concrete
		return (es_SteelStrain_tensile * nutral_axis) / (nutral_axis - nutral_axis);
	}
	//
	double Ast_min__compressionSteel(double b, double xu)
	{
		// wrt to compression area
		return 0.04 * b * xu;
	}
	double Ast_min__compressionSteel_2(double b, double d)
	{
		// wrt whole area
		return 0.02 * b * d;
	}
	double Ast_max__compressionSteel(double b, double D)
	{
		// [26.5.1.2]
		// D ? total 
		return 0.04 * b * D;
	}
	//
	double Interpolation_funtion(double x, double x1, double y1, double x2, double y2)
	{
		return (((y2 - y1) / (x2 - x1))*(x - x1)) + y1;
	}
	double Compessive_Stress_steel(double es_SteelStrain_compression, double fst)
	{
		// else stress 0 if exceeds limits fails
		// SP 16 pno 4 5
		if (fst == 250)
		{
			if ((es_SteelStrain_compression > 0) && (es_SteelStrain_compression <= 0.0001))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0, 0, 0.0001, 212.5);
			}
			else if ((es_SteelStrain_compression > 0.0001) && (es_SteelStrain_compression <= 0.0003))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0.0001, 212.5, 0.0003, 225);
			}
			else if ((es_SteelStrain_compression > 0.0003) && (es_SteelStrain_compression <= 0.0007))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0.0003, 225, 0.0007, 237.5);
			}
			else if ((es_SteelStrain_compression > 0.0007) && (es_SteelStrain_compression <= 0.0010))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0.0007, 237.5, 0.001, 243.75);
			}
			else if ((es_SteelStrain_compression > 0.001) && (es_SteelStrain_compression <= 0.002))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0.001, 243.75, 0.002, 250);
			}
			else return 0;

		}
		else if (fst == 415)
		{
			if ((es_SteelStrain_compression > 0) && (es_SteelStrain_compression <= 0.00144))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0, 0, 0.00144, 288);
			}
			else if ((es_SteelStrain_compression > 0.00144) && (es_SteelStrain_compression <= 0.00163))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0.00144, 288, 0.00163, 306);
			}
			else if ((es_SteelStrain_compression > 0.00163) && (es_SteelStrain_compression <= 0.00192))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0.00163, 306, 0.00192, 324);
			}
			else if ((es_SteelStrain_compression > 0.00192) && (es_SteelStrain_compression <= 0.00241))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0.00192, 324, 0.00241, 342);
			}
			else if ((es_SteelStrain_compression > 0.00241) && (es_SteelStrain_compression <= 0.00276))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0.00241, 342, 0.00276, 351);
			}

			else if ((es_SteelStrain_compression > 0.00276) && (es_SteelStrain_compression <= 0.00386))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0.00276, 351, 0.00380, 360);
			}
			else return 0;

		}
		else if (fst == 500)
		{
			if ((es_SteelStrain_compression > 0) && (es_SteelStrain_compression <= 0.00174))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0, 0, 0.00174, 347);
			}
			else if ((es_SteelStrain_compression > 0.00174) && (es_SteelStrain_compression <= 0.00195))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0.00174, 347, 0.00195, 369);
			}
			else if ((es_SteelStrain_compression > 0.00195) && (es_SteelStrain_compression <= 0.00226))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0.00195, 369, 0.00226, 391);
			}
			else if ((es_SteelStrain_compression > 0.00226) && (es_SteelStrain_compression <= 0.00277))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0.00226, 391, 0.00277, 413);
			}
			else if ((es_SteelStrain_compression > 0.00277) && (es_SteelStrain_compression <= 0.00312))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0.00277, 413, 0.00312, 423);
			}
			else if ((es_SteelStrain_compression > 0.00312) && (es_SteelStrain_compression <= 0.00417))
			{
				return Interpolation_funtion(es_SteelStrain_compression, 0.00312, 423, 0.00417, 434);
			}
			else return 0;

		}
	}
	double compression_force_steel(double Compessive_Stress_steel, double Ast_req_compression)
	{
		return Compessive_Stress_steel*Ast_req_compression;
	}

	double Moment_Mu2_double_rect_Tensile_only_EQsteel(double fst, double Ast2_req_double_reinforcment, double D, double Compessive_Stress_steel, double Top_effective_cover, double Bot_effective_cover)
	{
		// 0.87*fst*Ast2*(d-d');
		return 0.87 * fst*Ast2_req_double_reinforcment*(D - (Top_effective_cover - Bot_effective_cover));
	}
	double Moment_Mu2_double_rect_Comp_only_EQsteel(double Ast_req_compression, double D, double Compessive_Stress_steel, double Top_effective_cover, double Bot_effective_cover)
	{
		// D = total depth incleding effective covers
		//Effective_cover= clearcover + dia/2
		return Compessive_Stress_steel*Ast_req_compression*(D - (Top_effective_cover - Bot_effective_cover));
	}
	double Moment_Mu2_double_rect_(double Total_given_Mu_Moment,double Moment_Mu1_double_rect_)
	{
		return Total_given_Mu_Moment - Moment_Mu1_double_rect_;
	}
	double Moment_Mu1_double_rect_Comp_only_EQconc(double xu_d_lim, double b, double d, double fck)
	{
		// d= total depth - effective cover
		return Moment_limit_single_rect(xu_d_lim, b, d, fck);
	}
	double Moment_Mu1_double_rect_Tens_only_EQsteel(double fst, double fck, double Ast1, double d, double b)
	{
		// d= total depth - effective cover
		// Ast1= Ast limit ..
		// Cross check Mu1 comp = Mu1 tension
		return Moment_single_rect_EQsteel(fst, fck, Ast1, d, b);
	}
	double Ast1_req_double_reinforcment(double fck, double fst, double xu_by_d_limit)
	{
		return	Ast_req_P_Balanced(fck, fst, xu_by_d_limit);
	}
	double Ast2_req_double_reinforcment(double Total_given_Mu_Moment, double Moment_Mu1_double_rect_Comp_only_EQconc, double fst, double total_depth, double Top_effective_cover, double Bot_effective_cover)
	{
		// Total_given_Mu_Moment from USER
		double Mu2 = Total_given_Mu_Moment - Moment_Mu1_double_rect_Comp_only_EQconc;
		return	Mu2 / (0.87*fst*(total_depth - (Top_effective_cover - Bot_effective_cover)));
	}
	double Ast2_req_double_reinforcment(double Ast_total_tensile_double, double Ast1_req_double_reinforcment)
	{
		// For Know total Ast tensile steel in tension
		return	Ast_total_tensile_double - Ast1_req_double_reinforcment;
	}
	double Astc_req_double_reinforcment(double Compessive_Stress_steel, double fst, double Ast2)
	{
		// Acs area of steel in compression
		return  (Ast2*0.87*fst) / Compessive_Stress_steel;
	}
	double Moment_Muc_double_rect_total_compression(double Moment_Mu1_double_rect_Comp_only_EQconc, double  Moment_Mu2_double_rect_Comp_only_EQsteel)
	{
		// comperssion failure in        concrete comp + steel comp
			//Moment_Mu1_double_rect_Comp_only_EQcon = Limit
		return Moment_Mu1_double_rect_Comp_only_EQconc + Moment_Mu2_double_rect_Comp_only_EQsteel;
	}
	double Moment_Mut_double_rect_total_tension(double Moment_Mu1_double_rect_Tens_only_EQsteel, double   Moment_Mu2_double_rect_Tensile_only_EQsteel)
	{
		// comperssion failure in        concrete comp + steel comp
		//  Moment_Mu1_double_rect_Tens_only_EQsteel == Moment_Mu1_double_rect_Comp_only_EQcon
		return Moment_Mu1_double_rect_Tens_only_EQsteel + Moment_Mu2_double_rect_Tensile_only_EQsteel;
	}


	//
	double Service_Stress_fs(double Beta, double fst, double Ast_req, double Ast_provided)
	{
		//	Beta = (moment_after_redistribution)/(moment_Before redistribution);
		Beta = 1; // Assuming Beta as 1 if not assume as above;
		return (0.58 * fst* Ast_req / (Ast_provided))*Beta;
	}
	double Service_Stress_fs(double fst)
	{
		//	Beta = (moment_after_redistribution)/(moment_Before redistribution);
		return (0.58 * fst);
	}
	double Modification_Factor_F1(double Service_Stress_fs, double Pt)
	{
		// pt =100*Ast/(b*d)
		// Modification factor for tension Reinforcement
		return	1 / (0.225 + 0.00322*Service_Stress_fs + 0.625*log10(Pt));  // Must be <= 2
	}
	double Modification_Factor_F2(double Ast_compression)
	{
		// pc =100*Ast_compression/(b*d)
		// Modification factor for Compression Reinforcement
		return	(1.6*Ast_compression) / (Ast_compression + 0.275); //Must be <= 1.5
	}
	double Modification_Factor_F3(double bw, double bf)
	{

		// Modification factor for T L beam F3
		if ((bw / bf) <= 0.3)
		{
			return 0.8;
		}
		else if (((bw / bf) > 0.3) && ((bw / bf) <= 0.8))
		{
			return 0.8 + (2 / 7)*((bw / bf) - 0.3);
		}
		else
			// bw/bf =1  then F3=1.0 ;
			return 1;
		{

		}
	}
	double Span_by_Depth_L_by_d(double span, double depth, double Modification_Factor_F1, double Modification_Factor_F2, double Modification_Factor_F3)
	{
		return (span / depth)*Modification_Factor_F1*Modification_Factor_F2*Modification_Factor_F3;
	}
	double Span_by_Depth_L_by_d_code(double x, double span)
	{
		// 1 Cantilever   2 Simply supported  3 Contineous // 
		if (x == 1)  // cantilever
			if (span <= 10)
				return 7;
			else
				return 7 * 10 / span;
		if (x == 2)  // Simply supported
			if (span <= 10)
				return 20;
			else
				return 20 * 10 / span;
		if (x == 3)  // Contineous
			if (span <= 10)
				return 26;
			else
				return 26 * 10 / span;

	}
	//
	double flexural_strength_fcr(double fck)
	{
		// fck = Characterstick comressive strength !
		return 0.7*sqrt(fck);
	}
	double cracking_Moment_Mcr(double flexural_strength_fcr, double b, double d)
	{
		//Yt = d/2 generally Ig= b*d^3/12
		return flexural_strength_fcr*b*d*d / 6;
	}
	double cracking_Moment_Mcr(double flexural_strength_fcr, double Ig, int Yt)
	{
		//Yt = d/2 generally Ig= b*d^3/12
		return flexural_strength_fcr*Ig / Yt;
	}
	double cracking_Moment_Mcr_dsh(double cracking_Moment_Mcr)
	{
		//Yt = d/2 generally Ig= b*d^3/12
		// Page no 379 devdas
		// For safety
		return 0.7 * cracking_Moment_Mcr;
	}
	//
	double creep_cofficient_theta(double Age)
	{
		//[6.2.5.1]
		if (Age == 7)
			return 2.2;
		else if (Age == 28)
			return 1.6;
		else
			return 1.1;
	}
	double Moment_of_inertia(double b, double d)
	{
		// Igr = gross moment of inertia
		return (b*pow(d, 3)) / 12;
	}
	double Moment_of_inertia_Igr(double b, double d)
	{
		// Igr = gross moment of inertia
		return (b*pow(d, 3)) / 12;
	}
	double Moment_of_inertia_Icr_single(double modular_ratio_m,double nutral_axis_k_single_or_double ,double Ast,double b,double d)
	{
		double m = modular_ratio_m;
		double k = nutral_axis_k_single_or_double;
		return ((b*pow(k*d, 3)) / 3 + m*Ast*(pow((d - (k*d)), 2)));
	}
	double Moment_of_inertia_Ieff(double  Moment_of_inertia_Icr_single_or_double,double cracking_Moment_Mcr,double Moment_input,double b,double bw,double nutral_axis_k_single_or_double)
	{
		// This is applicable only in the M is grater than Mcr !!!
		// [C -2-14]
	double Icr = Moment_of_inertia_Icr_single_or_double;
	double	Mcr = cracking_Moment_Mcr;
	double M = Moment_input;
	double k = nutral_axis_k_single_or_double;
	double n = (1-(k/3))*(1-k)*(bw/b); // j =1-k/3  //381 Devdas 380
	return Icr/(1.2 - ((Mcr/M)*n))  ;  //381 Devdas
	}

	//
	double deflection_simply_sup(double w, double l, double E, double I)
	{
		// In terms of weight or florce applied
		return (5 * w*pow(l, 4)) / (384 * E*I);
	}
	double deflection_simply_sup_M(double M, double l, double E, double I)
	{
		return (5 * M * pow(l,2)) / (48 * E*I);
	}
	double deflection_ShortTerm(double load_w, double l, double E, double Ieff)
	{
		// In terms of weight or florce applied
		//M <= Mcr+   // Ieff = Igr 
		// //M > Mcr   Calculate  M 
		// deflection due to dead load as well as live load  D L D+L
		double w = load_w;
		return (5 * w*pow(l, 4)) / (384 * E*Ieff);
	}
	//
	double nutral_axis_k_single(double Ast,double b,double d,double modular_ratio_m)
	{
		// x=kd;  // Nutral axis value k___ WSM   
		// Single reinforcement P no 115 Devdas 
		double rho = Ast/(b*d);
		double m = modular_ratio_m;
		return sqrt((2 * rho*m) + pow((rho*m),2)) - (rho*m)  ;
	}
	//
	double Equivalent_ShearForce_Ve(double Vu, double Tu, double b) 
	{
		// Vu shear
		// Tu Torsional moment
		// [41.3.1]
		return (Vu + 1.6 * (Tu/b));
	}
	double Equivalent_ShearStress_Te(double Equivalent_ShearForce_Ve, double d, double b)
	{

		// [41.3.1]
		return (Equivalent_ShearForce_Ve / (b*d));
	}
	double Equivalent_Bending_Moment_Me(double Mu,double Tu,double D,double b)
	{
		// Bending Moment at cross Section
		// [41.4.2]
		// Tu Torsional moment
		//Pno 400 Pc varghis
		double Mt = Tu*((1 + (D / b)) / 1.70);
		if (Mu >= Mt)
			return Mu + Mt;
		else
			// MUST PROVEDE COMPRESSION STEEL FOR REMAINGIN MU PNO 409 PC VARGHESE
			return Mu - Mt;

	}
	double Asv_by_sv_Equv_Torsion_1(double Tu,double Vu, double fst,double b,double d,double effective_cover)
	{
		// Tu torsion moment from analysis NOt Tve
		// [41.4.3]
		double b1 = b - 2 * effective_cover;
		double d1 = d - effective_cover; // if D then d1 = D - 2*effective_cover  !!!
		return (Tu / (b1*d1*0.87*fst)) + (Vu / (2.5*d1*0.87*fst));
		// call Shear_spacing funcition for spacing
	}
	double Asv_by_sv_Equv_Torsion_2(double Equivalent_ShearStress_Ve, double Tc, double fst, double b, double sv, double effective_cover)
	{
		// Tu torsion moment
		// 
		// [41.4.3]
		// Tc ShearStress 
		double Tve = Equivalent_ShearStress_Ve;
		return ((Tve - Tc)*b)/(0.87*fst);
		// call Shear_spacing funcition for spacing
	}
	double Shear_spacing(double Asv_by_sv,double n_legs,double d)
	{
		// Asv_by_sv = Asv_by_sv_Equv_Torsion_1 or Asv_by_sv_Equv_Torsion_2
		return  (n_legs*3.143*pow(d, 2))/ (Asv_by_sv * 4) ;
	}
	
	//
	double Shear_stess_nominal_Tv(double Vu, double b, double d)
	{
		//[40.1]
		// For T L beam b is rib width not flang width
		return Vu / (b*d);
	}
	double Shear_force(double Shear_stess_nominal_Tv,double b,double d)
	{
		return Shear_stess_nominal_Tv*b*d;
	}
	double shear_stress_steel_Ts(double Shear_stress_strength_concrete_Tc,double Shear_stess_nominal_Tv)
	{
		//Tc <= Tc_man
		return Shear_stess_nominal_Tv - Shear_stress_strength_concrete_Tc ;
	}
	double shear_force_steel_Vs(double shear_stress_steel_Ts, double b,double d)
	{
		//Tc <= Tc_man
		return shear_stress_steel_Ts*b*d;
	}
	double Shear_stress_strength_concrete_Tc(double fck, double fst, double Ast_P)
	{
		//max shere can take by concrete
		//Shear_stress_concrete_Tc = Shear_strength_concrete_Tc
		//[is 456 Table 19]
		//[sp 24 -Section 39.2]
		double Pt = Ast_P * 100;   // pt of rib width only even in L T
		double beta = (0.8*fck) / (6.89*Pt); // beta must be less than 1 !!??!!
		return (0.85*sqrt(0.8*fck)*sqrt(1 + 5 * beta - 1)) / (6 * beta);
	}
	double Shear_stress_strength_concrete_Tc(double fck,double fst,double Ast,double b,double d)
	{
		//max shere can take by concrete
		//Shear_stress_concrete_Tc = Shear_strength_concrete_Tc
		//[is 456 Table 19]
		//[sp 24 -Section 39.2]
		double Pt=Ast/(b*d) * 100;   // pt of rib width only even in L T
		double beta = (0.8*fck)/(6.89*Pt); // beta must be less than 1 !!??!!
		return (0.85*sqrt(0.8*fck)*sqrt(1 + 5 * beta - 1)) / (6 * beta);
	}
	double Shear_stress_allowable_concrete_Tc_max(double fck)
	{
		//	[is 456 Table 20]
		//	[sp 24 pno 127	]
		return 0.62*sqrt(fck);
	}
	double Asv_by_sv_required(double b, double d,double fst,double Shear_stess_nominal_Tv, double Shear_stress_strength_concrete_Tc)
	{
		double Tv = Shear_stess_nominal_Tv	;
		double Tc = Shear_stress_strength_concrete_Tc;
		return (b*d*(Tv-Tc)) / (0.87*fst*d)	;		
		// call Shear_spacing funcition for spacing
	}
	double Asv_by_sv_Min_nominal_Shear_reinforcement(double b, double fst)
	{
		//[26.5.1.6]
		// Asv/bSv  >= 0.4/0.87fy
		return (0.4*b) / (0.4 / (0.87*fst));
		// Use Shear_spacing function to get spacing !
	}
	double shear_steel_Max_allowable_spacing(double d,double sv_req)
	{
		//	[26.5.5.1]
		//  sv_req  =[Shear_spacing] of [Asv_by_sv_required]
		double s1 = 300; 
		double s2 = 0.75*d;
		if ((sv_req <= s1) && (sv_req <= s2))
		{
			return sv_req;
		}
		else if ((s2 <= sv_req) && (s2 <= s1))
		{
			return s2;
		}
		else
			return s1;
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
					(Moment_single_rect_EQsteel(fst,fck, Ast, d,b) < Moment_limit_single_rect(Xu_by_d_limit(0.0035,e_s_max(415)),b, d,fck))*1.0 ,
					deflection_simply_sup_M(21.63*1000000, 4160, 25000  , 6.6667*100000000) });
			}

			for (double i = 0; i <= matrix2d.size() - 1; i++)
			{
				cout << endl << " Breadth:" << matrix2d[i][1]
					<< " Ast: " << matrix2d[i][0] * 1.0
					<< " Nutral Axis:" << matrix2d[i][2]
					<< " Moment_from_steelEQ:" << matrix2d[i][3]
					<< " Moment_from_concEQ:" << matrix2d[i][4]
					<< " Limiting_Moment:" << matrix2d[i][5]
					<< " Mu<MuLimt_check [0 concrete failed 1 steel yield ]:" << matrix2d[i][6]
					<< " Deflection:  " << matrix2d[i][7];
			}
		}

	};
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
	cout<<"nutral axis Workstress "<<beam1.nutral_axis_single_rect(250, 2945, 15, 600);
	getchar();
	double fst = 415, Ast = 1000, fck = 20, b = 300;
	beam1.nutral_axis_single_rect(fst, Ast, fck, b); 
	beam1.Prdouble_values(matrix2d);


	getchar();
    return 0;
}

