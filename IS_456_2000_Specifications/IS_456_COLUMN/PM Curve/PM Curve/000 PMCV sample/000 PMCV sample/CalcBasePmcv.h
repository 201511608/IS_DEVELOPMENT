// CalcBasePmcv.h: interface for the CCalcBasePmcv class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CALCBASEPMCV_H__825CF04D_530E_40F4_B15F_BCDD24E58BCB__INCLUDED_)
#define AFX_CALCBASEPMCV_H__825CF04D_530E_40F4_B15F_BCDD24E58BCB__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <afxtempl.h>

struct T_MDGN_PM_UNIT
{
	double dPhi; // strength reduction factor.
	double dXn;  // distance ratio from tensile edge to neutral axis.
	double dPn;  // nominal axial strength.
	double dMny; // nominal flexural strength about strong/major axis.
	double dMnz; // nominal flexural strength about weak/minor axis.
	double desi; // strain.	

	void Init() {}
	T_MDGN_PM_UNIT() { Init(); }
	T_MDGN_PM_UNIT& operator= (const T_MDGN_PM_UNIT &src) { return *this; }
};

struct T_MDGN_PM_CURV
{
	double dAngleNA;  // 중립축 회전각(Radian)
	double dDmax;     // 회전된 중립축에 대한 단면의 높이.
	double dDeff;
	double dPhiPnmax;
	T_MDGN_PM_UNIT eminPM;
	CArray<T_MDGN_PM_UNIT, T_MDGN_PM_UNIT&> aPmUnit;

	void Init() {}
	T_MDGN_PM_CURV() { Init(); }
	T_MDGN_PM_CURV& operator= (const T_MDGN_PM_CURV &src) { return *this; }
};

struct T_MDGN_PM_RBAR
{
	double dy;    // y coordinate.
	double dz;    // z coordinate.
	double dh;    // distance from neutral axis.
	double dAs;   // rebar area.
	double desi;  // strain.
	double dfsi;  // stress.	
	double dFsi;  // force.

	void Init() {}
	T_MDGN_PM_RBAR() { Init(); }
	T_MDGN_PM_RBAR& operator= (const T_MDGN_PM_RBAR &src) { return *this; }
};

struct T_MDGN_PM_SECT
{
	double dAngleNA;
	double dCb;  // distance ratio from compression edge to neutral axis.
	double dXn;  // distance ratio from tensile edge to neutral axis.
	double dPn;
	double dMny;
	double dMnz;
	double desiMax;
	double dPncc;
	double dPnsc;
	double dPnst;
	double dMnycc;
	double dMnysc;
	double dMnyst;
	double dMnzcc;
	double dMnzsc;
	double dMnzst;
	CArray<T_MDGN_PM_RBAR, T_MDGN_PM_RBAR&> arRbarRes;
	
	void Init() {}
	T_MDGN_PM_SECT() { Init(); }
	T_MDGN_PM_SECT& operator= (const T_MDGN_PM_SECT &src) { return *this; }	
	double GetTension()     { return dPnst; }
	double GetCompression() { return dPncc + dPnsc; }	
};

struct T_MDGN_CONC_POS_UNIT
{
	double dyc;
	double dzc;
	double dy[4];
	double dz[4];
	double dAc;

	void Init() {}
	T_MDGN_CONC_POS_UNIT() { Init(); }
	T_MDGN_CONC_POS_UNIT& operator= (const T_MDGN_CONC_POS_UNIT &src) { return *this; }
};

struct T_MDGN_RBAR_POS_UNIT
{
	double dy;
	double dz;
	double dAs;
	double dDia;

	void Init() {}
	T_MDGN_RBAR_POS_UNIT() { Init(); }
	T_MDGN_RBAR_POS_UNIT& operator= (const T_MDGN_RBAR_POS_UNIT &src) { return *this; }
};

struct T_MDGN_PM_RES // per load combination.
{
	double dPhi     ;
	double dPhiPn   ;
	double dPhiMn   ;
	double dPhiMny  ;
	double dPhiMnz  ;
	double dPhiPnmax;
	double dRatioP  ;
	double dRatioM  ;

	T_MDGN_PM_CURV PmcvData;	
	T_MDGN_PM_SECT CbPmData;
  T_MDGN_PM_SECT PuPnData;

	void Init();
	T_MDGN_PM_RES() { Init(); }
	T_MDGN_PM_RES& operator= (const T_MDGN_PM_RES &src) { return *this; }
};

struct T_MDGN_PM_CURV_3D
{
	int nDivNum90Deg;
	int nSectSymmType;
	CArray<T_MDGN_PM_CURV, T_MDGN_PM_CURV&> PmcvData;
	CArray<T_MDGN_PM_SECT, T_MDGN_PM_SECT&> CbPmData;
	CArray<T_MDGN_PM_SECT, T_MDGN_PM_SECT&> PuPnData;

	void Init();
	T_MDGN_PM_CURV_3D() { Init(); } 
	T_MDGN_PM_CURV_3D& operator= (const T_MDGN_PM_CURV_3D &src) { return *this; }
};


class CCalcBasePmcv  
{
public:
	CCalcBasePmcv();
	virtual ~CCalcBasePmcv();

	BOOL Get2DPMCurveData(T_MDGN_PM_CURV &rPmcv);
	BOOL Get3DPMCurveData(T_MDGN_PM_CURV_3D &rPmcv);

	BOOL GetPMCurveData(double dAngle, T_MDGN_PM_CURV &rPmData);
	// calculate nominal strength (Pn, Mn*)
	BOOL CalNomStrength(double dXn, double dAngleNA, T_MDGN_PM_SECT &rPmSect);
	// calculate design strength (phiPn, phiMn*) : Find design strength for same eccentricity (Mu/Pu).
	BOOL CalDgnStrengthByEccnRatio(double dPu, double dMuy, double dMuz, const T_MDGN_PM_CURV &PmData, T_MDGN_PM_RES &rPmRes);
	// calculate design strength (phiPn, phiMn*) : Find design strength for same Pu.
	BOOL CalDgnStrengthByPu(double dPu, double dMuy, double dMuz, const T_MDGN_PM_CURV &PmData, T_MDGN_PM_RES &rPmRes);
	
protected:
	double CalDmax(double dAngleNA, UINT unSectShape, double dH, double dB);
	BOOL IsPureComp(double dXn);
	BOOL IsPureTens(double dXn);
	BOOL CalRebarDetail(double dXn, double dAngleNA, CArray<T_MDGN_PM_RBAR,T_MDGN_PM_RBAR&> raPmRbar);
	BOOL CalDistanceFromSectionToAngleLine(double dx, double dy, double dAngle, double& dDistMax, double& dDistMin);
	double CalStrengthRatioByEccnRatio(double dPu, double dMuy, double dMuz, double dPhiPnmax, const T_MDGN_PM_CURV &PmData, double &rdPhi, double &rdPhiPn, double &rdPhiMny, double &rdPhiMnz);
	double CalStrengthRatioByPu(double dPu, double dMuy, double dMuz, double dPhiPnmax, const T_MDGN_PM_CURV &PmData, double &rdPhi, double &rdPhiPn, double &rdPhiMny, double &rdPhiMnz);
	
protected:
	double CalRotateAngle(double dMuy, double dMuz);
	double CalMaxDistFromRebarToNeutralAxis(double dXn, double dRota, double dDmax);

protected: // need to override.
	virtual double CalFc()  { return 0.0; }   // design compressive strength of concrete.
	virtual double CalFyr() { return 0.0; }   // design yield strength of rebar.
	virtual double CalAlpha() { return 0.0; } //
	virtual double CalBeta()  { return 0.0; } // 
	virtual double Calesi()   { return 0.0; } // 
	virtual double GetPnmaxReductionFactor() { return 0.0; }

protected:
	// section input data.
	CArray<T_MDGN_RBAR_POS_UNIT, T_MDGN_RBAR_POS_UNIT&> m_aRbarPos;
	CArray<T_MDGN_CONC_POS_UNIT, T_MDGN_CONC_POS_UNIT&> m_aConcPos;
	
};

#endif // !defined(AFX_CALCBASEPMCV_H__825CF04D_530E_40F4_B15F_BCDD24E58BCB__INCLUDED_)
