// CalcBasePmcv.cpp: implementation of the CCalcBasePmcv class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Pmcv.h"
#include "CalcBasePmcv.h"
#include <math.h>
#include "MathFunc.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

const double cZeroLimit = 1.0E-07;
const double cUpOnLimit = 1.0E+28;
const double cP0Tor = 0.001;

#define DGN_ITER_LIMIT 30

#define DGN_SECT_SHAPE_INDEX_REG_L     0
#define DGN_SECT_SHAPE_INDEX_REG_C     1
#define DGN_SECT_SHAPE_INDEX_REG_H     2
#define DGN_SECT_SHAPE_INDEX_REG_T     3
#define DGN_SECT_SHAPE_INDEX_REG_B     4
#define DGN_SECT_SHAPE_INDEX_REG_P     5
#define DGN_SECT_SHAPE_INDEX_REG_SR    6
#define DGN_SECT_SHAPE_INDEX_REG_SB    7
#define DGN_SECT_SHAPE_INDEX_REG_2L    8
#define DGN_SECT_SHAPE_INDEX_REG_2C    9
#define DGN_SECT_SHAPE_INDEX_REG_CC    10
#define DGN_SECT_SHAPE_INDEX_REG_URIB  11
#define DGN_SECT_SHAPE_INDEX_REG_OCT   12
#define DGN_SECT_SHAPE_INDEX_REG_SOCT  13
#define DGN_SECT_SHAPE_INDEX_REG_TRK   14
#define DGN_SECT_SHAPE_INDEX_REG_STRK  15
#define DGN_SECT_SHAPE_INDEX_REG_HTRK  16
#define DGN_SECT_SHAPE_INDEX_REG_CL    17  // Cross Angle // Only Tower
#define DGN_SECT_SHAPE_INDEX_REG_ROCT  18
#define DGN_SECT_SHAPE_INDEX_REG_BSTF  19  // Box with stiffener
#define DGN_SECT_SHAPE_INDEX_REG_PSTF  20  // Pipe with stiffener
#define DGN_SECT_SHAPE_INDEX_REG_GEN   21  // General Section (Value에서만 사용)
#define DGN_SECT_SHAPE_INDEX_REG_CR    22  // Cross (+자 단면) // Only Building
#define DGN_SECT_SHAPE_INDEX_REG_UDT   23  // Upside-Down T

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CCalcBasePmcv::CCalcBasePmcv()
{

}

CCalcBasePmcv::~CCalcBasePmcv()
{

}

BOOL CCalcBasePmcv::CalNomStrength(double dXn, double dAngleNA, T_MDGN_PM_SECT &rPmSect)
{
	rPmSect.Init();

	double dPn  = 0.0;	// Comp(+), Tens(-).
	double dMny = 0.0;
  double dMnz = 0.0;
  double desiMax = 0.0;
  double dTForce = 0.0; 
	double dCForce = 0.0;

	double dYbar = 0.0;  // (input) centroid coordinate y.
	double dZbar = 0.0;  // (input) centroid coordinate z.
	double dFc	 = CalFc();   // design compressive strength of concrete. (according to CODE)
	double dFyr	 = CalFyr();  // design yielding strength of rebar. (according to CODE)

	double dAlpha = CalAlpha();	
	double dH = 0.0;
	double dB = 0.0;
	UINT   unSectShape = 0;
	double dDmax  = CalDmax(dAngleNA, unSectShape, dH, dB);

	if(IsPureComp(dXn) || IsPureTens(dXn))
	{
		double dCb   = 0.0;
		double dPncc = 0.0; // concrete + compression.
		double dPnsc = 0.0; // rebar + compression.
		double dPnst = 0.0; // rebar + tension.
		double dPrsc = 0.0;	// reduction value (rebar + compression)
		
		// by Concrete.
    double dAg = 0.0;  // (input) section gross area.
		if(IsPureTens(dXn))	
		{// Tens.
			dPncc = 0.0;		
			dCb   = 0.0;
			desiMax = cUpOnLimit;
		}
		else if(IsPureTens(dXn))	
		{// Comp.
			dPncc = dAlpha*dFc*dAg;	
			dCb   = dDmax;
			desiMax = -cUpOnLimit;
		}

		// by Rebar.
		int nSizeRbar = m_aRbarPos.GetSize();
		for(int i=0; i<nSizeRbar; ++i)
		{
			double dAs = m_aRbarPos[i].dAs;
			if(IsPureTens(dXn))
			{
				dPnst += (-1.0)*dAs*dFyr;  // (-) sign = Tension.
			}
			else if(IsPureComp(dXn))
			{
				dPnsc += dAs*dFyr;
				dPrsc += dAs*dAlpha*dFc;
			}
		}
		
		if(/*m_bReduceRbar ||*/ IsPureComp(dXn))
		{// 합성단면 계산시 단면공제 여부(FALSE=공제안함, TRUE=공제함), 순수 압축시에는 무조건 공제
			dPncc  -= dPrsc;
		}
		dPn = dPncc + dPnsc + dPnst;
    dTForce = fabs(dPnst);
    dCForce = fabs(dPncc + dPnsc);    	
		
		rPmSect.dAngleNA = dAngleNA;
		rPmSect.dCb    = dCb;
		rPmSect.dXn    = dXn;
    rPmSect.dPn    = dPn;
    rPmSect.dPncc  = dPncc;
    rPmSect.dPnsc  = dPnsc;
    rPmSect.dPnst  = dPnst;
    rPmSect.dMny   = 0.0;
    rPmSect.dMnycc = 0.0;
    rPmSect.dMnysc = 0.0;
    rPmSect.dMnyst = 0.0;
    rPmSect.dMnz   = 0.0;
    rPmSect.dMnzcc = 0.0;
    rPmSect.dMnzsc = 0.0;
    rPmSect.dMnzst = 0.0; 
		rPmSect.desiMax= desiMax;  
		
		/*{ if(Get_Rebar_DatilData(dXn, dInRota*CMathFunc::Get_trang(), ResData.arRbarRes, dAsThk)) ResData.nType = 3; }*/
	}
	else
	{
		double dPncc=0.0, dPnsc=0.0, dPnst=0.0, dPrsc=0.0;	// c/s(Concrete/Steel), c/t(Compression/Tension).
		double dMnccy=0.0, dMnscy=0.0, dMnsty=0.0, dMrscy=0.0;
    double dMnccz=0.0, dMnscz=0.0, dMnstz=0.0, dMrscz=0.0;
		double dTotNnccDistY=0.0, dTotNnccDistZ=0.0;
		double dTotNnscDistY=0.0, dTotNnscDistZ=0.0;
		double dTotNnstDistY=0.0, dTotNnstDistZ=0.0;
		double dTotNrscDistY=0.0, dTotNrscDistZ=0.0;	// Nrsc in Compressive Concrete.
		
		double dBeta	= CalBeta();
		double dRota	= dAngleNA;
		double dDcomp	= dDmax*(1.0-dXn);
		double dDcent = dDmax*0.5;
		double dMaxDistTensRbar = CalMaxDistFromRebarToNeutralAxis(dXn, dRota, dDmax);

		if(fabs(dAngleNA) < cZeroLimit) // 0 Degree.
		{
			double dDz1 = dDmax*dXn;
      double dDy2 = 1.0/cZeroLimit;
			double dDz2 = (dDmax-dBeta*dDcomp);
// 			
// 			//////////////////////////
// 			// by Concrete.
// 			_UMD_RC_CONC_UNIT	ConcUnit;
// 			
// 			double dArea = 0.0;
// 			int iIndexC=0;
// 			POSITION PosC = m_ConcData.arConcUnit.GetStartPosition();
// 			while(PosC)
// 			{
// 				ConcUnit.Initialize();
// 				m_ConcData.arConcUnit.GetNextAssoc(PosC,iIndexC,ConcUnit);
// 				// Check if Concrete Compression Region is OK.
// 				double dCalcAc=0.0, dyInc=0.0, dzInc=0.0;
// 				if(Is_ConcComp(iIndexC,dRota,ConcUnit,dDy2,dDz2,dCalcAc,dyInc,dzInc))
// 				{
// 					double dCalcZc = ConcUnit.dyz[1] + dzInc;
// 					double dPncci = dAlpha*dFc * dCalcAc;
//           dArea += dCalcAc;
// 					dPncc += dPncci;
// 					dTotNnccDistZ += fabs(dPncci * dCalcZc);
// 				}
// 			}
// 			//////////////////////////
// 			// by Rebar.
// 			_UMD_RC_RBAR_UNIT	RbarUnit;
// 						
// 			int iIndexR=0;
// 			POSITION PosR = m_RbarData.arRbarUnit.GetStartPosition();
// 			while(PosR)
// 			{
// 				RbarUnit.Initialize();
// 				m_RbarData.arRbarUnit.GetNextAssoc(PosR,iIndexR,RbarUnit);
// 				double dAs=0.0;
// 				if(IsRCDgn())	dAs = (RbarUnit.iDgnType==0) ? RbarUnit.dLen * dAsThk : RbarUnit.dDgnArea;
// 				else					dAs = RbarUnit.dArea;
// 				double dDist = fabs(dDz1-RbarUnit.dyz[1]);
//         
// 				// Check if compression or not.
// 				if(dDz1 < RbarUnit.dyz[1])	// Compression.
// 				{
//           double desi = Calc_esi(dDcomp, dDist);
// 					double dPnsci = dAs * Calc_Fsi(TRUE, desi);
// 					dPnsc += dPnsci;
// 					dTotNnscDistZ += fabs(dPnsci * RbarUnit.dyz[1]);
// 				}
// 				else	// Tension.
// 				{
//           double desi = Calc_esi(dDcomp, dDist);
//           desiMax = max(desi, desiMax);
// 					double dPnsti = dAs * Calc_Fsi(FALSE, desi);
// 					dPnst += dPnsti;
// 					dTotNnstDistZ += fabs(dPnsti * RbarUnit.dyz[1]);
// 				}
// 				
// 				if(dDz2 < RbarUnit.dyz[1])	// Compression.
// 				{
// 					double dPrsci = dAlpha*dFc * dAs;
//           dPrsc += dPrsci;
// 					dTotNrscDistZ += fabs(dPrsci * RbarUnit.dyz[1]);
// 				}
// 			}
// 			// Mncc, Mnsc, Mnst, Mrsc.
// 			if(fabs(dPncc) > 0.0)	dMnccy = fabs(dPncc * (dZbar-fabs(dTotNnccDistZ/dPncc)));
// 			if(fabs(dPnsc) > 0.0)	dMnscy = fabs(dPnsc * (dZbar-fabs(dTotNnscDistZ/dPnsc)));
// 			if(fabs(dPnst) > 0.0)	dMnsty = fabs(dPnst * (dZbar-fabs(dTotNnstDistZ/dPnst)));
// 			if(fabs(dPrsc) > 0.0)	dMrscy = fabs(dPrsc * (dZbar-fabs(dTotNrscDistZ/dPrsc)));
		}
		else if(fabs(dAngleNA-CMathFunc::m_pi/2.0) < cZeroLimit) // 90 Degree.
		{

		}
		else // 0~90 Degree.
		{

		}
	}

	return TRUE;
}

BOOL CCalcBasePmcv::CalDgnStrengthByEccnRatio(double dPu, double dMuy, double dMuz, const T_MDGN_PM_CURV &PmData, T_MDGN_PM_RES &rPmRes)
{
	T_MDGN_PM_SECT rPmSect;

	double dDmax = PmData.dDmax;
	double dPhiPnmax = 0.0;//!!!CalPhiPnmax();

	double dPhiMn = 0.0;
  double dc = 0.0, dCForce = 0.0, dTForce = 0.0;
  double dRatP = 0.0, dRatMy = 0.0, dRatMz = 0.0, dRatMn = 0.0;
	
	BOOL bComp = dPu > 0.0 ? TRUE : FALSE;
  double dMu = sqrt(pow(dMuy, 2.0) + pow(dMuz, 2.0));
	
	if(fabs(dMu) < cZeroLimit)
	{
		if(fabs(dPu) > cZeroLimit) // pure tension or compression.
		{
			int nSizeUnit = PmData.aPmUnit.GetSize();
			T_MDGN_PM_UNIT appPmUnit;
			if(PmData.aPmUnit.GetSize() > 0) 
			{
				appPmUnit = (dPu > 0.0) ? PmData.aPmUnit[0] : PmData.aPmUnit[nSizeUnit-1];
			}			

			double dPhiMny = 0.0;
			double dPhiMnz = 0.0;
			double dPhiPn  = appPmUnit.dPhi * appPmUnit.dPn;
			if(dPu > 0.0) dPhiPn *= GetPnmaxReductionFactor();

			double dRatio  = dPhiPn==0.0 ? 0.0 : fabs(dPu/dPhiPn);

			rPmSect.dAngleNA = PmData.dAngleNA;
			rPmSect.dCb      = bComp ? PmData.dDmax : 0.0;
			rPmSect.dXn      = appPmUnit.dXn;
			rPmSect.dPn      = appPmUnit.dPn;
			rPmSect.dMny     = 0.0;
			rPmSect.dMnz     = 0.0;
			rPmSect.desiMax  = appPmUnit.desi;
		}
		else // not applied force.
		{
			//!!! set available result.
		}
	}
	else
	{
    double dPXn  = 0.0;
		double dNXn  = 0.0;
	  double dEcc1 = ( 1.01) * cUpOnLimit;
	  double dEcc2 = (-1.01) * cUpOnLimit;

		// Get Eccu.
		double dEccu = dPu/sqrt(pow(dPu,2.0)+pow(dMuy,2.0)+pow(dMuz,2.0));
    if(dEccu == 0.0)  dEccu = cZeroLimit;

		int nSizeUnit = PmData.aPmUnit.GetSize();
		for(int i=0 ; i<nSizeUnit ; i++)
		{
			T_MDGN_PM_UNIT PmUnit = PmData.aPmUnit[i];

			double dPhiPn  = PmUnit.dPhi * PmUnit.dPn;
      double dPhiMny = PmUnit.dPhi * PmUnit.dMny;
      double dPhiMnz = PmUnit.dPhi * PmUnit.dMnz;
      double dPhiMn  = sqrt(dPhiMny*dPhiMny + dPhiMnz*dPhiMnz);
			// Get Eccn.
			double dPhiPMM = sqrt(pow(dPhiPn,2.0)+pow(dPhiMny,2.0)+pow(dPhiMnz,2.0));
      double dEccn   = dPhiPMM==0.0 ? 0.0 : dPhiPn/dPhiPMM;
			if(dEccn >= dEccu && dEccn < dEcc1)	{dEcc1 = dEccn; dPXn = PmUnit.dXn;}
			if(dEccn <= dEccu && dEccn > dEcc2)	{dEcc2 = dEccn; dNXn = PmUnit.dXn;}
		}

		double dPhiPnmax = 0.0;
		double dPhi    = 0.0;
		double dPhiPn  = 0.0;
		double dPhiMny = 0.0;
		double dPhiMnz = 0.0;
		double dRatio  = 0.0;

    if(dNXn == 0.0 || dNXn == 1.0)
    {
      double dRatio = CalStrengthRatioByEccnRatio(dPu, dMuy, dMuz, dPhiPnmax, PmData, dPhi, dPhiPn, dPhiMny, dPhiMnz);
      dc = dDmax;
	    // Get Ratios.
	    dRatP  = (dPhiPn ==0.0 ? 0.0 : fabs(dPu/dPhiPn));
	    dRatMy = (dPhiMny==0.0 ? 0.0 : fabs(dMuy/dPhiMny));
      dRatMz = (dPhiMnz==0.0 ? 0.0 : fabs(dMuz/dPhiMnz));
      dPhiMn = sqrt(pow(dPhiMny,2.0)+pow(dPhiMnz,2.0));
      dRatMn = (dPhiMn==0.0 ? 0.0 : fabs(dMu/dPhiMn));
      double dRatM = max(dRatMn, max(dRatMy, dRatMz));
	    dRatio = max(dRatP, dRatM);			
    }
    else
    {
	    double dPn=0.0, dMn=0.0, dMny=0.0, dMnz=0.0, desiMax=0.0;
      double dPb=0.0, dMb=0.0;
            
			int nSize = PmData.aPmUnit.GetSize();
			for(int i=0 ; i<nSize ; i++)
			{
				T_MDGN_PM_UNIT PmUnit = PmData.aPmUnit.GetAt(i);
				dPn  = PmUnit.dPn;
		    dMny = PmUnit.dMny;
        dMnz = PmUnit.dMnz;
        dMn = sqrt(dMny*dMny + dMnz*dMnz);
				// Check Pb, Mb (Balanced condition).
				if(dMn > dMb)	{dPb = dPn;	dMb = dMn;}
			}
      double dStaXn = dPXn;
      double dEndXn = dNXn;
      double dXn = 0.0;
      int nIterNum = 0;
	    while(TRUE)
	    {
        dXn = (dStaXn + dEndXn)/2.;
				CalNomStrength(dXn, PmData.dAngleNA, rPmSect);
	      
        //dPhi = Calc_Phi(ResData.dPn, dPb, ResData.desiMax); !!!
        dPhiPn  = dPhi * rPmSect.dPn;
        dPhiMny = dPhi * rPmSect.dMny;
        dPhiMnz = dPhi * rPmSect.dMnz;
        double dEccn = dPhiPn/sqrt(pow(dPhiPn,2.0)+pow(dPhiMny,2.0)+pow(dPhiMnz,2.0));

        if(fabs(dEccn/dEccu - 1.0) <= 0.01)
          break;
        else if(nIterNum > DGN_ITER_LIMIT)
          break;
        if(dEccn > dEccu) dStaXn = dXn;
        else              dEndXn = dXn;

        nIterNum++;
      }
  
      // Compare pPn with pPnmax.
	    if(dPhiPn > dPhiPnmax)
			{
				double dRatPnmax = dPhiPn==0.0 ? 0.0 : dPhiPnmax/dPhiPn;
				dPhiMny *= dRatPnmax;
				dPhiMnz *= dRatPnmax;
				dPhiPn  *= dRatPnmax;
			}

      dc = dDmax*(1.0-dXn);
	    // Get Ratios.
	    dRatP  = (dPhiPn ==0.0 ? 0.0 : fabs(dPu/dPhiPn));
	    dRatMy = (dPhiMny==0.0 ? 0.0 : fabs(dMuy/dPhiMny));
      dRatMz = (dPhiMnz==0.0 ? 0.0 : fabs(dMuz/dPhiMnz));
      dPhiMn = sqrt(pow(dPhiMny,2)+pow(dPhiMnz,2));
      dRatMn = (dPhiMn==0.0 ? 0.0 : fabs(dMu/dPhiMn));
      double dRatM = max(dRatMn, max(dRatMy, dRatMz));
	    dRatio = max(dRatP,dRatM);
    }
	}

	return TRUE;
}

BOOL CCalcBasePmcv::CalDgnStrengthByPu(double dPu, double dMuy, double dMuz, const T_MDGN_PM_CURV &PmData, T_MDGN_PM_RES &rPmRes)
{
	double dPhi = 0.0;
	double dPhiPn = 0.0; 
	double dPhiMny = 0.0; 
	double dPhiMnz = 0.0; 
	double dRatio = 0.0;

	double dDmax = PmData.dDmax;
  double dPhiPnMax = 0.0;//Get_PhiPnmax(PmData);  
  double dPhiMn = 0.0;
  double dc = 0.0, dCForce = 0.0, dTForce = 0.0;
  double dRatP = 0.0, dRatMy = 0.0, dRatMz = 0.0, dRatMn = 0.0;
	
  double dMu = sqrt(pow(dMuy,2)+pow(dMuz,2));
	
	T_MDGN_PM_UNIT	PmUnit;
  T_MDGN_PM_SECT ResData;  ResData.Init();

  double dPXn=0.0, dNXn=0.0;
	double dPpPn=cUpOnLimit, dNpPn= -1.0*cUpOnLimit;

	int nSize = PmData.aPmUnit.GetSize();
	for(int i=0 ; i<nSize ; i++)
	{
		PmUnit = PmData.aPmUnit.GetAt(i);	
		double dpPn = PmUnit.dPhi * PmUnit.dPn;
		if(dpPn >= dPu && dpPn < dPpPn)	{dPpPn = dpPn; dPXn = PmUnit.dXn;}
		if(dpPn <= dPu && dpPn > dNpPn)	{dNpPn = dpPn; dNXn = PmUnit.dXn;}
	}

  if(dNXn == 0.0 || dNXn == 1.0)
  {
    CalStrengthRatioByPu(dPu, dMuy, dMuz, dPhiPnMax, PmData, dPhi, dPhiPn, dPhiMny, dPhiMnz);
    dc = dDmax;
	  // Get Ratios.
	  dRatP = (dPhiPn==0.0 ? 0.0 : fabs(dPu/dPhiPn));
	  dRatMy = (dPhiMny==0.0 ? 0.0 : fabs(dMuy/dPhiMny));
    dRatMz = (dPhiMnz==0.0 ? 0.0 : fabs(dMuz/dPhiMnz));
    dPhiMn = sqrt(pow(dPhiMny,2)+pow(dPhiMnz,2));
    dRatMn = (dPhiMn==0.0 ? 0.0 : fabs(dMu/dPhiMn));
    double dRatM = max(dRatMn, max(dRatMy, dRatMz));
	  dRatio = max(dRatP,dRatM);			
  }
  else
  {
	  double dPn=0.0, dMn=0.0, dMny=0.0, dMnz=0.0, desiMax=0.0;
    double dPb=0.0, dMb=0.0;
          
		int nSize = PmData.aPmUnit.GetSize();
		for(int i=0 ; i<nSize ; i++)
		{
			PmUnit = PmData.aPmUnit.GetAt(i);
			dPn  = PmUnit.dPn;
		  dMny = PmUnit.dMny;
      dMnz = PmUnit.dMnz;
      dMn = sqrt(dMny*dMny + dMnz*dMnz);
			// Check Pb, Mb (Balanced condition).
		  if(dMn > dMb)	{dPb = dPn;	dMb = dMn;}
		}
    double dStaXn = dPXn;
    double dEndXn = dNXn;
    double dXn = 0.0;
    int nIterNum = 0;
	  while(TRUE)
	  {
      dXn = (dStaXn + dEndXn)/2.;
	    CalNomStrength(dXn, PmData.dAngleNA, ResData);

      dPhi = 0.0;//Calc_Phi(ResData.dPn, dPb, ResData.desiMax);
      dPhiPn  = dPhi * ResData.dPn;
      dPhiMny = dPhi * ResData.dMny;
      dPhiMnz = dPhi * ResData.dMnz;


      if(fabs(dPhiPn - dPu) <= cP0Tor*dPhiPnMax)
        break;
      else if(nIterNum > DGN_ITER_LIMIT)
        break;
      if(dPhiPn > dPu)
        dStaXn = dXn;
      else
        dEndXn = dXn;

      nIterNum++;
    }

    // Compare pPn with pPnmax.
	  if(dPhiPn > dPhiPnMax)
		{
			dPhiMny = 0.0;
			dPhiMnz = 0.0;
			dPhiPn  = dPhiPnMax;
		}

    dc = dDmax*(1.0-dXn);
	  // Get Ratios.
	  dRatP = (dPhiPn==0.0 ? 0.0 : fabs(dPu/dPhiPn));
	  dRatMy = (dPhiMny==0.0 ? 0.0 : fabs(dMuy/dPhiMny));
    dRatMz = (dPhiMnz==0.0 ? 0.0 : fabs(dMuz/dPhiMnz));
    dPhiMn = sqrt(pow(dPhiMny,2)+pow(dPhiMnz,2));
    dRatMn = (dPhiMn==0.0 ? 0.0 : fabs(dMu/dPhiMn));
    double dRatM = max(dRatMn, max(dRatMy, dRatMz));
	  dRatio = max(dRatP,dRatM);
  }
	
	return TRUE;
}

double CCalcBasePmcv::CalDmax(double dAngleNA, UINT unSectShape, double dH, double dB)
{
	if(unSectShape==5/*SHAPE_SR*/ || unSectShape==6/*SHAPE_P*/) 
	{
		return dH;
	}
	else
	{
		return dB*sin(dAngleNA) + dH*cos(dAngleNA);		
	}	
}

BOOL CCalcBasePmcv::IsPureComp(double dXn)
{
	return (dXn < 0.0) ? TRUE : FALSE;
}

BOOL CCalcBasePmcv::IsPureTens(double dXn)
{
	return (fabs(dXn-1.0) < cZeroLimit) ? TRUE : FALSE;
}

BOOL CCalcBasePmcv::CalRebarDetail(double dXn, double dAngleNA, CArray<T_MDGN_PM_RBAR,T_MDGN_PM_RBAR&> raPmRbar)
{
	raPmRbar.RemoveAll();

	double dyc = 0.0; // centroid coordinate y.
	double dzc = 0.0; // centroid coordinate z.
	double dymax = dyc;
	double dymin = dyc;
	double dzmax = dzc;
	double dzmin = dzc;

	BOOL bGenSect = FALSE;
	if(bGenSect)
	{
		ASSERT(0);
	}
	else
	{
		if(!CalDistanceFromSectionToAngleLine(0, 0, dAngleNA, dzmax, dzmin)) return FALSE;
	}

// 	// 중립축 위치계산
// 	double dNAxisZ = dZmin + dXn*(dZmax-dZmin);
// 	double dDcomp  = dZmax - dNAxisZ;
// 	// Rebar
// 	nSize = m_RbarData.arRbarUnit.GetCount();
// 	arRbarRes.SetSize(nSize);
// 	
// 	_UMD_RC_RBAR_UNIT    RbarUnit;
// 	_UMD_RC_COL_RBAR_RES RbarRes;
// 	double dDist, desi;
// 	BOOL bComp;
// 	int iIndex=0;
// 	POSITION Pos = m_RbarData.arRbarUnit.GetStartPosition();
// 	int nCount = 0;
// 	while(Pos)
// 	{
// 		RbarUnit.Initialize();
// 		m_RbarData.arRbarUnit.GetNextAssoc(Pos, iIndex, RbarUnit);
// 		
// 		dx = RbarUnit.dyz[0];
// 		dy = RbarUnit.dyz[1];
//     dz = 0.0;
// 		
//     // 중립축을 반시계 방향으로 돌릴려면 단면좌표를 시계방향으로 돌려야함
//     CMathFunc::mathRotate(-dRotation, dCenter_y, dCenter_z, 0.0, 0.0, 0.0, 1.0, dx, dy, dz); 
// 		
// 		dDist = fabs(dy-dNAxisZ);
// 		bComp = (dy >= dNAxisZ);
// 		desi  = Calc_esi(dDcomp, dDist);
// 		
//     RbarRes.Initialize();
// 		RbarRes.dyz[0]  = RbarUnit.dyz[0];
// 		RbarRes.dyz[1]  = RbarUnit.dyz[1];
//     RbarRes.dds     = dDist;   // Add Seungjun '20120203 for eGen
//     RbarRes.desi    = desi; // Add Seungjun '20120203 for eGen
// 		if(IsRCDgn()) RbarRes.dAs = (RbarUnit.iDgnType==0) ? dAsThk*RbarUnit.dLen : RbarUnit.dDgnArea;
// 		else          RbarRes.dAs = RbarUnit.dArea;		
// 		RbarRes.dStress = Calc_Fsi(bComp, desi);
// 		RbarRes.dForce  = RbarRes.dAs * RbarRes.dStress;
// 		
// 		arRbarRes.SetAt(nCount, RbarRes);
// 		nCount++;
// 		if(nCount >= nSize)
// 			break;
// 	}	

	return TRUE;
}

BOOL CCalcBasePmcv::CalDistanceFromSectionToAngleLine(double dx, double dy, double dAngle, double& dDistMax, double& dDistMin)
{
	/*
  int iSectShp = m_iSectShp;

  // 형상에 따른 좌표수 정의
  int iPointSize=0, iCirclePointSize=0;
	if     (iSectShp==  4 || iSectShp==  7)	{iPointSize =  4;  iCirclePointSize =  0;}	// Box, Solid Box
	else if(iSectShp==  5 || iSectShp==  6)	{iPointSize =  0;  iCirclePointSize =  1;}	// Pipe, Solid Round
  else if(iSectShp== 12 || iSectShp== 13)	{iPointSize =  8;  iCirclePointSize =  0;}	// Octagon, Solid Octagon
  else if(iSectShp== 14 || iSectShp== 15)	{iPointSize =  4;  iCirclePointSize =  2;}	// Track, Solid Track
  else if(iSectShp== 16)                  {iPointSize =  4;  iCirclePointSize =  1;}	// Half Track 
	else ASSERT(0);                                                                     // 기타 해석불가 단면
  
	DgnPoint*  pPointData       = new DgnPoint[iPointSize];       //좌표점
  DgnPoint*  pCirclePointData = new DgnPoint[iCirclePointSize]; //원의 중심좌표 
  double* dCircleRadius    = new double[iCirclePointSize];//원의 반경

  double dHc, dB1, dH, dB, da, db, dShiftY, dShiftZ;
  // 좌표값 설정
  switch(iSectShp)
  {
  case 4:  // Box  
  case 7:  // Solid Box  
    dHc = m_dSize[0];
		dB1 = m_dSize[1];    
		// Point(4).
		pPointData[0].id = 1;	pPointData[0].uv[0] = 0.0;	pPointData[0].uv[1] = 0.0;
		pPointData[1].id = 2;	pPointData[1].uv[0] = dB1;	pPointData[1].uv[1] = 0.0;
		pPointData[2].id = 3;	pPointData[2].uv[0] = dB1;	pPointData[2].uv[1] = dHc;
		pPointData[3].id = 4;	pPointData[3].uv[0] = 0.0;	pPointData[3].uv[1] = dHc;    
    break;
  case 5:  // Pipe
  case 6:  // Solid Round
    dHc = m_dSize[0];
    // CirclePoint(1).
    pCirclePointData[0].id = 1;	pCirclePointData[0].uv[0] = dHc/2.0;	pCirclePointData[0].uv[1] = dHc/2.0;
    dCircleRadius[0] = dHc/2.0;
    break;
  case 12: // Octagon
  case 13: // Solid Octagon
    dH	= m_dSize[0];
		dB	= m_dSize[1];
		da	= m_dSize[2];
		db	= m_dSize[3];
    dShiftY = m_dYbar;
    dShiftZ = 0.0;
		// Point(8).
		pPointData[0].id = 1;	pPointData[0].uv[0] = -dB/2.+da+dShiftY;		pPointData[0].uv[1] = 0.0  +dShiftZ;
		pPointData[1].id = 2;	pPointData[1].uv[0] =  dB/2.-da+dShiftY;		pPointData[1].uv[1] = 0.0  +dShiftZ;
		pPointData[2].id = 3;	pPointData[2].uv[0] =  dB/2.   +dShiftY;		pPointData[2].uv[1] = db   +dShiftZ;
		pPointData[3].id = 4;	pPointData[3].uv[0] =  dB/2.   +dShiftY;		pPointData[3].uv[1] = dH-db+dShiftZ;
		pPointData[4].id = 5;	pPointData[4].uv[0] =  dB/2.-da+dShiftY;		pPointData[4].uv[1] = dH   +dShiftZ;
		pPointData[5].id = 6;	pPointData[5].uv[0] = -dB/2.+da+dShiftY;		pPointData[5].uv[1] = dH   +dShiftZ;
		pPointData[6].id = 7;	pPointData[6].uv[0] = -dB/2.   +dShiftY;		pPointData[6].uv[1] = dH-db+dShiftZ;
		pPointData[7].id = 8;	pPointData[7].uv[0] = -dB/2.   +dShiftY;		pPointData[7].uv[1] = db   +dShiftZ;
    break;
  case 14: // Track
  case 15: // Solid Track    
		dH	= m_dSize[0];
		dB	= m_dSize[1];
    dShiftY = m_dYbar;
		dShiftZ = 0.0;
    // Point(4).
    pPointData[0].id = 1;	pPointData[0].uv[0] = -dB/2.+dH/2.+dShiftY;		pPointData[0].uv[1] = 0.0  +dShiftZ;
    pPointData[1].id = 2;	pPointData[1].uv[0] =  dB/2.-dH/2.+dShiftY;		pPointData[1].uv[1] = 0.0  +dShiftZ;
    pPointData[2].id = 4;	pPointData[2].uv[0] =  dB/2.-dH/2.+dShiftY;		pPointData[2].uv[1] = dH   +dShiftZ;
    pPointData[3].id = 3;	pPointData[3].uv[0] = -dB/2.+dH/2.+dShiftY;		pPointData[3].uv[1] = dH   +dShiftZ;
    
    // CirclePoint(2).
    pCirclePointData[0].id = 1;  pCirclePointData[0].uv[0] = (pPointData[0].uv[0]+pPointData[3].uv[0])/2.0;  pCirclePointData[0].uv[1] = (pPointData[0].uv[1]+pPointData[3].uv[1])/2.0;
    pCirclePointData[1].id = 2;  pCirclePointData[1].uv[0] = (pPointData[1].uv[0]+pPointData[2].uv[0])/2.0;  pCirclePointData[1].uv[1] = (pPointData[1].uv[1]+pPointData[2].uv[1])/2.0;
    dCircleRadius[0] = dCircleRadius[1] = dH/2.0;    
		break;
  case 16: // Half Track        
		dH	= m_dSize[0];
		dB	= m_dSize[1];
    dShiftY = m_dYbar;
    // Point(4).
		pPointData[0].id = 1;	pPointData[0].uv[0] = 0.0;	      pPointData[0].uv[1] = 0.0;
		pPointData[1].id = 2;	pPointData[1].uv[0] = dB-dH/2.; 	pPointData[1].uv[1] = 0.0;
		pPointData[2].id = 3;	pPointData[2].uv[0] = dB-dH/2.;	  pPointData[2].uv[1] = dH;
		pPointData[3].id = 4;	pPointData[3].uv[0] = 0.0;	      pPointData[3].uv[1] = dH;    
    // CirclePoint(1).
    pCirclePointData[0].id = 1;  pCirclePointData[0].uv[0] = (pPointData[1].uv[0]+pPointData[2].uv[0])/2.0;  pCirclePointData[0].uv[1] = (pPointData[1].uv[1]+pPointData[2].uv[1])/2.0;
    dCircleRadius[0] = dH/2.0;      
    break;
  default: // 기타 해석불가 단면
    ASSERT(0);
    return FALSE;
    break;
  }

  double dYbar, dZbar;
  dYbar = m_dYbar;
  dZbar = m_dZbar;  
  
  // 기준 Line 정의
  double dLineI[3]={ dx,  dy, 0.0};
  double dLineJ[3]={ dx-cos(dAngle), dy+sin(dAngle), 0.0};	  
  double dPoint[3]={ dYbar,  dZbar, 0.0};
	double dDist = CMathFunc::mathDistanceFromIntersectPointToLine(dLineI,dLineJ,dPoint);
  dDistMax = dDistMin = dDist;//초기화
  int i;
  for(int i=0 ; i<iPointSize ; i++)
  {
    dPoint[0] = pPointData[i].uv[0];  dPoint[1] = pPointData[i].uv[1];  dPoint[2] = 0.0;
    dDist = CMathFunc::mathDistanceFromIntersectPointToLine(dLineI,dLineJ,dPoint);
    dDistMax = max(dDistMax, dDist);
    dDistMin = min(dDistMin, dDist);
  }
  for(int i=0 ; i<iCirclePointSize ; i++)
  {
    dPoint[0] = pCirclePointData[i].uv[0];  dPoint[1] = pCirclePointData[i].uv[1];  dPoint[2] = 0.0;
    dDist = CMathFunc::mathDistanceFromIntersectPointToLine(dLineI,dLineJ,dPoint);
    if(iSectShp != 16)
    {
      dDistMax = max(dDistMax, dDist+dCircleRadius[i]);
      dDistMin = min(dDistMin, dDist-dCircleRadius[i]);
    }
    else
    {// Half Track 일경우
      dDistMax = max(dDistMax, dDist+dCircleRadius[i]); // 역방향 고려시에는 dDistMin만 고려해야 함      
    }
  }

  delete [] pPointData;
  delete [] pCirclePointData;
  delete [] dCircleRadius;

  if(iPointSize <= 0 && iCirclePointSize <= 0) return FALSE;
	*/
  return TRUE;
}

double CCalcBasePmcv::CalStrengthRatioByEccnRatio(double dPu, double dMuy, double dMuz, double dPhiPnmax, const T_MDGN_PM_CURV &PmData, 
																		 double &rdPhi, double &rdPhiPn, double &rdPhiMny, double &rdPhiMnz)
{
	//ResData.dNARotate = PmData.dRotate;

	double dRatio=0.0;	// max(RatP, RatM).

	double dEcc1 = ( 1.01) * cUpOnLimit;
	double dEcc2 = (-1.01) * cUpOnLimit;
	double dPhi1=0.0, dPhi2=0.0;
	double dpPn1=0.0, dpPn2=0.0;
	double dpMny1=0.0, dpMny2=0.0;
  double dpMnz1=0.0, dpMnz2=0.0;	
  double dMu = sqrt(pow(dMuy,2)+pow(dMuz,2));

	T_MDGN_PM_UNIT PmUnit, PmUnit1, PmUnit2;
	BOOL bPureMoment = FALSE;

	if(fabs(dPu) > cZeroLimit && fabs(dMu) < cZeroLimit)	// Pure Tension or Compression.
	{
		double dPhicmax, dPhitmax;
    double dpPcmax=0.0, dpPtmax=0.0, dpMnmax=0.0;
    double dpPcMny=0.0, dpPtMny=0.0;
    double dpPcMnz=0.0, dpPtMnz=0.0;
	  // 0 Deg.
		int nSize = PmData.aPmUnit.GetSize();
	  for(int i=0 ; i<nSize ; i++)
	  {
		  PmUnit = PmData.aPmUnit.GetAt(i);
		  double dpPn = PmUnit.dPhi * PmUnit.dPn;
      double dMn = sqrt(PmUnit.dMny*PmUnit.dMny+PmUnit.dMnz*PmUnit.dMnz);
		  double dpMn = PmUnit.dPhi * dMn;
		  if(dpMn >= dpMnmax)	dpMnmax = dpMn;
			//if(i==0)				{dPhicmax = PmUnit.dPhi;  dpPcmax = dpPn; dpPcMny = PmUnit.dPhi*PmUnit.dMny; dpPcMnz = PmUnit.dPhi*PmUnit.dMnz;}
      if(PmUnit.dXn==0.0)				{dPhicmax = PmUnit.dPhi;  dpPcmax = dpPn; dpPcMny = PmUnit.dPhi*PmUnit.dMny; dpPcMnz = PmUnit.dPhi*PmUnit.dMnz;}
			if(i==nSize-1)	{dPhitmax = PmUnit.dPhi;  dpPtmax = dpPn; dpPtMny = PmUnit.dPhi*PmUnit.dMny; dpPtMnz = PmUnit.dPhi*PmUnit.dMnz;}
	  }
    //Edit Code By RSH 2002.08.24
	  //dpPcmax *= Get_MaxPnReductionFactor();	// Assume Ties. !!!
		
		rdPhi = (dPu > 0.0 ? dPhicmax : dPhitmax);
    if(!bPureMoment)
    {
		  rdPhiMny = 0.0;
      rdPhiMnz = 0.0;
    }
    else if(dPu > 0.0)
    {
      rdPhiMny = dpPcMny;
      rdPhiMnz = dpPcMnz;
    }
    else
    {
      rdPhiMny = dpPtMny;
      rdPhiMnz = dpPtMnz;
    }
		rdPhiPn = (dPu > 0.0 ? dpPcmax : dpPtmax);
		dRatio  = (rdPhiPn==0.0 ? 0.0 : fabs(dPu/rdPhiPn));

		if(nSize > 0) PmUnit = (dPu > 0.0) ? PmData.aPmUnit.GetAt(0) : PmData.aPmUnit.GetAt(nSize-1);
		else          PmUnit.Init();

// 		.dCb       = (dPu > 0.0) ? PmData.dDmax : 0.0; 
// 		ResData.dXn       = PmUnit.dXn;
// 		ResData.dPn       = PmUnit.dPn;
// 		ResData.dMny      = 0.0; //PmUnit.dMny;
// 		ResData.dMnz      = 0.0; //PmUnit.dMnz;
// 		ResData.desiMax   = PmUnit.desi;
	}
  else if(fabs(dPu) < cZeroLimit && fabs(dMu) < cZeroLimit)	// NO applied force.
	{
// 		dPhi    = 0.0;
// 		dPhiMny = 0.0;
//     dPhiMnz = 0.0;
// 		dPhiPn = 0.0;
// 		dRatio = 0.0;
// 
// 		ResData.nType = 0;
	}
	else
	{
		// Get Eccu.
    double dPuMu = sqrt(pow(dPu,2)+pow(dMuy,2)+pow(dMuz,2));
    double dEccu = (dPuMu==0.0)? 0.0 : dPu/dPuMu;

		int nSize = PmData.aPmUnit.GetSize();
		for(int i=0 ; i<nSize ; i++)
		{
			PmUnit = PmData.aPmUnit.GetAt(i);
			double dpPn  = PmUnit.dPhi * PmUnit.dPn;
      double dpMny = PmUnit.dPhi * PmUnit.dMny;
      double dpMnz = PmUnit.dPhi * PmUnit.dMnz;
      double dpMn = sqrt(dpMny*dpMny+dpMnz*dpMnz);
			// Get Eccn.
      double dpPMM = sqrt(pow(dpPn,2)+pow(dpMny,2)+pow(dpMnz,2));
      double dEccn = (dpPMM==0.0)? 0.0 : dpPn/dpPMM;
			if(dEccn >= dEccu && dEccn < dEcc1)	{dEcc1 = dEccn; dPhi1=PmUnit.dPhi; dpPn1 = dpPn; dpMny1 = dpMny; dpMnz1 = dpMnz; PmUnit1 = PmUnit;}
			if(dEccn <= dEccu && dEccn > dEcc2)	{dEcc2 = dEccn; dPhi2=PmUnit.dPhi; dpPn2 = dpPn; dpMny2 = dpMny; dpMnz2 = dpMnz; PmUnit2 = PmUnit;}
		}
		// Search phiPn, phiMn.
    double dMagnitude = sqrt(pow(dpPn1-dpPn2,2)+pow(dpMny1-dpMny2,2)+pow(dpMnz1-dpMnz2,2));
    if(dMagnitude < cZeroLimit)
    {
			rdPhi    = dPhi2;
      rdPhiPn  = dpPn2;
      rdPhiMny = dpMny2;
      rdPhiMnz = dpMnz2;
			
// 			ResData.dCb       = min(PmData.dDmax, (1.0-PmUnit2.dXn)*PmData.dDmax); 
// 			ResData.dXn       = PmUnit2.dXn;
// 			ResData.dPn       = PmUnit2.dPn;
// 			ResData.dMny      = PmUnit2.dMny;
// 			ResData.dMnz      = PmUnit2.dMnz;
// 			ResData.desiMax   = PmUnit2.desi;
    }
    else
    {
			double dUPhi = (dPhi1-dPhi2) / dMagnitude;
      double dUp = (dpPn1-dpPn2) / dMagnitude;
      double dUmy = (dpMny1-dpMny2) / dMagnitude;
      double dUmz = (dpMnz1-dpMnz2) / dMagnitude;
      double dAmp = 0.0;
      if(fabs(dUp) < cZeroLimit) ASSERT(0);
      if(fabs(dEccu) < cZeroLimit)
      {
        dAmp = -dpPn2/dUp;
				rdPhi    = dPhi2 + dUPhi*dAmp;
        rdPhiPn  = dpPn2 + dUp*dAmp;
        rdPhiMny = dpMny2 + dUmy*dAmp;
        rdPhiMnz = dpMnz2 + dUmz*dAmp;
				
// 				ResData.dXn       = PmUnit2.dXn  + (PmUnit1.dXn -PmUnit2.dXn )/dMagnitude;
// 				ResData.dCb       = min(PmData.dDmax, (1.0-ResData.dXn)*PmData.dDmax); 
// 				ResData.dPn       = PmUnit2.dPn  + (PmUnit1.dPn -PmUnit2.dPn )/dMagnitude;
// 				ResData.dMny      = PmUnit2.dMny + (PmUnit1.dMny-PmUnit2.dMny)/dMagnitude;
// 				ResData.dMnz      = PmUnit2.dMnz + (PmUnit1.dMnz-PmUnit2.dMnz)/dMagnitude;
// 				ResData.desiMax   = PmUnit2.desi + (PmUnit1.desi-PmUnit2.desi)/dMagnitude;
      }
      else
      {
        double dA = (1-1./pow(dEccu,2))*pow(dUp,2) + pow(dUmy,2) + pow(dUmz,2);
        double dB = 2*((1-1./pow(dEccu,2))*dpPn2*dUp + dpMny2*dUmy + dpMnz2*dUmz);
        double dC = (1-1./pow(dEccu,2))*pow(dpPn2,2) + pow(dpMny2,2) + pow(dpMnz2,2);
        double dParam1 = (pow(dB,2)-4*dA*dC < 0.0 ? 0.0 : pow(dB,2)-4*dA*dC);
        double dAmp1 = (-dB+sqrt(dParam1))/(2*dA);
        double dAmp2 = (-dB-sqrt(dParam1))/(2*dA);
				double dPhi_1   = dPhi2 + dUPhi*dAmp1;
        double dPhiPn1  = dpPn2 + dUp*dAmp1;
        double dPhiMny1 = dpMny2 + dUmy*dAmp1;
        double dPhiMnz1 = dpMnz2 + dUmz*dAmp1;
				double dPhi_2   = dPhi2 + dUPhi*dAmp2;
        double dPhiPn2  = dpPn2 + dUp*dAmp2;
        double dPhiMny2 = dpMny2 + dUmy*dAmp2;
        double dPhiMnz2 = dpMnz2 + dUmz*dAmp2;

        if(fabs(dPu) > cZeroLimit)
        {
          if(dPu*dPhiPn1 > 0.0 &&
            dPhiPn1 >= min(dpPn1,dpPn2) && dPhiPn1 <= max(dpPn1,dpPn2) &&
            dPhiMny1 >= min(dpMny1,dpMny2) && dPhiMny1 <= max(dpMny1,dpMny2) &&
            dPhiMnz1 >= min(dpMnz1,dpMnz2) && dPhiMnz1 <= max(dpMnz1,dpMnz2))
          {
						rdPhi    = dPhi_1;
            rdPhiPn  = dPhiPn1;
            rdPhiMny = dPhiMny1;
            rdPhiMnz = dPhiMnz1;
						
// 						ResData.dCb       = min(PmData.dDmax, (1.0-PmUnit1.dXn)*PmData.dDmax); 
// 						ResData.dXn       = PmUnit1.dXn;
// 						ResData.dPn       = PmUnit1.dPn;
// 						ResData.dMny      = PmUnit1.dMny;
// 						ResData.dMnz      = PmUnit1.dMnz;
// 						ResData.desiMax   = PmUnit1.desi;
          }
          else
          {
						rdPhi    = dPhi_2;
            rdPhiPn  = dPhiPn2;
            rdPhiMny = dPhiMny2;
            rdPhiMnz = dPhiMnz2;
						
// 						ResData.dCb       = min(PmData.dDmax, (1.0-PmUnit2.dXn)*PmData.dDmax); 
// 						ResData.dXn       = PmUnit2.dXn;
// 						ResData.dPn       = PmUnit2.dPn;
// 						ResData.dMny      = PmUnit2.dMny;
// 						ResData.dMnz      = PmUnit2.dMnz;
// 						ResData.desiMax   = PmUnit2.desi;
          }
        }
        else
        {
          if(dPhiMny1 >= min(dpMny1,dpMny2) && dPhiMny1 <= max(dpMny1,dpMny2) &&
            dPhiMnz1 >= min(dpMnz1,dpMnz2) && dPhiMnz1 <= max(dpMnz1,dpMnz2))
          {
						rdPhi    = dPhi_1;
            rdPhiPn  = dPhiPn1;
            rdPhiMny = dPhiMny1;
            rdPhiMnz = dPhiMnz1;
						
// 						ResData.dCb       = min(PmData.dDmax, (1.0-PmUnit1.dXn)*PmData.dDmax); 
// 						ResData.dXn       = PmUnit1.dXn;
// 						ResData.dPn       = PmUnit1.dPn;
// 						ResData.dMny      = PmUnit1.dMny;
// 						ResData.dMnz      = PmUnit1.dMnz;
// 						ResData.desiMax   = PmUnit1.desi;
          }
          else
          {
						rdPhi    = dPhi_2;
            rdPhiPn  = dPhiPn2;
            rdPhiMny = dPhiMny2;
            rdPhiMnz = dPhiMnz2;
						
// 						ResData.dCb       = min(PmData.dDmax, (1.0-PmUnit2.dXn)*PmData.dDmax); 
// 						ResData.dXn       = PmUnit2.dXn;
// 						ResData.dPn       = PmUnit2.dPn;
// 						ResData.dMny      = PmUnit2.dMny;
// 						ResData.dMnz      = PmUnit2.dMnz;
// 						ResData.desiMax   = PmUnit2.desi;
          }        
        }
      }
    }

		// Compare pPn with pPnmax.
		double dpPnmax = dPhiPnmax;

		if(rdPhiPn > dpPnmax)	
		{
			double dRatPnmax = rdPhiPn==0.0 ? 0.0 : dpPnmax/rdPhiPn;
			rdPhiMny *= dRatPnmax;
			rdPhiMnz *= dRatPnmax;
			rdPhiPn *= dRatPnmax;
		}
    double dpMn = sqrt(pow(rdPhiMny,2)+pow(rdPhiMnz,2));
		// Get Ratios.
		double dRatP =  (rdPhiPn==0.0 ? 0.0 : fabs(dPu/rdPhiPn));
		double dRatMy = (rdPhiMny==0.0 ? 0.0 : fabs(dMuy/rdPhiMny));
    double dRatMz = (rdPhiMnz==0.0 ? 0.0 : fabs(dMuz/rdPhiMnz));
    double dRatMn = (dpMn==0.0 ? 0.0 : fabs(dMu/dpMn));
    double dRatM = max(dRatMn, max(dRatMy, dRatMz));
		dRatio = max(dRatP,dRatM);
	}
	return dRatio;	
}

double CCalcBasePmcv::CalStrengthRatioByPu(double dPu, double dMuy, double dMuz, double dPhiPnmax, const T_MDGN_PM_CURV &PmData, 
																		 double &rdPhi, double &rdPhiPn, double &rdPhiMny, double &rdPhiMnz)
{
	return 0.0;
}

double CCalcBasePmcv::CalRotateAngle(double dMuy, double dMuz)
{
	double dABSMy = fabs(dMuy);
	double dABSMz = fabs(dMuz);

	double dAngle = 0.0;
	if(dABSMy > 0.0 && dABSMz > 0.0) dAngle = atan(dMuy/dMuz);
	else if(dABSMy == 0.0)           dAngle = 0.0;                 // y-y axis. ( 0 Degree)
	else if(dABSMz == 0.0)           dAngle = CMathFunc::m_pi/2.0; // z-z axis. (90 Degree)
	
	return dAngle;
}

double CCalcBasePmcv::CalMaxDistFromRebarToNeutralAxis(double dXn, double dRota, double dDmax)
{
	// Calculate maximum distance from Rebar to neutral axis.
	double dMaxDist=0.0;
	
// 	T_MDGN_RBAR_POS_UNIT RbarUnit;	
// 	if(fabs(dRota-0.0) < cZeroLimit)	// 0 Deg.
// 	{		
// 		int iIndexR=0;
// 		POSITION PosR = m_RbarData.arRbarUnit.GetStartPosition();
// 		while(PosR)
// 		{
// 			RbarUnit.Initialize();
// 			m_RbarData.arRbarUnit.GetNextAssoc(PosR,iIndexR,RbarUnit);
// 			if(bFixDgnRebarOnly)
// 			{ if(bDgn && RbarUnit.iDgnType==0) continue; }
// 			// Check if compression or not.
// 			double dDz1	= dDmax*dXn;
// 			if(dDz1 > RbarUnit.dyz[1])	// Only Tension.
// 			{
// 				double dDist = dDz1 - RbarUnit.dyz[1];
// 				if(dDist > dMaxDist)	dMaxDist = dDist;
// 			}
// 		}
// 	}
// 	else if(fabs(dRota-Get_Pi()/2.0) < cZeroLimit)	// 90 Deg.
// 	{		
// 		int iIndexR=0;
// 		POSITION PosR = m_RbarData.arRbarUnit.GetStartPosition();
// 		while(PosR)
// 		{
// 			RbarUnit.Initialize();
// 			m_RbarData.arRbarUnit.GetNextAssoc(PosR,iIndexR,RbarUnit);
// 			if(bFixDgnRebarOnly)
// 			{ if(bDgn && RbarUnit.iDgnType==0) continue; }
// 			// Check if compression or not.
// 			double dDy1	= dDmax*dXn;
// 			if(dDy1 > RbarUnit.dyz[0])	// Only Tension.
// 			{
// 				double dDist = dDy1 - RbarUnit.dyz[0];
// 				if(dDist > dMaxDist)	dMaxDist = dDist;
// 			}
// 		}
// 	}
// 	else
// 	{		
// 		int iIndexR=0;
// 		POSITION PosR = m_RbarData.arRbarUnit.GetStartPosition();
// 		while(PosR)
// 		{
// 			RbarUnit.Initialize();
// 			m_RbarData.arRbarUnit.GetNextAssoc(PosR,iIndexR,RbarUnit);
// 			if(bFixDgnRebarOnly)
// 			{ if(bDgn && RbarUnit.iDgnType==0) continue; }
// 			// Check if compression or not.
// 			double dDy1		 = dDmax*dXn / sin(dRota);
// 			double dDz1		 = dDmax*dXn / cos(dRota);
// 			double dCoordZ = (-1)*tan(dRota)*RbarUnit.dyz[0] + dDz1;	// by Neutral axis.
// 			if(dCoordZ > RbarUnit.dyz[1])	// Only Tension.
// 			{
//         if(dDy1 == 0.0 && dDz1 == 0.0)
//         {
// 					dDy1 = cZeroLimit / sin(dRota);
// 					dDz1 = cZeroLimit / cos(dRota);
//         }
// 				
// 				double dLineI[3]={dDy1,  0.0, 0.0};
// 				double dLineJ[3]={ 0.0, dDz1, 0.0};
// 				double dPoint[3]={RbarUnit.dyz[0], RbarUnit.dyz[1], 0.0};
// 				double dDist = CMathFunc::mathDistanceFromIntersectPointToLine(dLineI,dLineJ,dPoint);
// 				if(dDist > dMaxDist)	dMaxDist = dDist;
// 			}
// 		}
// 	}
	return dMaxDist;
}