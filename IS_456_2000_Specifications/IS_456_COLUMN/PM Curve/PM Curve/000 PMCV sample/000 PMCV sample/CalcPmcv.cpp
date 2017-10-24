// CalcPmcv.cpp: implementation of the CCalcPmcv class.
//
//////////////////////////////////////////////////////////////////////



#include "stdafx.h"
#include "Pmcv.h"
#include "CalcPmcv.h"
#include "MdgnPmTool.h"
#include "MdgnSectTool.h"
#include "PolyMaker.h"

#include <math.h>
#include "MathFunc.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CCalcPmcv::CCalcPmcv()
{

}

CCalcPmcv::~CCalcPmcv()
{

}

BOOL CCalcPmcv::SetCommInputData()
{
	m_iPmDivNum = 24;
	m_iHoopType = 1;
	m_bFindMmax = FALSE;
	m_bReduceRbar = FALSE;
	m_iIterNum  = 30;

	m_dPhi[0] = 0.85;
	m_dPhi[1] = 0.85;
	m_dPhi[2] = 0.70;
	m_dPhi[3] = 0.65;
	m_dPhi[4] = 0.75;

	return TRUE;
}

BOOL CCalcPmcv::SetSectInputData(int nSectShap, double dSize[8], const _UMD_RC_COL_MAINRBAR &MainRbar)
{
	m_iSectShp = nSectShap;
	for(int i=0; i<8; ++i)
	{
		m_dSize[i] = dSize[i];
	}

	if(nSectShap==DGN_SECT_SHAPE_INDEX_REG_SB)
	{
		// 우선은 사각형만 Test.
		m_dArea = m_dSize[0]*m_dSize[1];
		m_dYbar = m_dSize[1]*0.5;
		m_dZbar = m_dSize[0]*0.5;	
	}
	else if(nSectShap==DGN_SECT_SHAPE_INDEX_REG_SR)
	{
		m_dArea = pow(m_dSize[0],2.0)*CMathFunc::m_pi/4.0;
		m_dYbar = m_dSize[0]*0.5;
		m_dZbar = m_dSize[0]*0.5;
	}
	
	CMdgnSectTool SectTool;
	SectTool.MakeRbarDataReg(nSectShap, dSize, MainRbar, m_RbarData);
	SectTool.MakeConDataReg(nSectShap, dSize, m_ConcData, 72);

	m_iDivisionNumber_90deg = 6;
	
	return TRUE;
}

BOOL CCalcPmcv::SetSectInputDataGen(int nSectShap, double dSize[8], _UMD_RC_CON_POLY_DATA &ConPoly, _UMD_RC_COL_MAINRBAR &MainRbar)
{
	m_iSectShp    = nSectShap;
	for(int i=0; i<8; ++i)
	{
		m_dSize[i] = dSize[i];
	}	
	m_ConPolyData = ConPoly;
	
	CMdgnSectTool SectTool;	
	SectTool.MakeRbarDataGen(nSectShap, MainRbar, m_RbarData);

	if(m_iPmDivNum%2==0) m_iPmDivNum++;
	m_iTypeStressStrain = 2;

	int nSliceNum    = 50;
	int nDivNum90Deg = 6;

	_UMD_RC_GSEC_POLYGON Poly = ConPoly.OutPoly;
	CArray<_UMD_RC_GSEC_POLYGON, _UMD_RC_GSEC_POLYGON&> aOutPoly;
	aOutPoly.RemoveAll();
	aOutPoly.Add(Poly);

	SectTool.Calc_SectGeneral(aOutPoly, ConPoly.arInnPoly, m_dArea, m_dYbar, m_dZbar);

	m_iSliceNum             = nSliceNum;
	m_iDivisionNumber_90deg = nDivNum90Deg;

	if(m_iTypeStressStrain==3)
		SectTool.Get_RotateSectData(1, m_dYbar, m_dZbar, nSliceNum, nDivNum90Deg, m_ConPolyData, m_RbarData, m_RotSectData);
	else
		SectTool.Get_RotateSectData(2, m_dYbar, m_dZbar, nSliceNum, nDivNum90Deg, m_ConPolyData, m_RbarData, m_RotSectData);

	m_iSymmetryType = m_RotSectData.iSymmetryType;
	

	
	return TRUE;
}

BOOL CCalcPmcv::SetMatlInputData(double dFck, double dEc, double dFyr, double dEs)
{
	m_dFc  = dFck;
	m_dEc  = dEc;
	m_dFyr = dFyr;
	m_dEsr = dEs;

	return TRUE;
}

BOOL CCalcPmcv::SetLoadInputData(double dPu, double dMuy, double dMuz, double dMu)
{
	m_arLoad.RemoveAll();

	_UMD_RC_FORCE Load; Load.Init();
	Load.bUse = TRUE;
	Load.dP  = dPu;
	Load.dMy = dMuy;
	Load.dMz = dMuz;

	m_arLoad.Add(Load);

	return TRUE;
}

BOOL CCalcPmcv::ExecuteCheck(BOOL bAllAxis)
{
  FILE *fout;
  errno_t err = fopen_s(&fout, GetFileName(), "w");
  if(err != 0 || fout==NULL)
	{
		AfxMessageBox("Failed in creating the file!");
		return FALSE;
	}
	//ofstream fout(GetFileName(), ios::trunc);

	SetDesignFlag(FALSE);
	m_3DPmData.Init();
	m_arPmmRes.RemoveAll();
	
	int nLoadSize = m_arLoad.GetSize();
	
	m_arPmmRes.SetSize(nLoadSize);
	
	_UMD_RC_FORCE         LoadData;
	T_PM_CHK_RES   PmmRes;	
	_UMD_3DPM_DATA        PMData3D;        

	BOOL   bIsGenType = IsGeneralType();
	double dAg = m_dArea;
	double dAsCur = GetRebarArea();
	double dRhoMax = GetMaxRebarRatio(m_dMaxRhoUser);
	double dRhoMin = GetMinRebarRatio();

	T_PMCV_2D                  PmcvData;	
	T_PM_SECT_RES        CbPmData;
	
	if(!bIsGenType)
	{	
		//DB형식은 전체 3D가 필요할때에만 생성함		
		if(bAllAxis) Get3DPMCurveDataReg(0.0, m_3DPmData);
	}
	else
	{//Gen형식은 Checking전에 3D PMcurve가 생성되어 있어야 함
		Get3DPMCurveDataGen(0.0, m_3DPmData);

		// Phi를 고려하지 않은 CDgnCalcBase_3DPM_Tool용 3D PM상관도 좌표정보 생성
		CMdgnPmTool PMTool;
		PMTool.Get_ConvertToPMTool(1, m_3DPmData, PMData3D);
	}

	for(int i=0 ; i<nLoadSize ; i++)
	{
		PmmRes.Init();
		// 하중입력///////////////////////////
		LoadData = m_arLoad.GetAt(i);
		if(LoadData.bUse == FALSE) 
		{
			PmmRes.bUse = FALSE;
			m_arPmmRes.SetAt(i, PmmRes);
			continue;
		}
		m_dPu  = LoadData.dP ;
		m_dMu  = LoadData.dM ;
		m_dVu  = LoadData.dV ;
		m_dMuy = LoadData.dMy; 
		m_dMuz = LoadData.dMz;
		m_dVuy = LoadData.dVy;
		m_dVuz = LoadData.dVz;
		m_bEqSpecial = LoadData.bEqSpecial;   //(KCI기준 적용하지 않음)
		m_iSeismicTypeLcom = 0; //(KCI기준 적용하지 않음)

		//휨.압축강도 계산////////////////////
		
		// - Check Rho.
    double dPhiPnmax = 0.0;
	  double dRhoCur = dAsCur/dAg;
	  int iAsRes=0;
	  if(dRhoCur < dRhoMin)				iAsRes = 1;
	  else if(dRhoCur > dRhoMax)	iAsRes = 2;
    
	  double dRatio;
		if(!bIsGenType)
		{	
			// Get PM Data.
			GetRealPMCurveData(0.0, PmcvData, CbPmData, PmmRes.dPhiPnmax);
			// Get phiPn, phiMn, Ratio.
			Check_PnMn(0.0, PmcvData, PmmRes.dPhi, PmmRes.dPhiPn, PmmRes.dPhiMny, PmmRes.dPhiMnz, dRatio);
			PmmRes.dPhiMn = sqrt(pow(PmmRes.dPhiMny,2)+pow(PmmRes.dPhiMnz,2));
			PmmRes.dRatioM = PmmRes.dRatioP = dRatio;			
		}
		else 
		{
			Check_PnMn_Gen(PMData3D, PmmRes.dPhi, PmmRes.dPhiPn, PmmRes.dPhiMny, PmmRes.dPhiMnz, dRatio, PmcvData, CbPmData);
			PmmRes.dPhiMn = sqrt(pow(PmmRes.dPhiMny,2)+pow(PmmRes.dPhiMnz,2));
			PmmRes.dRatioM = PmmRes.dRatioP = dRatio;
			PmmRes.dPhiPnmax = Get_PhiPnmax(PmcvData);
		}
		PmmRes.bUse = TRUE;
		PmmRes.dAs = dAsCur;
		PmmRes.iAsRes = iAsRes;
		PmmRes.dRho      = dAsCur/dAg;
		PmmRes.dRhoMax   = dRhoMax;
		PmmRes.dRhoMin   = dRhoMin;
		PmmRes.bPmcvData = TRUE;
		PmmRes.bCbData   = TRUE;
		PmmRes.PmcvData  = PmcvData;
		PmmRes.CbPmData  = CbPmData;

		//ofstream fout(GetFileName(), ios::app);		
		int nSizePm = PmmRes.PmcvData.arPmUnit.GetSize();
		for(int p=0; p<nSizePm; ++p)
		{
			CString strOut = _T("");
			strOut.Format("%.3f    %.3f    %.3f\n", PmmRes.PmcvData.arPmUnit[p].dPn, PmmRes.PmcvData.arPmUnit[p].dMny, PmmRes.PmcvData.arPmUnit[p].dMnz);
      fprintf(fout, strOut);
			//fout<<strOut<<endl;
		}		
		// - Save Result Data
		m_arPmmRes.SetAt(i, PmmRes);
  }
  fclose(fout);
	//fout.close();

	return TRUE;
}

BOOL CCalcPmcv::IsGeneralType()
{
	if(m_iSectShp == DGN_SECT_SHAPE_INDEX_REG_GEN) return TRUE;
	return FALSE;
}

double CCalcPmcv::GetRebarArea()
{
	double dAs=0.0;
	
	_UMD_RC_RBAR_UNIT RbarUnit;
	int iIndex=0;
	POSITION Pos = m_RbarData.arRbarUnit.GetStartPosition();
	while(Pos)
	{
		RbarUnit.Init();
		m_RbarData.arRbarUnit.GetNextAssoc(Pos, iIndex, RbarUnit);
		dAs += RbarUnit.dArea;
	}
	return dAs;
}

BOOL CCalcPmcv::Get3DPMCurveDataReg(double dAsThk, T_PMCV_3D& PmData3D)
{
	PmData3D.Init();
	if(IsGeneralType()) return FALSE;
	
  m_dTolM0 = min(0.01, m_dTolM0);
	
	double dRotate, dPhiPnmaxUnit;
	T_PMCV_2D PmData;
	T_PM_SECT_RES CbPmData;
	int iRotaNum = m_iDivisionNumber_90deg+1; // 0~90도 까지
	
	PmData3D.b3DData = TRUE;
	PmData3D.bCbData = TRUE;
	PmData3D.bP0Data = FALSE;
	PmData3D.iSymmetryType = 2; //DB단면은 모두 좌우상하 대칭임(HalfTrack은 비대칭이지만 대책이 나올때 까지 보류:Civil초기부터의 문제임)
	PmData3D.iDivisionNumber_90deg = m_iDivisionNumber_90deg;
	PmData3D.arCbPmData.SetSize(iRotaNum);
	PmData3D.arPmcvData.SetSize(iRotaNum);	
	for(int k=0; k<iRotaNum; k++)
	{
		dRotate = (CMathFunc::m_pi/2.0) * (double)(k)/(double)m_iDivisionNumber_90deg;		
		Get_PMCurve(dRotate, PmData, CbPmData, dPhiPnmaxUnit, dAsThk);		
		if(k==0) PmData3D.dPhiPnmax = dPhiPnmaxUnit;
		PmData3D.arPmcvData.SetAt(k, PmData);
		PmData3D.arCbPmData.SetAt(k, CbPmData);
	}
	return TRUE;
}

BOOL CCalcPmcv::Get3DPMCurveDataGen(double dAsThk, T_PMCV_3D& PmData3D)
{
  PmData3D.Init();
	
  m_dTolM0 = min(0.01, m_dTolM0);
	
	int iDivNum = m_iPmDivNum-1;
  if(iDivNum==0)		return FALSE;
	if(iDivNum%2!=0)	return FALSE;	// Always even.
	
	CMdgnSectTool SectTool;
	CPolyMaker PolyMaker;
	
	// 단면의 회전정보가 없으면 값을 채워 넣음
	if(m_iTypeStressStrain==3)
	{// Slice정보는 없어도 됨
		if(m_RotSectData.nDataType != 1 && m_RotSectData.nDataType != 2)
		{	SectTool.Get_RotateSectData(1, m_dYbar, m_dZbar, m_iSliceNum, m_iDivisionNumber_90deg, m_ConPolyData, m_RbarData, m_RotSectData); }
	}
	else
	{// Slice정보는 필수
		if(m_RotSectData.nDataType != 2)
		{ SectTool.Get_RotateSectData(2, m_dYbar, m_dZbar, m_iSliceNum, m_iDivisionNumber_90deg, m_ConPolyData, m_RbarData, m_RotSectData); }
	}
	
	// 3D Pm-curve작성
	int iHoldNumberToPmcvData = Get_HoldNumberToPmcvData();
	PmData3D.b3DData = TRUE;
	PmData3D.bCbData = TRUE;
	PmData3D.bP0Data = TRUE;
	PmData3D.iSymmetryType = m_RotSectData.iSymmetryType;
	PmData3D.iDivisionNumber_90deg = m_iDivisionNumber_90deg;
	PmData3D.arCbPmData.SetSize(iHoldNumberToPmcvData);
	PmData3D.arPmcvData.SetSize(iHoldNumberToPmcvData);
	PmData3D.arP0PmData.SetSize(iHoldNumberToPmcvData);
	int i=0;
	T_PMCV_2D     PmcvData;
	T_PM_SECT_RES CbPmData;
	T_PM_SECT_RES P0PmData;
	double dPhiPnmax;
  for(int i=0 ; i<iHoldNumberToPmcvData ; i++)
  { 
		Get_PMCurve_Gen(i, PmcvData, CbPmData, P0PmData, dPhiPnmax, dAsThk);
		if(i==0) PmData3D.dPhiPnmax = dPhiPnmax;
		PmData3D.arCbPmData.SetAt(i, CbPmData);
		PmData3D.arP0PmData.SetAt(i, P0PmData);	
		PmData3D.arPmcvData.SetAt(i, PmcvData);		
  }	
	
  return TRUE;
}

BOOL CCalcPmcv::GetRealPMCurveData(double dAsThk, T_PMCV_2D& PmcvData, T_PM_SECT_RES& CbPmData, double& dPhiPnmax, int nCheckType/*=0*/)
{
	PmcvData.Init();
	CbPmData.Init();
	dPhiPnmax = 0.0;
	
  double dRotaByPnMn = 0.0, dRotaByNeutral = 0.0;
  double dRotaByForce = Calc_Rota(-1.0);
	double dRotaCal = dRotaByForce;
	double dRotaSta = 0.0;
	double dRotaEnd = CMathFunc::m_pi/2.0;
  double dMu = sqrt(pow(m_dMuy,2)+pow(m_dMuz,2));
	int iRotaCount=0;
	
	T_PM_SECT_RES ResData;  ResData.Init();
  if(dMu > cUMDRC_Zero)
  {
		while(1)	// Unlimited looping.
		{
			// Calculated by Finding Pn,Mny,Mnz.
      Get_PMCurve(dRotaCal, PmcvData, CbPmData, dPhiPnmax, dAsThk);
			
      // Rota by Pn,Mny,Mnz.
      double dPhi=0.0, dPhiPn=0.0, dPhiMn=0.0, dPhiMny=0.0, dPhiMnz=0.0;
			if(nCheckType==0) Calc_Ratio(m_dPu, m_dMuy, m_dMuz, dPhiPnmax, PmcvData, dPhi, dPhiPn, dPhiMny, dPhiMnz);
			else              Calc_Ratio_ToPu(m_dPu, m_dMuy, m_dMuz, dPhiPnmax, PmcvData, dPhi, dPhiPn, dPhiMny, dPhiMnz, ResData);
			if(dPhiMny==0.0 && dPhiMnz==0.0)	dRotaByPnMn = dRotaByForce;
			else if(dPhiMnz==0.0)			        dRotaByPnMn = 0.0;
			else if(dPhiMny==0.0)						  dRotaByPnMn = CMathFunc::m_pi/2.0;
			else														  dRotaByPnMn = atan(fabs(dPhiMnz/dPhiMny));
			
			iRotaCount++;
			// Change by ZINU.('03.01.09). Compare by Tangent value.
			double dTanF = tan(dRotaByForce);
			double dTanR = tan(dRotaByPnMn);
      if(dTanF==0.0)  break;
			if(fabs((dTanF-dTanR)/dTanF) < CONST_UMDRC_dRATLIM3)	break;	// Satisfy.
      if(iRotaCount > m_iIterNum)			       		break;	// Not satisfy.
			// Use by Bi-Section Method.
			if(dRotaByPnMn > dRotaByForce)	dRotaEnd = dRotaCal;
			else												    dRotaSta = dRotaCal;
			dRotaCal = (dRotaSta + dRotaEnd) / 2.0;      
		}
  }
	else
	{
		Get_PMCurve(0.0, PmcvData, CbPmData, dPhiPnmax, dAsThk);
	}
	return TRUE;
}


BOOL CCalcPmcv::Get_PMCurve(double dInRota, T_PMCV_2D& PmcvData, T_PM_SECT_RES& CbPmData, double& dPhiPnmax, double dAsThk)
{
	PmcvData.Init();
	CbPmData.Init();

	BOOL bFixDgnRebarOnly = (m_bDesign && dAsThk<=0.0);

	dPhiPnmax=0.0;

  double dFyr = m_dFyr;
	T_PM_UNIT	PmUnit;

	double dPn=0.0, dMn=0.0;
	double dPb=0.0, dMb=0.0;
	double dPcmax=0.0, dPtmax=0.0;
	int    iPcmaxID=1, iPtmaxID=1;
  double dXn = 0.0;
	int iDivNum = m_iPmDivNum;
	PmcvData.arPmUnit.SetSize(iDivNum+1);
	T_PM_SECT_RES ResData;
	for(int i=0; i<=iDivNum; i++)
	{
		// IRC_21_00일 경우에는 재정의 필요
    if(i==0)	dXn = -0.5;
    else
    {
		  double dX	 = (double)(i-1)/(iDivNum-1);
		  dXn = 2.0202*pow(dX,3) - 3.13131*pow(dX,2) + 2.11111*dX;
    }
    
    if(!Calc_PnMn(dXn, dInRota, dAsThk, ResData)) return FALSE;//상속받은 Class(기준별Class)에서 재정의 되어 있으면 해당 Calc_PnMn함수로 넘어감

		PmUnit.Init();
		PmUnit.dXn	= dXn;
		PmUnit.dPn	= ResData.dPn;
		PmUnit.dMny	= ResData.dMny;
    PmUnit.dMnz	= ResData.dMnz;
    PmUnit.desi = ResData.desiMax;
		dPn = PmUnit.dPn;
    dMn = sqrt(ResData.dMny*ResData.dMny + ResData.dMnz*ResData.dMnz);
		// Check Pb, Mb (Balanced condition).
		if(dMn > dMb)	{dPb = dPn;	dMb = dMn;}
		// Calculate Pcmax, Ptmax.
		if(dPcmax < ResData.dPn)	{dPcmax = ResData.dPn; iPcmaxID=i;}
		if(dPtmax > ResData.dPn)	{dPtmax = ResData.dPn; iPtmaxID=i;}
		PmcvData.arPmUnit.SetAt(i, PmUnit);
	}
  // Set Phi, Mn(Not calculated).
	for(int i=0; i<=iDivNum; i++)
	{
		PmUnit = PmcvData.arPmUnit.GetAt(i);
		PmUnit.dPhi = Calc_Phi(PmUnit.dPn, dPb, PmUnit.desi);
		PmcvData.arPmUnit.SetAt(i, PmUnit);
	}
	
	double dDmax  = Calc_Dmax(dInRota);
	double dDcomp	= dDmax*0.5;
	double dDeff	= dDcomp + Calc_MaxDistTensColmRbar(0.5, dInRota, dDmax, bFixDgnRebarOnly);		
	PmcvData.dRotate = dInRota;
	PmcvData.dDmax   = dDmax;
	PmcvData.dDeff   = dDeff;
	if(!m_bDesign)
	{// Checking상태이고, 작용부재력에 대한 PM-Curve산출시
		if(m_bFindMmax)
		{
			BOOL bStop = FALSE;
			double dXnTol	 = 1.0-7*m_dTolM0;	
			double dXn0Tol = 1.0-1*m_dTolM0;				
			double dMn_y=0.0, dMn_z=0.0;
			double dXnPrev=0.0, dXnNext=1.0;
			double dMnPrev=0.0, dMnNext=0.0;
			
			int iIterCount=0;
			// Bi-section Method.
			while(!bStop)
			{
				// 1st Pmcv.
				dXn = (dXnPrev+dXnNext) / 2.0;
				if(!Calc_PnMn(dXn, dInRota, dAsThk, CbPmData)) return FALSE;				
				dPn = CbPmData.dPn;   dMn_y = CbPmData.dMny;   dMn_z = CbPmData.dMnz;
				dMn = sqrt(dMn_y*dMn_y+dMn_z*dMn_z);
				// Calculate Only if Balanced Condition required.
				
				// Change by ZINU.('03.12.18). Calculate by 3 Points.
				//  M
				//  |--
				//  |   3rd
				//  |     1st
				//  |       2nd
				//  |        |
				//  +------+------ P
				//  |   __/
				//  |--'
				// 2nd Pmcv for Xn Direction.
				double dXn2 = dXn / dXnTol;
				T_PM_SECT_RES resultData2;
				double dPn2=0.0, dMn2=0.0, dMn2_y=0.0, dMn2_z=0.0;
				if(!Calc_PnMn(dXn2, dInRota, dAsThk, resultData2)) return FALSE;		
				dPn2 = resultData2.dPn;   dMn2_y = resultData2.dMny;   dMn2_z = resultData2.dMnz;
				dMn2 = sqrt(dMn2_y*dMn2_y+dMn2_z*dMn2_z);
				// 3rd Pmcv for Xn Direction.
				double dXn3 = dXn * dXnTol;
				T_PM_SECT_RES resultData3;
				double dPn3=0.0, dMn3=0.0, dMn3_y=0.0, dMn3_z=0.0;
				if(!Calc_PnMn(dXn3, dInRota, dAsThk, resultData3)) return FALSE;
				dPn3 = resultData3.dPn;   dMn3_y = resultData3.dMny;   dMn3_z = resultData3.dMnz;
				dMn3 = sqrt(dMn3_y*dMn3_y+dMn3_z*dMn3_z);
				// Compare 1st' with 2nd' and 3rd'.
				if(dMn2 < dMn && dMn3 < dMn)	bStop = TRUE;	// Satisfy.
				else if(dMn2 < dMn && dMn3 > dMn)	dXnNext = dXn;	// Move Compressive Region.       
				else if(dMn2 > dMn && dMn3 < dMn)	dXnPrev = dXn;	// Move Tensile Region.
				else
				{
					if(dMn3 > dMn2)	dXnNext = dXn;	// Move Compressive Region.
					else						dXnPrev = dXn;	// Move Tensile Region.
				}
				iIterCount++;
				if(iIterCount > m_iIterNum)	bStop = TRUE;	// Not Satisfy.
				
				if(bStop)
				{
					// SHIN : dXn의 값이 꼭 최대값을 가지는 것은 아니므로 3개의 평균값을 취함
					//        단 Mn값을 평균을 취하면 최대값보다 작은 값을 가지게 되므로 평균값을 취하지 않음
					dXn = (dXn+dXn2+dXn3)/3.0;
					dPn = (dPn+dPn2+dPn3)/3.0;
          if(!Calc_PnMn(dXn, dInRota, dAsThk, CbPmData)) return FALSE;
				}			      
			}
		}
		else
		{			
			double dXb    = Get_BalancedXb(dDmax, dDeff);
			if(!Calc_PnMn(dXb, dInRota, dAsThk, CbPmData)) return FALSE;
		}
	}
	else
	{
		if(!Calc_PnMn(0.5, dInRota, dAsThk, CbPmData)) return FALSE;// Dgn상태에서는 속도확보를 위해 중립축이 중앙에 있을때의 값을 사용하며 이는 3D PM상관도 망에서 검색시 지표가 되는 값으로만 활용한다.
	}
  
	PmUnit.Init();
	PmUnit = PmcvData.arPmUnit.GetAt(iPcmaxID);
	dPhiPnmax = Get_MaxPnReductionFactor() * PmUnit.dPhi * PmUnit.dPn;

	// Calc. Pn_emin
	T_PM_UNIT Pmemin, PmUnit1, PmUnit2;	
	double dPn1, dPn2;
	int nSize = PmcvData.arPmUnit.GetSize();
	if(nSize > 2)
	{
		PmUnit1 = PmcvData.arPmUnit.GetAt(0);
		for(int i=1; i<nSize; i++)
		{
			PmUnit2 = PmcvData.arPmUnit.GetAt(i);
			dPn1 = PmUnit1.dPhi*PmUnit1.dPn;
			dPn2 = PmUnit2.dPhi*PmUnit2.dPn;
			if(dPhiPnmax <= max(dPn1,dPn2)+cUMDRC_Zero && dPhiPnmax >= min(dPn1,dPn2)-cUMDRC_Zero)
			{
				double dWeight1, dWeight2;
				dWeight1 = dPn1-dPn2==0.0 ? 1.0 : (dPhiPnmax-dPn2)/(dPn1-dPn2);
				dWeight2 = 1.0 - dWeight1;
				Pmemin.dPhi = dWeight1*PmUnit1.dPhi + dWeight2*PmUnit2.dPhi;
				Pmemin.dPn  = dWeight1*PmUnit1.dPn  + dWeight2*PmUnit2.dPn;
				Pmemin.dMny = dWeight1*PmUnit1.dMny + dWeight2*PmUnit2.dMny;
				Pmemin.dMnz = dWeight1*PmUnit1.dMnz + dWeight2*PmUnit2.dMnz;
				Pmemin.dXn  = dWeight1*PmUnit1.dXn  + dWeight2*PmUnit2.dXn;
				Pmemin.desi = dWeight1*PmUnit1.desi + dWeight2*PmUnit2.desi;
				break;
			}
			PmUnit1 = PmUnit2;
		}
		PmcvData.Pmemin = Pmemin;
	}
	return TRUE;
}

BOOL CCalcPmcv::Check_PnMn(double dAsCur, T_PMCV_2D& PmcvData, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, double& dRatio, T_PM_SECT_RES* pResData)
{		
	dPhi = dPhiPn = 0.0; dPhiMny = 0.0; dPhiMnz = 0.0; dRatio = 0.0;

	double dDmax = PmcvData.dDmax;
  double dPhiPnMax = Get_PhiPnmax(PmcvData);  
  double dPhiMn = 0.0;
  double dc = 0.0, dCForce = 0.0, dTForce = 0.0;
  double dRatP = 0.0, dRatMy = 0.0, dRatMz = 0.0, dRatMn = 0.0;

  double dMu = sqrt(pow(m_dMuy,2)+pow(m_dMuz,2));

	T_PM_UNIT	PmUnit;
	T_PM_SECT_RES ResData;  ResData.Init();  ResData.nType=0;

  if(fabs(dMu) < cUMDRC_Zero)
  {
	  if(fabs(m_dPu) > cUMDRC_Zero)	// Pure Tension or Compression.
	  {
			// 0 Deg.
      double dpPcmax=0.0, dpPtmax=0.0, dpMnmax=0.0;
      int nSize=PmcvData.arPmUnit.GetSize();
      for(int i=0 ; i<nSize ; i++)
			{
				PmUnit = PmcvData.arPmUnit.GetAt(i);
				dPhi = PmUnit.dPhi;
		    double dpPn = PmUnit.dPhi * PmUnit.dPn;
        double dMn = sqrt(PmUnit.dMny*PmUnit.dMny+PmUnit.dMnz*PmUnit.dMnz);
		    double dpMn = PmUnit.dPhi * dMn;
		    if(dpMn >= dpMnmax)	dpMnmax = dpMn;
		    if(i==0)				dpPcmax = dpPn;
		    if(i==nSize-1)	dpPtmax = dpPn;
			}	    
      //Edit Code By RSH 2002.08.24
	    dpPcmax *= Get_MaxPnReductionFactor();	// Assume Ties.

		  dPhiMny = 0.0;
      dPhiMnz = 0.0;
		  dPhiPn = (m_dPu > 0.0 ? dpPcmax : dpPtmax);
		  dRatio = (dPhiPn==0.0 ? 0.0 : fabs(m_dPu/dPhiPn));


			if(nSize > 0) PmUnit = (m_dPu > 0.0) ? PmcvData.arPmUnit.GetAt(0) : PmcvData.arPmUnit.GetAt(nSize-1);
			else          PmUnit.Init();
			ResData.nType = 2;
			ResData.dNARotate = PmcvData.dRotate;
			ResData.dCb       = (m_dPu > 0.0) ? PmcvData.dDmax : 0.0; 
			ResData.dXn       = PmUnit.dXn;
			ResData.dPn       = PmUnit.dPn;
			ResData.dMny      = 0.0; //PmUnit.dMny;
			ResData.dMnz      = 0.0; //PmUnit.dMnz;
			ResData.desiMax   = PmUnit.desi;
    }
    else if(fabs(m_dPu) < cUMDRC_Zero && fabs(dMu) < cUMDRC_Zero)	// NO applied force.
	  {
			dPhi    = 0.0;
		  dPhiMny = 0.0;
      dPhiMnz = 0.0;
		  dPhiPn = 0.0;
		  dRatio = 0.0;

			ResData.Init();
			ResData.nType = 2;
	  }

    dRatP = (dPhiPn == 0.0 ? 0.0 :fabs(m_dPu/dPhiPn));
    dRatMy = (dPhiMny == 0.0 ? 0.0 :fabs(m_dMuy/dPhiMny));
    dRatMz = (dPhiMnz == 0.0 ? 0.0 :fabs(m_dMuz/dPhiMnz));
    dPhiMn = sqrt(pow(dPhiMny,2)+pow(dPhiMnz,2));
    dRatMn = (dPhiMn==0.0 ? 0.0 : fabs(dMu/dPhiMn));	
   
    return TRUE;
  }
	else
	{
    double dPXn=0.0, dNXn=0.0;
	  double dEcc1 = ( 1.01) * cUMDRC_Upon;
	  double dEcc2 = (-1.01) * cUMDRC_Upon;

		// Get Eccu.
		double dEccu=m_dPu/sqrt(pow(m_dPu,2)+pow(m_dMuy,2)+pow(m_dMuz,2));
    if(dEccu == 0.0)  dEccu = cUMDRC_Zero;

		int nSize = PmcvData.arPmUnit.GetSize();
		for(int i=0 ; i<nSize ; i++)
		{
			PmUnit = PmcvData.arPmUnit.GetAt(i);	
			double dpPn = PmUnit.dPhi * PmUnit.dPn;
      double dpMny = PmUnit.dPhi * PmUnit.dMny;
      double dpMnz = PmUnit.dPhi * PmUnit.dMnz;
      double dpMn = sqrt(dpMny*dpMny+dpMnz*dpMnz);
			// Get Eccn.
			double dpPMM = sqrt(pow(dpPn,2)+pow(dpMny,2)+pow(dpMnz,2));
      double dEccn = dpPMM==0.0 ? 0.0 : dpPn/dpPMM;
			if(dEccn >= dEccu && dEccn < dEcc1)	{dEcc1 = dEccn; dPXn = PmUnit.dXn;}
			if(dEccn <= dEccu && dEccn > dEcc2)	{dEcc2 = dEccn; dNXn = PmUnit.dXn;}
		}

    double dAsThk = 0.0;
    if(m_bDesign)
    {
      double dRbarLen	    = Get_RbarLenDgn();    
			double dRbarFixArea = Get_RbarFixDgnArea();

	    if(dRbarLen==0.0 && dRbarFixArea==0.0)	return FALSE;//0.0			
			if(dRbarLen==0.0) dAsThk = 0.0;
      else dAsThk = max(0.0, (dAsCur-dRbarFixArea)/dRbarLen);			
    }

    if(dNXn == 0.0 || dNXn == 1.0)
    {
      Calc_Ratio(m_dPu, m_dMuy, m_dMuz, dPhiPnMax, PmcvData, dPhi, dPhiPn, dPhiMny, dPhiMnz, ResData);
      dc = dDmax;
	    // Get Ratios.
	    dRatP = (dPhiPn==0.0 ? 0.0 : fabs(m_dPu/dPhiPn));
	    dRatMy = (dPhiMny==0.0 ? 0.0 : fabs(m_dMuy/dPhiMny));
      dRatMz = (dPhiMnz==0.0 ? 0.0 : fabs(m_dMuz/dPhiMnz));
      dPhiMn = sqrt(pow(dPhiMny,2)+pow(dPhiMnz,2));
      dRatMn = (dPhiMn==0.0 ? 0.0 : fabs(dMu/dPhiMn));
      double dRatM = max(dRatMn, max(dRatMy, dRatMz));
	    dRatio = max(dRatP,dRatM);			
    }
    else
    {
	    double dPn=0.0, dMn=0.0, dMny=0.0, dMnz=0.0, desiMax=0.0;
      double dPb=0.0, dMb=0.0;
            
			int nSize = PmcvData.arPmUnit.GetSize();
			for(int i=0 ; i<nSize ; i++)
			{
				PmUnit = PmcvData.arPmUnit.GetAt(i);
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
      int iIterNum = 0;
	    while(TRUE)
	    {
        dXn = (dStaXn + dEndXn)/2.;
	      Calc_PnMn(dXn, PmcvData.dRotate, dAsThk, ResData);

        dPhi = Calc_Phi(ResData.dPn, dPb, ResData.desiMax);
        dPhiPn  = dPhi * ResData.dPn;
        dPhiMny = dPhi * ResData.dMny;
        dPhiMnz = dPhi * ResData.dMnz;
        double dEccn = dPhiPn/sqrt(pow(dPhiPn,2)+pow(dPhiMny,2)+pow(dPhiMnz,2));

        if(fabs(dEccn/dEccu - 1) <= 0.01)
          break;
        else if(iIterNum > m_iIterNum)
          break;
        if(dEccn > dEccu)
          dStaXn = dXn;
        else
          dEndXn = dXn;

        iIterNum++;
      }
  
      // Compare pPn with pPnmax.
	    if(dPhiPn > dPhiPnMax)
			{
				double dRatPnmax = dPhiPn==0.0 ? 0.0 : dPhiPnMax/dPhiPn;
				dPhiMny *= dRatPnmax;
				dPhiMnz *= dRatPnmax;
				dPhiPn *= dRatPnmax;
			}

      dc = dDmax*(1.0-dXn);
	    // Get Ratios.
	    dRatP = (dPhiPn==0.0 ? 0.0 : fabs(m_dPu/dPhiPn));
	    dRatMy = (dPhiMny==0.0 ? 0.0 : fabs(m_dMuy/dPhiMny));
      dRatMz = (dPhiMnz==0.0 ? 0.0 : fabs(m_dMuz/dPhiMnz));
      dPhiMn = sqrt(pow(dPhiMny,2)+pow(dPhiMnz,2));
      dRatMn = (dPhiMn==0.0 ? 0.0 : fabs(dMu/dPhiMn));
      double dRatM = max(dRatMn, max(dRatMy, dRatMz));
	    dRatio = max(dRatP,dRatM);
    }
  }
	
	if(pResData != NULL) (*pResData) = ResData;
  return TRUE; 
}

BOOL CCalcPmcv::Check_PnMn_Gen(_UMD_3DPM_DATA& PmData3D, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, double& dRatio, T_PMCV_2D& PmcvData, T_PM_SECT_RES& CbPmData, T_PM_SECT_RES* pResData)
{
	CbPmData.Init();
	PmcvData.Init();
	dPhi = 0.0; dPhiPn = 0.0; dPhiMny = 0.0; dPhiMnz = 0.0; dRatio = 0.0;

	CMdgnPmTool PmTool;
	PmTool.Set_Tolerance(m_dTolAngle, m_dTolM0, m_dTolP0);
	double dMu = sqrt(pow(m_dMuy,2)+pow(m_dMuz,2));
	double dPn, dMny, dMnz, desi, dXn;	
	BOOL bCheck;
	int iLeftFitDirID;
	double dModul_DirID;
	
	T_PM_SECT_RES ResData;  ResData.Init();  ResData.nType=0;
  int nRebar = m_RbarData.arRbarUnit.GetCount();
  BOOL bOK = FALSE;
  if(nRebar>0)
  {
    bOK = PmTool.Get_PmcvUnitAt3DPmcvData(PmData3D, m_dPu, m_dMuy, m_dMuz, dPn, dMny, dMnz, desi, dXn, bCheck, iLeftFitDirID, dModul_DirID);
  }
  else
  {
    BOOL bResInStartPMM;
    bOK = PmTool.Get_PmcvUnitAt3DPmcvData(PmData3D, 0.000001, 0.0, 0.0, m_dPu, m_dMuy, m_dMuz, dPn, dMny, dMnz, desi, dXn, bCheck, iLeftFitDirID, dModul_DirID, FALSE, bResInStartPMM);
  }
	//if(PmTool.Get_PmcvUnitAt3DPmcvData(PmData3D, m_dPu, m_dMuy, m_dMuz, dPn, dMny, dMnz, desi, dXn, bCheck, iLeftFitDirID, dModul_DirID))
  if(bOK)
	{		
		CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&> arPmcv;
		_UMD_PMCV_UNIT_BALANCE PmBalance;
		T_PM_UNIT        PmUnit;
		_UMD_PMCV_UNIT_POLY    PmUnitPoly;
		PmTool.Get_PmcvDataToDirID(PmData3D, iLeftFitDirID, dModul_DirID, arPmcv, PmBalance);

		int nSize = arPmcv.GetSize();
		PmcvData.dRotate = PmBalance.dAngleb;
		PmcvData.dDmax   = PmBalance.dDmax;
		PmcvData.dDeff   = PmBalance.dDeffb;
		PmcvData.arPmUnit.SetSize(nSize);
		int i=0;		
		for(int i=0 ; i<nSize ; i++)
		{
			PmUnitPoly = arPmcv.GetAt(i);
			PmUnit.Init();
			PmUnit.dXn  = PmUnitPoly.dXn;
			PmUnit.dPn  = PmUnitPoly.dPn;
			PmUnit.dMny = PmUnitPoly.dMny;
			PmUnit.dMnz = PmUnitPoly.dMnz;
			PmUnit.desi = PmUnitPoly.desi;
			PmcvData.arPmUnit.SetAt(i, PmUnit);
		}

		CbPmData.nType = 2;
		CbPmData.dNARotate = PmBalance.dAngleb;
		CbPmData.dCb       = PmBalance.dCb;
		CbPmData.dXn       = PmBalance.dXb;
		CbPmData.dPn       = PmBalance.dPb;
		CbPmData.dMny      = PmBalance.dMby;
		CbPmData.dMnz      = PmBalance.dMbz;
		CbPmData.desiMax   = PmBalance.desib;

		dPhiPn  = dPn;
		dPhiMny = dMny;
		dPhiMnz = dMnz;

		//Phi 계산
		if(PmData3D.nType == 1)
		{// Phi 미고려 PmData3D사용
			for(int i=0 ; i<nSize ; i++)
			{
				PmUnit = PmcvData.arPmUnit.GetAt(i);
				PmUnit.dPhi = Calc_Phi(PmUnit.dPn, CbPmData.dPn, PmUnit.desi);
				PmcvData.arPmUnit.SetAt(i, PmUnit);
			}
			// CbPmData는 Phi 미고려값이므로 변동없음
			dPhi = Calc_Phi(dPn, CbPmData.dPn, desi);
			dPhiPn  = dPhi*dPn;
			dPhiMny = dPhi*dMny;
			dPhiMnz = dPhi*dMnz;
		}
		else
		{// Phi 고려 PmData3D사용
			for(int i=0 ; i<nSize ; i++)
			{
				PmUnit = PmcvData.arPmUnit.GetAt(i);
				PmUnit.dPhi = Calc_Phi2(PmUnit.dPn, CbPmData.dPn, PmUnit.desi);
				PmUnit.dPn  /= PmUnit.dPhi;
				PmUnit.dMny /= PmUnit.dPhi;
				PmUnit.dMnz /= PmUnit.dPhi;
				PmcvData.arPmUnit.SetAt(i, PmUnit);
			}
			dPhi = Calc_Phi2(CbPmData.dPn, CbPmData.dPn, CbPmData.desiMax);
			CbPmData.dPn  /= dPhi;
			CbPmData.dMny /= dPhi;
			CbPmData.dMnz /= dPhi;
			// dPn, dMny, dMnz는 Phi 고려값이므로 변동없음
		}

		double dPhiPnMax = Get_PhiPnmax(PmcvData);
		// Compare pPn with pPnmax.
	  if(dPhiPn > dPhiPnMax)
		{
			double dRatPnmax = dPhiPn==0.0 ? 0.0 : dPhiPnMax/dPhiPn;
			dPhiMny *= dRatPnmax;
			dPhiMnz *= dRatPnmax;
			dPhiPn *= dRatPnmax;
		}

		// Calc. Pn_emin
		T_PM_UNIT Pmemin, PmUnit1, PmUnit2;	
		double dPn1, dPn2;
		if(nSize > 2)
		{
			PmUnit1 = PmcvData.arPmUnit.GetAt(0);
			for(int i=1; i<nSize; i++)
			{
				PmUnit2 = PmcvData.arPmUnit.GetAt(i);
				dPn1 = PmUnit1.dPhi*PmUnit1.dPn;
				dPn2 = PmUnit2.dPhi*PmUnit2.dPn;
				if(dPhiPnMax <= max(dPn1,dPn2)+cUMDRC_Zero && dPhiPnMax >= min(dPn1,dPn2)-cUMDRC_Zero)
				{
					double dWeight1, dWeight2;
					dWeight1 = dPn1-dPn2==0.0 ? 1.0 : (dPhiPnMax-dPn2)/(dPn1-dPn2);
					dWeight2 = 1.0 - dWeight1;
					Pmemin.dPhi = dWeight1*PmUnit1.dPhi + dWeight2*PmUnit2.dPhi;
					Pmemin.dPn  = dWeight1*PmUnit1.dPn  + dWeight2*PmUnit2.dPn;
					Pmemin.dMny = dWeight1*PmUnit1.dMny + dWeight2*PmUnit2.dMny;
					Pmemin.dMnz = dWeight1*PmUnit1.dMnz + dWeight2*PmUnit2.dMnz;
					Pmemin.dXn  = dWeight1*PmUnit1.dXn  + dWeight2*PmUnit2.dXn;
					Pmemin.desi = dWeight1*PmUnit1.desi + dWeight2*PmUnit2.desi;
					break;
				}
				PmUnit1 = PmUnit2;
			}
			PmcvData.Pmemin = Pmemin;
		}

    // Get Ratios.
	  double dRatP = (dPhiPn==0.0 ? 0.0 : fabs(m_dPu/dPhiPn));
	  double dRatMy = (dPhiMny==0.0 ? 0.0 : fabs(m_dMuy/dPhiMny));
    double dRatMz = (dPhiMnz==0.0 ? 0.0 : fabs(m_dMuz/dPhiMnz));
    double dPhiMn = sqrt(pow(dPhiMny,2)+pow(dPhiMnz,2));
    double dRatMn = (dPhiMn==0.0 ? 0.0 : fabs(dMu/dPhiMn));
    double dRatM = max(dRatMn, max(dRatMy, dRatMz));
	  dRatio = max(dRatP,dRatM);

		ResData.nType = 2;
		ResData.dNARotate = PmBalance.dAngleb;
		ResData.dCb       = min(PmBalance.dDmax, (1.0-dXn)*PmBalance.dDmax); 
		ResData.dXn       = dXn;
		ResData.dPn       = dPn;
		ResData.dMny      = dMny;
		ResData.dMnz      = dMnz;
		ResData.desiMax   = desi;
		if(pResData != NULL) (*pResData) = ResData;
		return TRUE;
	}
	if(pResData != NULL) (*pResData) = ResData;
	return FALSE;//0.0
}

double CCalcPmcv::Calc_Rota(double dInRota/*-1.0*/)
{
	double dRota = dInRota;
	// dRota Unit is Radian.
	if(dRota < 0.0)	// Calculate Real Rotate Angle.
	{
		dRota = 0.0;
		
		double dez = (fabs(m_dPu) < cUMDRC_Zero ? fabs(m_dMuy) : fabs(m_dMuy/m_dPu));
		double dey = (fabs(m_dPu) < cUMDRC_Zero ? fabs(m_dMuz) : fabs(m_dMuz/m_dPu));
		if(dey > 0.0 && dez > 0.0)	dRota	= atan(dey/dez);
		else if(dey==0.0)						dRota = 0.0;					// y-y axis ( 0 Degree).
		else if(dez==0.0)						dRota = CMathFunc::m_pi/2.0;	// z-z axis (90 Degree).
	}	
	return dRota;
}

double CCalcPmcv::Calc_Dmax(double dInRota/*-1.0*/)
{
	double dDmax=0.0;
	// Get Required Data from MyDB.
	int iSectShp = m_iSectShp;
  double dH,dB;
  dH = m_dSize[0];
	dB = m_dSize[1];
  
	if(iSectShp==DGN_SECT_SHAPE_INDEX_REG_P || iSectShp==DGN_SECT_SHAPE_INDEX_REG_SR)	dB = dH;
	double dRota = Calc_Rota(dInRota);
	dDmax = dB*sin(dRota) + dH*cos(dRota);
	
	return dDmax;
}


double CCalcPmcv::Calc_Ratio(double dPu, double dMuy, double dMuz, double dPhiPnmax, T_PMCV_2D& PmData, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, BOOL bPureMonemt)
{
	T_PM_SECT_RES ResData;  ResData.Init();
	return Calc_Ratio(dPu, dMuy, dMuz, dPhiPnmax, PmData, dPhi, dPhiPn, dPhiMny, dPhiMnz, ResData, bPureMonemt);
}

double CCalcPmcv::Calc_Ratio(double dPu, double dMuy, double dMuz, double dPhiPnmax, T_PMCV_2D& PmData, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, T_PM_SECT_RES& ResData, BOOL bPureMonemt)
{
	ResData.Init();
	ResData.nType = 2;
	ResData.dNARotate = PmData.dRotate;

  int iCode = m_iDgnCode;

	double dRatio=0.0;	// max(RatP, RatM).

	double dEcc1 = ( 1.01) * cUMDRC_Upon;
	double dEcc2 = (-1.01) * cUMDRC_Upon;
	double dPhi1=0.0, dPhi2=0.0;
	double dpPn1=0.0, dpPn2=0.0;
	double dpMny1=0.0, dpMny2=0.0;
  double dpMnz1=0.0, dpMnz2=0.0;	
  double dMu = sqrt(pow(dMuy,2)+pow(dMuz,2));

	T_PM_UNIT	PmUnit, PmUnit1, PmUnit2;

	if(fabs(dPu) > cUMDRC_Zero && fabs(dMu) < cUMDRC_Zero)	// Pure Tension or Compression.
	{
		double dPhicmax, dPhitmax;
    double dpPcmax=0.0, dpPtmax=0.0, dpMnmax=0.0;
    double dpPcMny=0.0, dpPtMny=0.0;
    double dpPcMnz=0.0, dpPtMnz=0.0;
	  // 0 Deg.
		int nSize = PmData.arPmUnit.GetSize();
	  for(int i=0 ; i<nSize ; i++)
	  {
		  PmUnit = PmData.arPmUnit.GetAt(i);
		  double dpPn = PmUnit.dPhi * PmUnit.dPn;
      double dMn = sqrt(PmUnit.dMny*PmUnit.dMny+PmUnit.dMnz*PmUnit.dMnz);
		  double dpMn = PmUnit.dPhi * dMn;
		  if(dpMn >= dpMnmax)	dpMnmax = dpMn;
			//if(i==0)				{dPhicmax = PmUnit.dPhi;  dpPcmax = dpPn; dpPcMny = PmUnit.dPhi*PmUnit.dMny; dpPcMnz = PmUnit.dPhi*PmUnit.dMnz;}
      if(PmUnit.dXn==0.0)				{dPhicmax = PmUnit.dPhi;  dpPcmax = dpPn; dpPcMny = PmUnit.dPhi*PmUnit.dMny; dpPcMnz = PmUnit.dPhi*PmUnit.dMnz;}
			if(i==nSize-1)	{dPhitmax = PmUnit.dPhi;  dpPtmax = dpPn; dpPtMny = PmUnit.dPhi*PmUnit.dMny; dpPtMnz = PmUnit.dPhi*PmUnit.dMnz;}
	  }
    //Edit Code By RSH 2002.08.24
	  dpPcmax *= Get_MaxPnReductionFactor();	// Assume Ties.

		dPhi    = (dPu > 0.0 ? dPhicmax : dPhitmax);
    if(!bPureMonemt)
    {
		  dPhiMny = 0.0;
      dPhiMnz = 0.0;
    }
    else if(dPu > 0.0)
    {
      dPhiMny = dpPcMny;
      dPhiMnz = dpPcMnz;
    }
    else
    {
      dPhiMny = dpPtMny;
      dPhiMnz = dpPtMnz;
    }
		dPhiPn = (dPu > 0.0 ? dpPcmax : dpPtmax);
		dRatio = (dPhiPn==0.0 ? 0.0 : fabs(dPu/dPhiPn));

		if(nSize > 0) PmUnit = (m_dPu > 0.0) ? PmData.arPmUnit.GetAt(0) : PmData.arPmUnit.GetAt(nSize-1);
		else          PmUnit.Init();
		ResData.dCb       = (m_dPu > 0.0) ? PmData.dDmax : 0.0; 
		ResData.dXn       = PmUnit.dXn;
		ResData.dPn       = PmUnit.dPn;
		ResData.dMny      = 0.0; //PmUnit.dMny;
		ResData.dMnz      = 0.0; //PmUnit.dMnz;
		ResData.desiMax   = PmUnit.desi;
	}
  else if(fabs(dPu) < cUMDRC_Zero && fabs(dMu) < cUMDRC_Zero)	// NO applied force.
	{
		dPhi    = 0.0;
		dPhiMny = 0.0;
    dPhiMnz = 0.0;
		dPhiPn = 0.0;
		dRatio = 0.0;

		ResData.nType = 0;
	}
	else
	{
		// Get Eccu.
    double dPuMu = sqrt(pow(dPu,2)+pow(dMuy,2)+pow(dMuz,2));
    double dEccu = (dPuMu==0.0)? 0.0 : dPu/dPuMu;

		int nSize = PmData.arPmUnit.GetSize();
		for(int i=0 ; i<nSize ; i++)
		{
			PmUnit = PmData.arPmUnit.GetAt(i);
			double dpPn = PmUnit.dPhi * PmUnit.dPn;
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
    if(dMagnitude < cUMDRC_Zero)
    {
			dPhi    = dPhi2;
      dPhiPn  = dpPn2;
      dPhiMny = dpMny2;
      dPhiMnz = dpMnz2;
			
			ResData.dCb       = min(PmData.dDmax, (1.0-PmUnit2.dXn)*PmData.dDmax); 
			ResData.dXn       = PmUnit2.dXn;
			ResData.dPn       = PmUnit2.dPn;
			ResData.dMny      = PmUnit2.dMny;
			ResData.dMnz      = PmUnit2.dMnz;
			ResData.desiMax   = PmUnit2.desi;
    }
    else
    {
			double dUPhi = (dPhi1-dPhi2) / dMagnitude;
      double dUp = (dpPn1-dpPn2) / dMagnitude;
      double dUmy = (dpMny1-dpMny2) / dMagnitude;
      double dUmz = (dpMnz1-dpMnz2) / dMagnitude;
      double dAmp = 0.0;
      if(fabs(dUp) < cUMDRC_Zero) ASSERT(0);
      if(fabs(dEccu) < cUMDRC_Zero)
      {
        dAmp = -dpPn2/dUp;
				dPhi    = dPhi2 + dUPhi*dAmp;
        dPhiPn  = dpPn2 + dUp*dAmp;
        dPhiMny = dpMny2 + dUmy*dAmp;
        dPhiMnz = dpMnz2 + dUmz*dAmp;
				
				ResData.dXn       = PmUnit2.dXn  + (PmUnit1.dXn -PmUnit2.dXn )/dMagnitude;
				ResData.dCb       = min(PmData.dDmax, (1.0-ResData.dXn)*PmData.dDmax); 
				ResData.dPn       = PmUnit2.dPn  + (PmUnit1.dPn -PmUnit2.dPn )/dMagnitude;
				ResData.dMny      = PmUnit2.dMny + (PmUnit1.dMny-PmUnit2.dMny)/dMagnitude;
				ResData.dMnz      = PmUnit2.dMnz + (PmUnit1.dMnz-PmUnit2.dMnz)/dMagnitude;
				ResData.desiMax   = PmUnit2.desi + (PmUnit1.desi-PmUnit2.desi)/dMagnitude;
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

        if(fabs(m_dPu) > cUMDRC_Zero)
        {
          if(m_dPu*dPhiPn1 > 0.0 &&
            dPhiPn1 >= min(dpPn1,dpPn2) && dPhiPn1 <= max(dpPn1,dpPn2) &&
            dPhiMny1 >= min(dpMny1,dpMny2) && dPhiMny1 <= max(dpMny1,dpMny2) &&
            dPhiMnz1 >= min(dpMnz1,dpMnz2) && dPhiMnz1 <= max(dpMnz1,dpMnz2))
          {
						dPhi    = dPhi_1;
            dPhiPn  = dPhiPn1;
            dPhiMny = dPhiMny1;
            dPhiMnz = dPhiMnz1;
						
						ResData.dCb       = min(PmData.dDmax, (1.0-PmUnit1.dXn)*PmData.dDmax); 
						ResData.dXn       = PmUnit1.dXn;
						ResData.dPn       = PmUnit1.dPn;
						ResData.dMny      = PmUnit1.dMny;
						ResData.dMnz      = PmUnit1.dMnz;
						ResData.desiMax   = PmUnit1.desi;
          }
          else
          {
						dPhi    = dPhi_2;
            dPhiPn  = dPhiPn2;
            dPhiMny = dPhiMny2;
            dPhiMnz = dPhiMnz2;
						
						ResData.dCb       = min(PmData.dDmax, (1.0-PmUnit2.dXn)*PmData.dDmax); 
						ResData.dXn       = PmUnit2.dXn;
						ResData.dPn       = PmUnit2.dPn;
						ResData.dMny      = PmUnit2.dMny;
						ResData.dMnz      = PmUnit2.dMnz;
						ResData.desiMax   = PmUnit2.desi;
          }
        }
        else
        {
          if(dPhiMny1 >= min(dpMny1,dpMny2) && dPhiMny1 <= max(dpMny1,dpMny2) &&
            dPhiMnz1 >= min(dpMnz1,dpMnz2) && dPhiMnz1 <= max(dpMnz1,dpMnz2))
          {
						dPhi    = dPhi_1;
            dPhiPn  = dPhiPn1;
            dPhiMny = dPhiMny1;
            dPhiMnz = dPhiMnz1;
						
						ResData.dCb       = min(PmData.dDmax, (1.0-PmUnit1.dXn)*PmData.dDmax); 
						ResData.dXn       = PmUnit1.dXn;
						ResData.dPn       = PmUnit1.dPn;
						ResData.dMny      = PmUnit1.dMny;
						ResData.dMnz      = PmUnit1.dMnz;
						ResData.desiMax   = PmUnit1.desi;
          }
          else
          {
						dPhi    = dPhi_2;
            dPhiPn  = dPhiPn2;
            dPhiMny = dPhiMny2;
            dPhiMnz = dPhiMnz2;
						
						ResData.dCb       = min(PmData.dDmax, (1.0-PmUnit2.dXn)*PmData.dDmax); 
						ResData.dXn       = PmUnit2.dXn;
						ResData.dPn       = PmUnit2.dPn;
						ResData.dMny      = PmUnit2.dMny;
						ResData.dMnz      = PmUnit2.dMnz;
						ResData.desiMax   = PmUnit2.desi;
          }        
        }
      }
    }

		// Compare pPn with pPnmax.
		double dpPnmax = dPhiPnmax;

		if(dPhiPn > dpPnmax)	
		{
			double dRatPnmax = dPhiPn==0.0 ? 0.0 : dpPnmax/dPhiPn;
			dPhiMny *= dRatPnmax;
			dPhiMnz *= dRatPnmax;
			dPhiPn *= dRatPnmax;
		}
    double dpMn = sqrt(pow(dPhiMny,2)+pow(dPhiMnz,2));
		// Get Ratios.
		double dRatP = (dPhiPn==0.0 ? 0.0 : fabs(dPu/dPhiPn));
		double dRatMy = (dPhiMny==0.0 ? 0.0 : fabs(dMuy/dPhiMny));
    double dRatMz = (dPhiMnz==0.0 ? 0.0 : fabs(dMuz/dPhiMnz));
    double dRatMn = (dpMn==0.0 ? 0.0 : fabs(dMu/dpMn));
    double dRatM = max(dRatMn, max(dRatMy, dRatMz));
		dRatio = max(dRatP,dRatM);
	}
	return dRatio;
}

double CCalcPmcv::Calc_Ratio_ToPu(double dPu, double dMuy, double dMuz, double dPhiPnmax, T_PMCV_2D& PmData, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, T_PM_SECT_RES& ResData)
{
	ResData.Init();
	dPhi = dPhiPn = dPhiMny = dPhiMnz = 0.0;
	ResData.nType = 0;
	ResData.dNARotate = PmData.dRotate;

  int iCode = m_iDgnCode;

	double dRatio=0.0;	// max(RatP, RatM).

  double dMu = sqrt(pow(dMuy,2)+pow(dMuz,2));
	
	T_PM_UNIT	PmUnit, PmUnit1, PmUnit2;
	int nSize = PmData.arPmUnit.GetSize();
	for(int i=0 ; i<nSize-1 ; i++)
	{
		PmUnit1 = PmData.arPmUnit.GetAt(i);
		PmUnit2 = PmData.arPmUnit.GetAt(i+1);	

		if(PmUnit1.dPn < PmUnit2.dPn)
		{
			PmUnit = PmUnit1;
			PmUnit1 = PmUnit2;
			PmUnit2 = PmUnit;
		}
		
		double dpPn1 = PmUnit1.dPhi * PmUnit1.dPn;
		double dpPn2 = PmUnit2.dPhi * PmUnit2.dPn;
		
		if(dPu <= dpPn1 && dPu >= dpPn2)
		{
			double dMagnitude = (dpPn1-dpPn2)==0.0 ? 1.0 : (dPu-dpPn2)/(dpPn1-dpPn2);
			
			dPhi    = dMagnitude*PmUnit1.dPhi + (1.0-dMagnitude)*PmUnit2.dPhi;
			dPhiPn  = dPhi * (dMagnitude*PmUnit1.dPn  + (1.0-dMagnitude)*PmUnit2.dPn );
			dPhiMny = dPhi * (dMagnitude*PmUnit1.dMny + (1.0-dMagnitude)*PmUnit2.dMny);
			dPhiMnz = dPhi * (dMagnitude*PmUnit1.dMnz + (1.0-dMagnitude)*PmUnit2.dMnz);
			
			ResData.dXn       = (dMagnitude*PmUnit1.dXn  + (1.0-dMagnitude)*PmUnit2.dXn );
			ResData.dCb       = min(PmData.dDmax, (1.0-ResData.dXn)*PmData.dDmax); 
			ResData.dPn       = (dMagnitude*PmUnit1.dPn  + (1.0-dMagnitude)*PmUnit2.dPn );
			ResData.dMny      = (dMagnitude*PmUnit1.dMny + (1.0-dMagnitude)*PmUnit2.dMny);
			ResData.dMnz      = (dMagnitude*PmUnit1.dMnz + (1.0-dMagnitude)*PmUnit2.dMnz);
			ResData.desiMax   = (dMagnitude*PmUnit1.desi + (1.0-dMagnitude)*PmUnit2.desi);

			// Compare pPn with pPnmax.
			double dpPnmax = dPhiPnmax;

			if(dPhiPn > dpPnmax)	
			{
				dPhiMny = 0.0;
				dPhiMnz = 0.0;
				dPhiPn  = dpPnmax;
				ResData.nType = 0;
			}
			else 
				ResData.nType = 2;
			double dpMn = sqrt(pow(dPhiMny,2)+pow(dPhiMnz,2));
			// Get Ratios.
			double dRatP  = (dPhiPn==0.0  ? 0.0 : fabs(dPu/dPhiPn));
			double dRatMy = (dPhiMny==0.0 ? 0.0 : fabs(dMuy/dPhiMny));
			double dRatMz = (dPhiMnz==0.0 ? 0.0 : fabs(dMuz/dPhiMnz));
			double dRatMn = (dpMn==0.0    ? 0.0 : fabs(dMu/dpMn));
			double dRatM = max(dRatMn, max(dRatMy, dRatMz));
			dRatio = max(dRatP,dRatM);
		}
	}

	return dRatio;
}

BOOL CCalcPmcv::Calc_PnMn_Gen(double dXn, double dInRota, double dAsThk, T_PM_SECT_RES& ResData, BOOL bRbarDetail)
{	
  CMdgnSectTool SectTool;
	CPolyMaker PolyMaker;

	double dRotation = dInRota*CMathFunc::m_trang;
	double dCenter_y  = m_dYbar;
	double dCenter_z  = m_dZbar;
	double dYmax, dYmin, dZmax, dZmin;

	BOOL bFixDgnRebarOnly = (m_bDesign && dAsThk<=0.0);

	
	/////////////////////////////////////////////
	// 단면의 회전정보 생성
	_UMD_RC_CON_PART_LIST arConcUnit; arConcUnit.Init();
	_UMD_RC_RBAR_LIST arRbarUnit;     arRbarUnit.Init();
	if(!SectTool.Get_RotateRebarData(dCenter_y, dCenter_z, dRotation, m_RbarData, arRbarUnit)) return FALSE;
	
	if(m_iTypeStressStrain==3)
	{ if(!SectTool.Get_RotateConPolyData(dCenter_y, dCenter_z, dRotation, m_ConPolyData, dYmax, dYmin, dZmax, dZmin, &PolyMaker)) return FALSE; }
	else
	{ if(!SectTool.Get_RotateConUnitPartList(2, dCenter_y, dCenter_z, dRotation, m_iSliceNum, m_ConPolyData, arConcUnit)) return FALSE; }

	dYmax = arConcUnit.dYmax;
	dYmin = arConcUnit.dYmin;
	dZmax = arConcUnit.dZmax;
	dZmin = arConcUnit.dZmin;

	double dDmax = dZmax - dZmin;      
  double dDeff = Calc_Deff_Gen(dDmax, dZmin, arRbarUnit, bFixDgnRebarOnly);	
	
	if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
	{ if(!Calc_PnMn_Gen_RectStress(dXn, dInRota, dAsThk, ResData, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit, bRbarDetail))	return FALSE;	}
	else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
	{	if(!Calc_PnMn_Gen_SliceMethod(dXn, dInRota, dAsThk, ResData, arConcUnit, arRbarUnit, bRbarDetail))	return FALSE;	}
	
	return TRUE;
}

BOOL CCalcPmcv::Is_PureComp(double dXn)	
{
	if(m_iSectShp == DGN_SECT_SHAPE_INDEX_REG_GEN)
		return (fabs(dXn) < cUMDRC_Zero ? TRUE : FALSE);
	else
		return (dXn       < 0.0         ? TRUE : FALSE);
}	// TRUE (dXn=0).
//Check
BOOL CCalcPmcv::Is_PureTens(double dXn)	{return (fabs(dXn-1.0) < cUMDRC_Zero ? TRUE : FALSE);}	// TRUE (dXn=1).



double CCalcPmcv::Calc_MaxDistTensColmRbar(double dXn, double dRota, double dDmax, BOOL bFixDgnRebarOnly)
{
	// Calculate maximum distance from Rebar to neutral axis.
	double dMaxDist=0.0;
	BOOL bDgn = m_bDesign;
	
	_UMD_RC_RBAR_UNIT	RbarUnit;
	if(fabs(dRota-0.0) < cUMDRC_Zero)	// 0 Deg.
	{		
		int iIndexR=0;
		POSITION PosR = m_RbarData.arRbarUnit.GetStartPosition();
		while(PosR)
		{
			RbarUnit.Init();
			m_RbarData.arRbarUnit.GetNextAssoc(PosR,iIndexR,RbarUnit);
			if(bFixDgnRebarOnly)
			{ if(bDgn && RbarUnit.iDgnType==0) continue; }
			// Check if compression or not.
			double dDz1	= dDmax*dXn;
			if(dDz1 > RbarUnit.dyz[1])	// Only Tension.
			{
				double dDist = dDz1 - RbarUnit.dyz[1];
				if(dDist > dMaxDist)	dMaxDist = dDist;
			}
		}
	}
	else if(fabs(dRota-CMathFunc::m_pi/2.0) < cUMDRC_Zero)	// 90 Deg.
	{		
		int iIndexR=0;
		POSITION PosR = m_RbarData.arRbarUnit.GetStartPosition();
		while(PosR)
		{
			RbarUnit.Init();
			m_RbarData.arRbarUnit.GetNextAssoc(PosR,iIndexR,RbarUnit);
			if(bFixDgnRebarOnly)
			{ if(bDgn && RbarUnit.iDgnType==0) continue; }
			// Check if compression or not.
			double dDy1	= dDmax*dXn;
			if(dDy1 > RbarUnit.dyz[0])	// Only Tension.
			{
				double dDist = dDy1 - RbarUnit.dyz[0];
				if(dDist > dMaxDist)	dMaxDist = dDist;
			}
		}
	}
	else
	{		
		int iIndexR=0;
		POSITION PosR = m_RbarData.arRbarUnit.GetStartPosition();
		while(PosR)
		{
			RbarUnit.Init();
			m_RbarData.arRbarUnit.GetNextAssoc(PosR,iIndexR,RbarUnit);
			if(bFixDgnRebarOnly)
			{ if(bDgn && RbarUnit.iDgnType==0) continue; }
			// Check if compression or not.
			double dDy1		 = dDmax*dXn / sin(dRota);
			double dDz1		 = dDmax*dXn / cos(dRota);
			double dCoordZ = (-1)*tan(dRota)*RbarUnit.dyz[0] + dDz1;	// by Neutral axis.
			if(dCoordZ > RbarUnit.dyz[1])	// Only Tension.
			{
        if(dDy1 == 0.0 && dDz1 == 0.0)
        {
			    dDy1 = cUMDRC_Zero / sin(dRota);
			    dDz1 = cUMDRC_Zero / cos(dRota);
        }

				double dLineI[3]={dDy1,  0.0, 0.0};
				double dLineJ[3]={ 0.0, dDz1, 0.0};
				double dPoint[3]={RbarUnit.dyz[0], RbarUnit.dyz[1], 0.0};
				double dDist = CMathFunc::mathDistanceFromIntersectPointToLine(dLineI,dLineJ,dPoint);
				if(dDist > dMaxDist)	dMaxDist = dDist;
			}
		}
	}
	return dMaxDist;
}

double CCalcPmcv::Calc_MaxDistTensColmRbar2(double dDy, double dDz, double dRota, BOOL bFixDgnRebarOnly)
{  
  // Calculate maximum distance from Rebar to neutral axis.
	double dMaxDist=0.0;
  // Get Required Data from MyDB.
  double dYbar = m_dYbar;
	double dZbar = m_dZbar;
  
	_UMD_RC_RBAR_UNIT	RbarUnit;
	BOOL bDgn = m_bDesign;

	if(fabs(dRota-0.0) < cUMDRC_Zero)	// 0 Deg.
	{
		int iIndexR=0;
		POSITION PosR = m_RbarData.arRbarUnit.GetStartPosition();
		while(PosR)
		{
			RbarUnit.Init();
			m_RbarData.arRbarUnit.GetNextAssoc(PosR,iIndexR,RbarUnit);
			if(bFixDgnRebarOnly)
			{ if(bDgn && RbarUnit.iDgnType==0) continue; }
			// Check if compression or not.
			double dDz1	= dDz;
			if(dDz1 > RbarUnit.dyz[1])	// Only Tension.
			{
				double dDist = dDz1 - RbarUnit.dyz[1];
				if(dDist > dMaxDist)	dMaxDist = dDist;
			}
		}
	}
	else if(fabs(dRota-CMathFunc::m_pi/2.0) < cUMDRC_Zero)	// 90 Deg.
	{
		int iIndexR=0;
		POSITION PosR = m_RbarData.arRbarUnit.GetStartPosition();
		while(PosR)
		{
			RbarUnit.Init();
			m_RbarData.arRbarUnit.GetNextAssoc(PosR,iIndexR,RbarUnit);
			if(bFixDgnRebarOnly)
			{ if(bDgn && RbarUnit.iDgnType==0) continue; }
			// Check if compression or not.
			double dDy1	= dDy;
			if(dDy1 > RbarUnit.dyz[0])	// Only Tension.
			{
				double dDist = dDy1 - RbarUnit.dyz[0];
				if(dDist > dMaxDist)	dMaxDist = dDist;
			}
		}
	}
	else
	{		
    int iIndexR=0;
		POSITION PosR = m_RbarData.arRbarUnit.GetStartPosition();
		while(PosR)
		{
			RbarUnit.Init();
			m_RbarData.arRbarUnit.GetNextAssoc(PosR,iIndexR,RbarUnit);
			if(bFixDgnRebarOnly)
			{ if(m_bDesign && RbarUnit.iDgnType==0) continue; }
			// Check if compression or not.
			double dDy1		 = dDy;
			double dDz1		 = dDz;
			double dCoordZ = (-1)*tan(dRota)*RbarUnit.dyz[0] + dDz1;	// by Neutral axis.
			if(dCoordZ > RbarUnit.dyz[1])	// Only Tension.
			{
        if(dDy1 == 0.0 && dDz1 == 0.0)
        {
			    dDy1 = cUMDRC_Zero / sin(dRota);
			    dDz1 = cUMDRC_Zero / cos(dRota);
        }

				double dLineI[3]={dDy1,  0.0, 0.0};
				double dLineJ[3]={ 0.0, dDz1, 0.0};
				double dPoint[3]={RbarUnit.dyz[0], RbarUnit.dyz[1], 0.0};
				double dDist = CMathFunc::mathDistanceFromIntersectPointToLine(dLineI,dLineJ,dPoint);
				if(dDist > dMaxDist)	dMaxDist = dDist;
			}
		}
	}
	return dMaxDist;
}
BOOL CCalcPmcv::Get_PMCurve_Gen(double dInRota, T_PMCV_2D& PmcvData, T_PM_SECT_RES& CbPmData, double& dPhiPnmax, double dAsThk, BOOL bRbarDetail)
{
  PmcvData.Init();
	CbPmData.Init();

	CMdgnSectTool SectTool;
	CPolyMaker PolyMaker;

	double dRotation = dInRota*CMathFunc::m_trang;
	double dCenter_y  = m_dYbar;
	double dCenter_z  = m_dZbar;
	double dYmax, dYmin, dZmax, dZmin;

	BOOL bFixDgnRebarOnly = (m_bDesign && dAsThk<=0.0);
	
	/////////////////////////////////////////////
	// 단면의 회전정보 생성
	_UMD_RC_CON_PART_LIST arConcUnit; arConcUnit.Init();
	_UMD_RC_RBAR_LIST arRbarUnit;     arRbarUnit.Init();
	if(!SectTool.Get_RotateRebarData(dCenter_y, dCenter_z, dRotation, m_RbarData, arRbarUnit)) return FALSE;
	
	if(m_iTypeStressStrain==3)
	{ if(!SectTool.Get_RotateConPolyData(dCenter_y, dCenter_z, dRotation, m_ConPolyData, dYmax, dYmin, dZmax, dZmin, &PolyMaker)) return FALSE; }
	else
	{ if(!SectTool.Get_RotateConUnitPartList(2, dCenter_y, dCenter_z, dRotation, m_iSliceNum, m_ConPolyData, arConcUnit)) return FALSE; }

	dYmax = arConcUnit.dYmax;
	dYmin = arConcUnit.dYmin;
	dZmax = arConcUnit.dZmax;
	dZmin = arConcUnit.dZmin;

	double dDmax = dZmax - dZmin;      
  double dDeff = Calc_Deff_Gen(dDmax, dZmin, arRbarUnit, bFixDgnRebarOnly);	

	dPhiPnmax=0.0;

  double dFyr = m_dFyr;
	T_PM_UNIT	PmUnit;

	double dPn=0.0, dMn=0.0;
	double dPb=0.0, dMb=0.0;
	double dPcmax=0.0, dPtmax=0.0;
	int    iPcmaxID=1, iPtmaxID=1;
  double dXn = 0.0;
	int iDivNum = m_iPmDivNum;
	PmcvData.arPmUnit.SetSize(iDivNum+1);
	T_PM_SECT_RES ResData;
  for(int i=0; i<=iDivNum; i++)
	{
		// IRC_21_00일 경우에는 재정의 필요
    if(i==0)	dXn = -0.5;
    else
    {
		  double dX	 = (double)(i-1)/(iDivNum-1);
		  dXn = 2.0202*pow(dX,3) - 3.13131*pow(dX,2) + 2.11111*dX;
    }

		if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
		{ if(!Calc_PnMn_Gen_RectStress(dXn, dInRota, dAsThk, ResData, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit))	return FALSE;	}
		else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
		{	if(!Calc_PnMn_Gen_SliceMethod(dXn, dInRota, dAsThk, ResData, arConcUnit, arRbarUnit))	return FALSE;	}
		
		PmUnit.Init();
		PmUnit.dXn	= dXn;
		PmUnit.dPn	= ResData.dPn;
		PmUnit.dMny	= ResData.dMny;
    PmUnit.dMnz	= ResData.dMnz;
    PmUnit.desi = ResData.desiMax;
		dPn = PmUnit.dPn;
    dMn = sqrt(ResData.dMny*ResData.dMny + ResData.dMnz*ResData.dMnz);
		// Check Pb, Mb (Balanced condition).
		if(dMn > dMb)	{dPb = dPn;	dMb = dMn;}
		// Calculate Pcmax, Ptmax.
		if(dPcmax < ResData.dPn)	{dPcmax = ResData.dPn; iPcmaxID=i;}
		if(dPtmax > ResData.dPn)	{dPtmax = ResData.dPn; iPtmaxID=i;}
		PmcvData.arPmUnit.SetAt(i, PmUnit);
	}
  // Set Phi, Mn(Not calculated).
	for(int i=0; i<=iDivNum; i++)
	{
		PmUnit = PmcvData.arPmUnit.GetAt(i);
		PmUnit.dPhi = Calc_Phi(PmUnit.dPn, dPb, PmUnit.desi);
		PmcvData.arPmUnit.SetAt(i, PmUnit);
	}
	
	PmcvData.dRotate = dInRota;
	PmcvData.dDmax   = dDmax;
	PmcvData.dDeff   = dDeff;
	
	if(m_bFindMmax)
	{
		BOOL bStop = FALSE;
		double dXnTol	 = 1.0-7*m_dTolM0;	
		double dXn0Tol = 1.0-1*m_dTolM0;				
		double dMn_y=0.0, dMn_z=0.0;
		double dXnPrev=0.0, dXnNext=1.0;
		double dMnPrev=0.0, dMnNext=0.0;
		
		int iIterCount=0;
		// Bi-section Method.
		while(!bStop)
		{
			// 1st Pmcv.
			dXn = (dXnPrev+dXnNext) / 2.0;
			
			if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
			{ if(!Calc_PnMn_Gen_RectStress(dXn, dInRota, dAsThk, CbPmData, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit))	return FALSE;	}
			else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
			{	if(!Calc_PnMn_Gen_SliceMethod(dXn, dInRota, dAsThk, CbPmData, arConcUnit, arRbarUnit))	return FALSE;	}

			dPn = CbPmData.dPn;   dMn_y = CbPmData.dMny;   dMn_z = CbPmData.dMnz;
			dMn = sqrt(dMn_y*dMn_y+dMn_z*dMn_z);
			// Calculate Only if Balanced Condition required.
			
			// Change by ZINU.('03.12.18). Calculate by 3 Points.
			//  M
			//  |--
			//  |   3rd
			//  |     1st
			//  |       2nd
			//  |        |
			//  +------+------ P
			//  |   __/
			//  |--'
			// 2nd Pmcv for Xn Direction.
			double dXn2 = dXn / dXnTol;
			T_PM_SECT_RES resultData2;
			double dPn2=0.0, dMn2=0.0, dMn2_y=0.0, dMn2_z=0.0;

			if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
			{ if(!Calc_PnMn_Gen_RectStress(dXn2, dInRota, dAsThk, resultData2, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit))	return FALSE;	}
			else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
			{	if(!Calc_PnMn_Gen_SliceMethod(dXn2, dInRota, dAsThk, resultData2, arConcUnit, arRbarUnit))	return FALSE;	}
			
			dPn2 = resultData2.dPn;   dMn2_y = resultData2.dMny;   dMn2_z = resultData2.dMnz;
			dMn2 = sqrt(dMn2_y*dMn2_y+dMn2_z*dMn2_z);
			// 3rd Pmcv for Xn Direction.
			double dXn3 = dXn * dXnTol;
			T_PM_SECT_RES resultData3;
			double dPn3=0.0, dMn3=0.0, dMn3_y=0.0, dMn3_z=0.0;

			if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
			{ if(!Calc_PnMn_Gen_RectStress(dXn3, dInRota, dAsThk, resultData3, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit))	return FALSE;	}
			else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
			{	if(!Calc_PnMn_Gen_SliceMethod(dXn3, dInRota, dAsThk, resultData3, arConcUnit, arRbarUnit))	return FALSE;	}

			dPn3 = resultData3.dPn;   dMn3_y = resultData3.dMny;   dMn3_z = resultData3.dMnz;
			dMn3 = sqrt(dMn3_y*dMn3_y+dMn3_z*dMn3_z);
			// Compare 1st' with 2nd' and 3rd'.
			if(dMn2 < dMn && dMn3 < dMn)	bStop = TRUE;	// Satisfy.
			else if(dMn2 < dMn && dMn3 > dMn)	dXnNext = dXn;	// Move Compressive Region.       
			else if(dMn2 > dMn && dMn3 < dMn)	dXnPrev = dXn;	// Move Tensile Region.
			else
			{
				if(dMn3 > dMn2)	dXnNext = dXn;	// Move Compressive Region.
				else						dXnPrev = dXn;	// Move Tensile Region.
			}
			iIterCount++;
			if(iIterCount > m_iIterNum)	bStop = TRUE;	// Not Satisfy.
			
			if(bStop)
			{
				// SHIN : dXn의 값이 꼭 최대값을 가지는 것은 아니므로 3개의 평균값을 취함
				//        단 Mn값을 평균을 취하면 최대값보다 작은 값을 가지게 되므로 평균값을 취하지 않음
				dXn = (dXn+dXn2+dXn3)/3.0;
				dPn = (dPn+dPn2+dPn3)/3.0;
        if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
		    { if(!Calc_PnMn_Gen_RectStress(dXn, dInRota, dAsThk, CbPmData, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit, bRbarDetail))	return FALSE;	}
		    else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
		    {	if(!Calc_PnMn_Gen_SliceMethod(dXn, dInRota, dAsThk, CbPmData, arConcUnit, arRbarUnit, bRbarDetail))	return FALSE;	}
			}			      
		}
	}
	else
	{			
		double dXb    = Get_BalancedXb(dDmax, dDeff);
		if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
		{ if(!Calc_PnMn_Gen_RectStress(dXb, dInRota, dAsThk, CbPmData, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit, bRbarDetail))	return FALSE;	}
		else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
		{	if(!Calc_PnMn_Gen_SliceMethod(dXb, dInRota, dAsThk, CbPmData, arConcUnit, arRbarUnit, bRbarDetail))	return FALSE;	}
	}
	
  
	PmUnit.Init();
	PmUnit = PmcvData.arPmUnit.GetAt(iPcmaxID);
	dPhiPnmax = Get_MaxPnReductionFactor() * PmUnit.dPhi * PmUnit.dPn;

	// Calc. Pn_emin
	T_PM_UNIT Pmemin, PmUnit1, PmUnit2;	
	double dPn1, dPn2;
	int nSize = PmcvData.arPmUnit.GetSize();
	if(nSize > 2)
	{
		PmUnit1 = PmcvData.arPmUnit.GetAt(0);
		for(int i=1; i<nSize; i++)
		{
			PmUnit2 = PmcvData.arPmUnit.GetAt(i);
			dPn1 = PmUnit1.dPhi*PmUnit1.dPn;
			dPn2 = PmUnit2.dPhi*PmUnit2.dPn;
			if(dPhiPnmax <= max(dPn1,dPn2)+cUMDRC_Zero && dPhiPnmax >= min(dPn1,dPn2)-cUMDRC_Zero)
			{
				double dWeight1, dWeight2;
				dWeight1 = dPn1-dPn2==0.0 ? 1.0 : (dPhiPnmax-dPn2)/(dPn1-dPn2);
				dWeight2 = 1.0 - dWeight1;
				Pmemin.dPhi = dWeight1*PmUnit1.dPhi + dWeight2*PmUnit2.dPhi;
				Pmemin.dPn  = dWeight1*PmUnit1.dPn  + dWeight2*PmUnit2.dPn;
				Pmemin.dMny = dWeight1*PmUnit1.dMny + dWeight2*PmUnit2.dMny;
				Pmemin.dMnz = dWeight1*PmUnit1.dMnz + dWeight2*PmUnit2.dMnz;
				Pmemin.dXn  = dWeight1*PmUnit1.dXn  + dWeight2*PmUnit2.dXn;
				Pmemin.desi = dWeight1*PmUnit1.desi + dWeight2*PmUnit2.desi;
				break;
			}
			PmUnit1 = PmUnit2;
		}
		PmcvData.Pmemin = Pmemin;
	}
	return TRUE;
}

BOOL CCalcPmcv::Get_PMCurve_Gen(int iAxisDir, T_PMCV_2D& PmcvData, T_PM_SECT_RES& CbPmData, T_PM_SECT_RES& P0PmData, 
																													double& dPhiPnmax, double dAsThk, BOOL bRbarDetail)
{	
	PmcvData.Init();
	CbPmData.Init();
	P0PmData.Init();

	int iDivNum = m_iPmDivNum-1;
  if(iDivNum==0)		return FALSE;
	if(iDivNum%2!=0)	return FALSE;	// Always even.

	CMdgnSectTool SectTool;
	CPolyMaker PolyMaker;

	double dRotation = Get_RotationAngle(iAxisDir);
	double dInRot    = dRotation*CMathFunc::m_trrad;
	double dCenter_y  = m_dYbar;
	double dCenter_z  = m_dZbar;
	double dYmax, dYmin, dZmax, dZmin;

	BOOL bFixDgnRebarOnly = (m_bDesign && dAsThk<=0.0);

	/////////////////////////////////////////////
	// 단면의 회전정보가 없으면 값을 채워 넣음
	if(m_iTypeStressStrain==3)
	{
		if(m_RotSectData.nDataType != 1 && m_RotSectData.nDataType != 2)
		{			
			SectTool.Get_RotateSectData(1, m_dYbar, m_dZbar, m_iSliceNum, m_iDivisionNumber_90deg, m_ConPolyData, m_RbarData, m_RotSectData);
		}
		SectTool.Get_RotateConPolyData(dCenter_y, dCenter_z, dRotation, m_ConPolyData, dYmax, dYmin, dZmax, dZmin, &PolyMaker);
	}
	else
	{
		if(m_RotSectData.nDataType != 2)
		{
			SectTool.Get_RotateSectData(2, m_dYbar, m_dZbar, m_iSliceNum, m_iDivisionNumber_90deg, m_ConPolyData, m_RbarData, m_RotSectData);
		}
	}
	_UMD_RC_CON_PART_LIST arConcUnit;
	_UMD_RC_RBAR_LIST arRbarUnit;

	if(!Cng_AxisDir_Gen(iAxisDir, arConcUnit, arRbarUnit))	return FALSE;
	dYmax = arConcUnit.dYmax;
	dYmin = arConcUnit.dYmin;
	dZmax = arConcUnit.dZmax;
	dZmin = arConcUnit.dZmin;

	double dDmax = dZmax - dZmin;      
  double dDeff = Calc_Deff_Gen(dDmax, dZmin, arRbarUnit, bFixDgnRebarOnly);
	
	int i=0;
  
	/////////////////////////////////
	// Search Balanced Condition.	  
  T_PM_SECT_RES resultData;
    
	double dXb, dPb;
  double dXnTol	 = 1.0-7*m_dTolM0;	
	double dXn0Tol = 1.0-1*m_dTolM0;	
	
	double dXn=0.0, dPn=0.0, dMn=0.0, dMn_y=0.0, dMn_z=0.0;
	double dXnPrev=0.0, dXnNext=1.0;
	double dMnPrev=0.0, dMnNext=0.0;
	BOOL bStop = FALSE;// 무조건 최대 모멘트 지점을 구함(m_bCalcPmcv ? FALSE : TRUE);
	
  int iIterCount=0;
  if(m_bFindMmax)
  {		  
		// Bi-section Method.
		while(!bStop)
		{
			// 1st Pmcv.
			dXn = (dXnPrev+dXnNext) / 2.0;
			
			if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
			{ if(!Calc_PnMn_Gen_RectStress(dXn, dInRot, dAsThk, resultData, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit))	return FALSE;	}
			else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
			{	if(!Calc_PnMn_Gen_SliceMethod(dXn, dInRot, dAsThk, resultData, arConcUnit, arRbarUnit))	return FALSE;	}
      
			dPn = resultData.dPn;   dMn_y = resultData.dMny;   dMn_z = resultData.dMnz;
      dMn = sqrt(dMn_y*dMn_y+dMn_z*dMn_z);
			// Calculate Only if Balanced Condition required.
			
			// Change by ZINU.('03.12.18). Calculate by 3 Points.
			//  M
			//  |--
			//  |   3rd
			//  |     1st
			//  |       2nd
			//  |        |
			//  +------+------ P
			//  |   __/
			//  |--'
			// 2nd Pmcv for Xn Direction.
			double dXn2 = dXn / dXnTol;
      T_PM_SECT_RES resultData2;
			double dPn2=0.0, dMn2=0.0, dMn2_y=0.0, dMn2_z=0.0;

			if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
			{ if(!Calc_PnMn_Gen_RectStress(dXn2, dInRot, dAsThk, resultData2, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit))	return FALSE;	}
			else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
			{	if(!Calc_PnMn_Gen_SliceMethod(dXn2, dInRot, dAsThk, resultData2, arConcUnit, arRbarUnit))	return FALSE;	}
      
			dPn2 = resultData2.dPn;   dMn2_y = resultData2.dMny;   dMn2_z = resultData2.dMnz;
      dMn2 = sqrt(dMn2_y*dMn2_y+dMn2_z*dMn2_z);
			// 3rd Pmcv for Xn Direction.
			double dXn3 = dXn * dXnTol;
      T_PM_SECT_RES resultData3;
			double dPn3=0.0, dMn3=0.0, dMn3_y=0.0, dMn3_z=0.0;

			if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
			{ if(!Calc_PnMn_Gen_RectStress(dXn3, dInRot, dAsThk, resultData3, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit))	return FALSE;	}
			else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
			{	if(!Calc_PnMn_Gen_SliceMethod(dXn3, dInRot, dAsThk, resultData3, arConcUnit, arRbarUnit))	return FALSE;	}
      
      dPn3 = resultData3.dPn;   dMn3_y = resultData3.dMny;   dMn3_z = resultData3.dMnz;
      dMn3 = sqrt(dMn3_y*dMn3_y+dMn3_z*dMn3_z);
			// Compare 1st' with 2nd' and 3rd'.
			if(dMn2 < dMn && dMn3 < dMn)	bStop = TRUE;	// Satisfy.
			else if(dMn2 < dMn && dMn3 > dMn)	dXnNext = dXn;	// Move Compressive Region.       
			else if(dMn2 > dMn && dMn3 < dMn)	dXnPrev = dXn;	// Move Tensile Region.
			else
			{
				if(dMn3 > dMn2)	dXnNext = dXn;	// Move Compressive Region.
				else						dXnPrev = dXn;	// Move Tensile Region.
			}
			iIterCount++;
			if(iIterCount > m_iIterNum)	bStop = TRUE;	// Not Satisfy.
			
      if(bStop)
			{
        // SHIN : dXn의 값이 꼭 최대값을 가지는 것은 아니므로 3개의 평균값을 취함
        //        단 Mn값을 평균을 취하면 최대값보다 작은 값을 가지게 되므로 평균값을 취하지 않음
				dXn = (dXn+dXn2+dXn3)/3.0;
				dPn = (dPn+dPn2+dPn3)/3.0;
        
		    if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
		    { if(!Calc_PnMn_Gen_RectStress(dXn, dInRot, dAsThk, resultData, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit, bRbarDetail))	return FALSE;	}
		    else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
		    {	if(!Calc_PnMn_Gen_SliceMethod(dXn, dInRot, dAsThk, resultData, arConcUnit, arRbarUnit, bRbarDetail))	return FALSE;	}       
			}			      
		}
  }
  else
  {
    dXn = Get_BalancedXb(dDmax, dDeff);

		if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
		{ if(!Calc_PnMn_Gen_RectStress(dXn, dInRot, dAsThk, resultData, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit, bRbarDetail))	return FALSE;	}
		else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
		{	if(!Calc_PnMn_Gen_SliceMethod(dXn, dInRot, dAsThk, resultData, arConcUnit, arRbarUnit, bRbarDetail))	return FALSE;	}        
  }
	// 계산된 값으로 채움
	CbPmData = resultData;
	dPb = resultData.dPn;
	dXb = dXn;

	/////////////////////////////////
	// Calculate M0.
	// Set dMn0 by m_dPu.
	double dP0 = 0.0;	// SHIN ('06.02.15) 원하는 하중값으로 설정할 생각이나 사용하지 않음 defult=0.0
	double dXn0=0.0, den0=0.0, dPn0=0.0, dMn0=0.0, dMn0_y=0.0, dMn0_z=0.0;
	double dXn0Prev=0.0, dXn0Next=1.0;
	bStop = FALSE;
	
  iIterCount=0;
  // Bi-section Method.
	while(!bStop)
	{
		dXn0 = (dXn0Prev+dXn0Next) / 2.0;
		if(fabs(dP0) < cUMDRC_Zero)	// Pure Bending.
		{
			// 2nd Pmcv for Xn Direction.
			double dXn2 = dXn0 / dXn0Tol;
      T_PM_SECT_RES resultData2;
			double dPn2=0.0, dMn2=0.0, dMn2_y=0.0, dMn2_z=0.0;

			if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
			{ if(!Calc_PnMn_Gen_RectStress(dXn2, dInRot, dAsThk, resultData2, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit))	return FALSE;	}
			else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
			{	if(!Calc_PnMn_Gen_SliceMethod(dXn2, dInRot, dAsThk, resultData2, arConcUnit, arRbarUnit))	return FALSE;	}
      
      dPn2 = resultData2.dPn;   dMn2_y = resultData2.dMny;   dMn2_z = resultData2.dMnz;
      dMn2 = sqrt(dMn2_y*dMn2_y+dMn2_z*dMn2_z);

			// 3rd Pmcv for Xn Direction.
			double dXn3 = dXn0 * dXn0Tol;
      T_PM_SECT_RES resultData3;
			double dPn3=0.0, dMn3=0.0, dMn3_y=0.0, dMn3_z=0.0;

			if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
			{ if(!Calc_PnMn_Gen_RectStress(dXn3, dInRot, dAsThk, resultData3, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit))	return FALSE;	}
			else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
			{	if(!Calc_PnMn_Gen_SliceMethod(dXn3, dInRot, dAsThk, resultData3, arConcUnit, arRbarUnit))	return FALSE;	}
 
		  dPn3 = resultData3.dPn;   dMn3_y = resultData3.dMny;   dMn3_z = resultData3.dMnz;
      dMn3 = sqrt(dMn3_y*dMn3_y+dMn3_z*dMn3_z);

			if(dPn2*dPn3 < 0.0)	bStop = TRUE;	// Satisfy.
			else if(dPn2 > 0.0 && dPn3 > 0.0)	dXn0Prev = dXn0;	// Move Tensile Region (Include dPn2=dPn3=0).
			else															dXn0Next = dXn0;	// Move Compressive Region.
			iIterCount++;
			if(iIterCount > m_iIterNum)	bStop = TRUE;	// Not Satisfy.
			if(bStop)
			{
				if(fabs(dPn2+dPn3) < cUMDRC_Zero)
        {
          dMn0_y = (dMn2_y+dMn3_y)/2.0;
          dMn0_z = (dMn2_z+dMn3_z)/2.0;
        }
				else if(fabs(dPn2-dPn3) < cUMDRC_Zero)	
        {
          dMn0_y = (dMn2_y+dMn3_y)/2.0;
          dMn0_z = (dMn2_z+dMn3_z)/2.0;
        }
				else
        {
          dMn0_y = (dMn2_y*dPn3-dMn3_y*dPn2)/(dPn3-dPn2);
          dMn0_z = (dMn2_z*dPn3-dMn3_z*dPn2)/(dPn3-dPn2);
        }
				dXn0 = (dXn2+dXn3)/2.0;
			}
		}
		else	// Axial + Bending.
		{
			if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
			{ if(!Calc_PnMn_Gen_RectStress(dXn0, dInRot, dAsThk, resultData, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit))	return FALSE;	}
			else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
			{	if(!Calc_PnMn_Gen_SliceMethod(dXn0, dInRot, dAsThk, resultData, arConcUnit, arRbarUnit))	return FALSE;	}
 
      dPn0 = resultData.dPn;   dMn0_y = resultData.dMny;   dMn0_z = resultData.dMnz;
			double dTol = fabs((dP0-dPn0)/dP0);
			if(dTol < m_dTolM0)	bStop = TRUE;
			else if(dPn0 > dP0)	dXn0Prev = dXn0;	// Move Tensile Region.
			else if(dPn0 < dP0)	dXn0Next = dXn0;	// Move Compressive Region.
			iIterCount++;
			if(iIterCount > m_iIterNum)	bStop = TRUE;	// Not Satisfy.
		}
	}
	// 계산된 값으로 채움
	P0PmData = resultData;

	/////////////////////////////////
	// Calculate Pn, Mn.   
	T_PM_UNIT	PmUnit;
	PmcvData.arPmUnit.SetSize(iDivNum+1);
	PmcvData.dRotate = dInRot;
	PmcvData.dDmax   = dDmax;
	PmcvData.dDeff   = dDeff;

	int iPcmaxID, iPtmaxID;
	double dPcmax, dPtmax;
	iPcmaxID = iPtmaxID = 0;
	dPcmax = dPtmax = 0.0;
	for(int i=0; i<=iDivNum; i++)
	{
		dXn = Get_Xn(i,iDivNum, dXb);

		if(m_iTypeStressStrain==3)//Whitney(직사각형 응력블럭)
		{ if(!Calc_PnMn_Gen_RectStress(dXn, dInRot, dAsThk, resultData, dYmax, dYmin, dZmax, dZmin, &PolyMaker, arRbarUnit))	return FALSE;	}
		else                      //Parabolic-Descending, Parabolic-Plateau, Whitney(직사각형 응력블럭:SliceMethod)
		{	if(!Calc_PnMn_Gen_SliceMethod(dXn, dInRot, dAsThk, resultData, arConcUnit, arRbarUnit))	return FALSE;	}        

		PmUnit.Init();		
		PmUnit.dXn	= dXn;
		PmUnit.dPn	= resultData.dPn;
		PmUnit.dMny	= resultData.dMny;
    PmUnit.dMnz	= resultData.dMnz;
    PmUnit.desi = resultData.desiMax;
		PmUnit.dPhi = Calc_Phi(PmUnit.dPn, dPb, PmUnit.desi);
    dMn = sqrt(resultData.dMny*resultData.dMny + resultData.dMnz*resultData.dMnz);		
		// Calculate Pcmax, Ptmax.
		if(dPcmax < resultData.dPn)	{dPcmax = resultData.dPn; iPcmaxID=i;}
		if(dPtmax > resultData.dPn)	{dPtmax = resultData.dPn; iPtmaxID=i;}
		PmcvData.arPmUnit.SetAt(i, PmUnit);	  
	}

	dPhiPnmax = Get_MaxPnReductionFactor() * PmUnit.dPhi * dPcmax;	

	// Calc. Pn_emin
	T_PM_UNIT Pmemin, PmUnit1, PmUnit2;	
	double dPn1, dPn2;
	int nSize = PmcvData.arPmUnit.GetSize();
	if(nSize > 2)
	{
		PmUnit1 = PmcvData.arPmUnit.GetAt(0);
		for(int i=1; i<nSize; i++)
		{
			PmUnit2 = PmcvData.arPmUnit.GetAt(i);
			dPn1 = PmUnit1.dPhi*PmUnit1.dPn;
			dPn2 = PmUnit2.dPhi*PmUnit2.dPn;
			if(dPhiPnmax <= max(dPn1,dPn2)+cUMDRC_Zero && dPhiPnmax >= min(dPn1,dPn2)-cUMDRC_Zero)
			{
				double dWeight1, dWeight2;
				dWeight1 = dPn1-dPn2==0.0 ? 1.0 : (dPhiPnmax-dPn2)/(dPn1-dPn2);
				dWeight2 = 1.0 - dWeight1;
				Pmemin.dPhi = dWeight1*PmUnit1.dPhi + dWeight2*PmUnit2.dPhi;
				Pmemin.dPn  = dWeight1*PmUnit1.dPn  + dWeight2*PmUnit2.dPn;
				Pmemin.dMny = dWeight1*PmUnit1.dMny + dWeight2*PmUnit2.dMny;
				Pmemin.dMnz = dWeight1*PmUnit1.dMnz + dWeight2*PmUnit2.dMnz;
				Pmemin.dXn  = dWeight1*PmUnit1.dXn  + dWeight2*PmUnit2.dXn;
				Pmemin.desi = dWeight1*PmUnit1.desi + dWeight2*PmUnit2.desi;
				break;
			}
			PmUnit1 = PmUnit2;
		}
		PmcvData.Pmemin = Pmemin;
	}
  return TRUE;
}

double CCalcPmcv::Get_RbarLenDgn()
{
	double dRbarLen=0.0;
	
	_UMD_RC_RBAR_UNIT RbarUnit;
	int iIndex=0;
	POSITION Pos = m_RbarData.arRbarUnit.GetStartPosition();
	while(Pos)
	{
		RbarUnit.Init();
		m_RbarData.arRbarUnit.GetNextAssoc(Pos, iIndex, RbarUnit);
		if(RbarUnit.iDgnType==0) dRbarLen += RbarUnit.dLen;
	}
	return dRbarLen;
}

double CCalcPmcv::Get_RbarFixDgnArea()
{
	double dDgnArea=0.0;
	
	_UMD_RC_RBAR_UNIT RbarUnit;
	int iIndex=0;
	POSITION Pos = m_RbarData.arRbarUnit.GetStartPosition();
	while(Pos)
	{
		RbarUnit.Init();
		m_RbarData.arRbarUnit.GetNextAssoc(Pos, iIndex, RbarUnit);
		if(RbarUnit.iDgnType==1) dDgnArea += RbarUnit.dDgnArea;
	}
	return dDgnArea;
}

BOOL CCalcPmcv::Calc_PnMn(double dXn, double dInRota, double dAsThk, T_PM_SECT_RES& ResData, BOOL bRbarDetail)
{
	ResData.Init();
	double dPn = 0.0;	// Comp(+), Tens(-).
	double dMny = 0.0;
  double dMnz = 0.0;
  double desiMax = 0.0;
  double dTForce = 0.0; 
	double dCForce = 0.0;
  
	// Get Required Data from MyDB.
	double dYbar = m_dYbar;
	double dZbar = m_dZbar;
	double dFc	 = m_dFc;
	double dFyr	 = Get_ReBarStrength();
  
  double dAlpha = Calc_Alpha();

	BOOL bFixDgnRebarOnly = (m_bDesign && dAsThk<=0.0);
  

	double dDmax = Calc_Dmax(dInRota);
	if(Is_PureComp(dXn) || Is_PureTens(dXn))
	{
		double dCb;
		double dPncc=0.0, dPnsc=0.0, dPnst=0.0, dPrsc=0.0;	// c/s(Concrete/Steel), c/t(Compression/Tension).
		
		// by Concrete.
    double dAg = m_dArea;
		if(Is_PureTens(dXn))	
		{// Tens.
			dPncc = 0.0;		
			dCb   = 0.0;
			desiMax = cUMDRC_Upon;
		}
		else if(Is_PureComp(dXn))	
		{// Comp.
			dPncc = dAlpha*dFc*dAg;	
			dCb   = dDmax;
			desiMax = -cUMDRC_Upon;
		}
		// by Rebar.
		_UMD_RC_RBAR_UNIT	RbarUnit;
		
		int iIndexR=0;
		POSITION PosR = m_RbarData.arRbarUnit.GetStartPosition();
		while(PosR)
		{
			RbarUnit.Init();
			m_RbarData.arRbarUnit.GetNextAssoc(PosR,iIndexR,RbarUnit);
			double dAs=0.0;
			if(m_bDesign)	dAs = (RbarUnit.iDgnType==0) ? RbarUnit.dLen * dAsThk : RbarUnit.dDgnArea;
			else					dAs = RbarUnit.dArea;
			//대칭단면의 경우 모멘트는 0이므로 dMnxxx값에 대한 계산을 수행하지 않음
			if(Is_PureTens(dXn))	dPnst += (-1)*dFyr*dAs;	// Tension.
			else if(Is_PureComp(dXn))	// Compression.
			{
				dPnsc		+= dFyr*dAs;
				dPrsc   += dAlpha*dFc*dAs;
			}
		}
		
		if(m_bReduceRbar || Is_PureComp(dXn))
		{// 합성단면 계산시 단면공제 여부(FALSE=공제안함, TRUE=공제함), 순수 압축시에는 무조건 공제
			dPncc  -= dPrsc;
		}
		dPn = dPncc + dPnsc + dPnst;
    dTForce = fabs(dPnst);
    dCForce = fabs(dPncc + dPnsc);    	

		ResData.nType  = 1;//상세
		ResData.dNARotate = dInRota;
		ResData.dCb    = dCb;
		ResData.dXn    = dXn;
    ResData.dPn    = dPn;
    ResData.dPncc  = dPncc;
    ResData.dPnsc  = dPnsc;
    ResData.dPnst  = dPnst;
    ResData.dMny   = 0.0;
    ResData.dMnccy = 0.0;
    ResData.dMnscy = 0.0;
    ResData.dMnsty = 0.0;
    ResData.dMnz   = 0.0;
    ResData.dMnccz = 0.0;
    ResData.dMnscz = 0.0;
    ResData.dMnstz = 0.0; 
		ResData.dCForce= dCForce;
		ResData.dTForce= dTForce;
		ResData.desiMax= desiMax;  
		if(bRbarDetail)
		{ if(Get_Rebar_DatilData(dXn, dInRota*CMathFunc::m_trang, ResData.arRbarRes, dAsThk)) ResData.nType = 3; }
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

		double dBeta	= Calc_Beta();
		double dRota	= Calc_Rota(dInRota);
		double dDcomp	= dDmax*(1.0-dXn);
		double dDcent = dDmax*0.5;
		double dMaxDistTensRbar = Calc_MaxDistTensColmRbar(dXn, dRota, dDmax, bFixDgnRebarOnly);

		if(fabs(dRota) < cUMDRC_Zero)	// 0 Deg.
		{
			double dDz1 = dDmax*dXn;
      double dDy2 = 1/(cUMDRC_Zero);
			double dDz2 = (dDmax-dBeta*dDcomp);

			//////////////////////////
			// by Concrete.
			_UMD_RC_CONC_UNIT	ConcUnit;
			
			double dArea = 0.0;
			int iIndexC=0;
			POSITION PosC = m_ConcData.arConcUnit.GetStartPosition();
			while(PosC)
			{
				ConcUnit.Init();
				m_ConcData.arConcUnit.GetNextAssoc(PosC,iIndexC,ConcUnit);
				// Check if Concrete Compression Region is OK.
				double dCalcAc=0.0, dyInc=0.0, dzInc=0.0;
				if(Is_ConcComp(iIndexC,dRota,ConcUnit,dDy2,dDz2,dCalcAc,dyInc,dzInc))
				{
					double dCalcZc = ConcUnit.dyz[1] + dzInc;
					double dPncci = dAlpha*dFc * dCalcAc;
          dArea += dCalcAc;
					dPncc += dPncci;
					dTotNnccDistZ += fabs(dPncci * dCalcZc);
				}
			}
			//////////////////////////
			// by Rebar.
			_UMD_RC_RBAR_UNIT	RbarUnit;
						
			int iIndexR=0;
			POSITION PosR = m_RbarData.arRbarUnit.GetStartPosition();
			while(PosR)
			{
				RbarUnit.Init();
				m_RbarData.arRbarUnit.GetNextAssoc(PosR,iIndexR,RbarUnit);
				double dAs=0.0;
				if(m_bDesign)	dAs = (RbarUnit.iDgnType==0) ? RbarUnit.dLen * dAsThk : RbarUnit.dDgnArea;
				else					dAs = RbarUnit.dArea;
				double dDist = fabs(dDz1-RbarUnit.dyz[1]);
        
				// Check if compression or not.
				if(dDz1 < RbarUnit.dyz[1])	// Compression.
				{
          double desi = Calc_esi(dDcomp, dDist);
					double dPnsci = dAs * Calc_Fsi(TRUE, desi);
					dPnsc += dPnsci;
					dTotNnscDistZ += fabs(dPnsci * RbarUnit.dyz[1]);
				}
				else	// Tension.
				{
          double desi = Calc_esi(dDcomp, dDist);
          desiMax = max(desi, desiMax);
					double dPnsti = dAs * Calc_Fsi(FALSE, desi);
					dPnst += dPnsti;
					dTotNnstDistZ += fabs(dPnsti * RbarUnit.dyz[1]);
				}

				if(dDz2 < RbarUnit.dyz[1])	// Compression.
				{
					double dPrsci = dAlpha*dFc * dAs;
          dPrsc += dPrsci;
					dTotNrscDistZ += fabs(dPrsci * RbarUnit.dyz[1]);
				}
			}
			// Mncc, Mnsc, Mnst, Mrsc.
			if(fabs(dPncc) > 0.0)	dMnccy = fabs(dPncc * (dZbar-fabs(dTotNnccDistZ/dPncc)));
			if(fabs(dPnsc) > 0.0)	dMnscy = fabs(dPnsc * (dZbar-fabs(dTotNnscDistZ/dPnsc)));
			if(fabs(dPnst) > 0.0)	dMnsty = fabs(dPnst * (dZbar-fabs(dTotNnstDistZ/dPnst)));
			if(fabs(dPrsc) > 0.0)	dMrscy = fabs(dPrsc * (dZbar-fabs(dTotNrscDistZ/dPrsc)));
		}
		else if(fabs(dRota-CMathFunc::m_pi/2.0) < cUMDRC_Zero)	// 90 Deg.
		{
			double dDy1 = dDmax*dXn;
			double dDy2 = (dDmax-dBeta*dDcomp);
      double dDz2 = 1./(cUMDRC_Zero);
			//////////////////////////
			// by Concrete.
			_UMD_RC_CONC_UNIT	ConcUnit;
			
			double dArea = 0.0;
			int iIndexC=0;
			POSITION PosC = m_ConcData.arConcUnit.GetStartPosition();
			while(PosC)
			{
				ConcUnit.Init();
				m_ConcData.arConcUnit.GetNextAssoc(PosC,iIndexC,ConcUnit);
				// Check if Concrete Compression Region is OK.
				double dCalcAc=0.0, dyInc=0.0, dzInc=0.0;
				if(Is_ConcComp(iIndexC,dRota,ConcUnit,dDy2,dDz2,dCalcAc,dyInc,dzInc))
				{
					double dCalcYc = ConcUnit.dyz[0] + dyInc;
					double dPncci = dAlpha*dFc * dCalcAc;
          dArea += dCalcAc;
					dPncc += dPncci;
					dTotNnccDistY += fabs(dPncci * dCalcYc);
				}
			}

			//////////////////////////
			// by Rebar.
			_UMD_RC_RBAR_UNIT	RbarUnit;
			
			int iIndexR=0;
			POSITION PosR = m_RbarData.arRbarUnit.GetStartPosition();
			while(PosR)
			{
				RbarUnit.Init();
				m_RbarData.arRbarUnit.GetNextAssoc(PosR,iIndexR,RbarUnit);
				double dAs=0.0;
				if(m_bDesign)	dAs = (RbarUnit.iDgnType==0) ? RbarUnit.dLen * dAsThk : RbarUnit.dDgnArea;
				else					dAs = RbarUnit.dArea;
				double dDist = fabs(dDy1-RbarUnit.dyz[0]);
        
				// Check if compression or not.
				if(dDy1 < RbarUnit.dyz[0])	// Compression.
				{
          double desi = Calc_esi(dDcomp, dDist);
					double dPnsci = dAs * Calc_Fsi(TRUE, desi);
					dPnsc += dPnsci;
					dTotNnscDistY += fabs(dPnsci * RbarUnit.dyz[0]);

				}
				else	// Tension.
				{
          double desi = Calc_esi(dDcomp, dDist);
          desiMax = max(desi, desiMax);
					double dPnsti = dAs * Calc_Fsi(FALSE, desi);
					dPnst += dPnsti;
					dTotNnstDistY += fabs(dPnsti * RbarUnit.dyz[0]);
				}

				if(dDy2 < RbarUnit.dyz[0])	// Compression.
				{
					double dPrsci = dAlpha*dFc * dAs;
          dPrsc += dPrsci;
					dTotNrscDistY += fabs(dPrsci * RbarUnit.dyz[0]);
				}
			}
			// Mncc, Mnsc, Mnst, Mrsc.
			if(fabs(dPncc) > 0.0)	dMnccz = fabs(dPncc * (dYbar-fabs(dTotNnccDistY/dPncc)));
			if(fabs(dPnsc) > 0.0)	dMnscz = fabs(dPnsc * (dYbar-fabs(dTotNnscDistY/dPnsc)));
			if(fabs(dPnst) > 0.0)	dMnstz = fabs(dPnst * (dYbar-fabs(dTotNnstDistY/dPnst)));
			if(fabs(dPrsc) > 0.0)	dMrscz = fabs(dPrsc * (dYbar-fabs(dTotNrscDistY/dPrsc)));
		}
		else	// 0 ~ 90 Deg.
		{
      // SHIN ('06.04.13) 단면형상에 따라 콘크리트와 철근의 중립축이 일치 하지 않아 수정함       
      double dDistMax=0.0, dDistMin=0.0;

      Calc_DistanceFromSectionToAngleLine(0, 0, dRota, dDistMax, dDistMin);
      dDmax = dDistMax-dDistMin;
      dDcomp = dDmax*(1.0-dXn);
      double dDy1 = (dDistMax-dDcomp) / sin(dRota);
			double dDz1 = (dDistMax-dDcomp) / cos(dRota);
			double dDy2 = (dDistMax-dBeta*dDcomp) / sin(dRota);
			double dDz2 = (dDistMax-dBeta*dDcomp) / cos(dRota);
      dMaxDistTensRbar = Calc_MaxDistTensColmRbar2(dDy1, dDz1, dRota, bFixDgnRebarOnly);
      
			//////////////////////////
			// by Concrete.
			_UMD_RC_CONC_UNIT	ConcUnit;
			
			double dArea = 0.0;
			int iIndexC=0;
			POSITION PosC = m_ConcData.arConcUnit.GetStartPosition();
			while(PosC)
			{
				ConcUnit.Init();
				m_ConcData.arConcUnit.GetNextAssoc(PosC,iIndexC,ConcUnit);
				// Check if Concrete Compression Region is OK.
				double dCalcAc=0.0, dyInc=0.0, dzInc=0.0;
				if(Is_ConcComp(iIndexC,dRota,ConcUnit,dDy2,dDz2,dCalcAc,dyInc,dzInc))
				{
					double dCalcYc = ConcUnit.dyz[0] + dyInc;
					double dCalcZc = ConcUnit.dyz[1] + dzInc;
					double dPncci = dAlpha*dFc * dCalcAc;
          dArea += dCalcAc;
					dPncc += dPncci;
					dTotNnccDistY += fabs(dPncci * dCalcYc);
					dTotNnccDistZ += fabs(dPncci * dCalcZc);
				}
			}

			//////////////////////////
			// by Rebar.
			_UMD_RC_RBAR_UNIT	RbarUnit;
			
			int iIndexR=0;
			POSITION PosR = m_RbarData.arRbarUnit.GetStartPosition();
			while(PosR)
			{
				RbarUnit.Init();
				m_RbarData.arRbarUnit.GetNextAssoc(PosR,iIndexR,RbarUnit);
				double dAs=0.0;
				if(m_bDesign)	dAs = (RbarUnit.iDgnType==0) ? RbarUnit.dLen * dAsThk : RbarUnit.dDgnArea;
				else					dAs = RbarUnit.dArea;

        if(dDy1 == 0.0 && dDz1 == 0.0)
        {
			    dDy1 = cUMDRC_Zero / sin(dRota);
			    dDz1 = cUMDRC_Zero / cos(dRota);
        }

				double dLineI[3]={dDy1,  0.0, 0.0};
				double dLineJ[3]={ 0.0, dDz1, 0.0};
				double dPoint[3]={RbarUnit.dyz[0], RbarUnit.dyz[1], 0.0};
				double dDist = CMathFunc::mathDistanceFromIntersectPointToLine(dLineI,dLineJ,dPoint);
				
				// Check if compression or not.
        double dCoordZ = (-1)*tan(dRota)*RbarUnit.dyz[0] + dDz1;	// by Neutral axis.
				if(dCoordZ < RbarUnit.dyz[1])	// Compression.
				{
          double desi = Calc_esi(dDcomp, dDist);
					double dPnsci = dAs * Calc_Fsi(TRUE, desi);
					dPnsc += dPnsci;
					dTotNnscDistY += fabs(dPnsci * RbarUnit.dyz[0]);
					dTotNnscDistZ += fabs(dPnsci * RbarUnit.dyz[1]);
				}
				else	// Tension.
				{
          double desi = Calc_esi(dDcomp, dDist);
          desiMax = max(desi, desiMax);
					double dPnsti = dAs * Calc_Fsi(FALSE, desi);
					dPnst += dPnsti;
					dTotNnstDistY += fabs(dPnsti * RbarUnit.dyz[0]);
					dTotNnstDistZ += fabs(dPnsti * RbarUnit.dyz[1]);
				}
				// Reduce Compression Rebar Area.        
        double dCoordZc = (-1)*tan(dRota)*RbarUnit.dyz[0] + dDz2;	// by Neutral axis.
				if(dCoordZc < RbarUnit.dyz[1])
				{
					double dPrsci = dAlpha*dFc * dAs;
					dPrsc += dPrsci;
					dTotNrscDistY += fabs(dPrsci * RbarUnit.dyz[0]);
					dTotNrscDistZ += fabs(dPrsci * RbarUnit.dyz[1]);
				}
			}

			if(fabs(dPncc) > 0.0)
			{
        double dDistY = dYbar - fabs(dTotNnccDistY/dPncc);
        double dDistZ = dZbar - fabs(dTotNnccDistZ/dPncc);
				dMnccy = fabs(dPncc * dDistZ);
        dMnccz = fabs(dPncc * dDistY);
			}
			// Mnsc (Steel Compression).
			if(fabs(dPnsc) > 0.0)
			{
        double dDistY = dYbar - fabs(dTotNnscDistY/dPnsc);
        double dDistZ = dZbar - fabs(dTotNnscDistZ/dPnsc);
				dMnscy = fabs(dPnsc * dDistZ);
        dMnscz = fabs(dPnsc * dDistY);
			}
			// Mnst (Steel Tension).
			if(fabs(dPnst) > 0.0)
			{
        double dDistY = dYbar - fabs(dTotNnstDistY/dPnst);
        double dDistZ = dZbar - fabs(dTotNnstDistZ/dPnst);
				dMnsty = fabs(dPnst * dDistZ);
        dMnstz = fabs(dPnst * dDistY);
			}
      // Mrsc (Steel Compression - Reduced dupulicated area).
			if(fabs(dPrsc) > 0.0)
			{
        double dDistY = dYbar - fabs(dTotNrscDistY/dPrsc);
        double dDistZ = dZbar - fabs(dTotNrscDistZ/dPrsc);
				dMrscy = fabs(dPrsc * dDistZ);
        dMrscz = fabs(dPrsc * dDistY);
			}
		}
		// Calculate Pn, Mn.
		if(m_bReduceRbar)
		{// 합성단면 계산시 단면공제 여부(FALSE=공제안함, TRUE=공제함)
			dPncc  -= dPrsc;
			dMnccy -= dMrscy;
			dMnccz -= dMrscz;
		}
		dPn = dPncc + dPnsc + dPnst;
		dMny = dMnccy + dMnscy + dMnsty;
    dMnz = dMnccz + dMnscz + dMnstz;
    dTForce = fabs(dPnst);
    dCForce = fabs(dPncc + dPnsc);
    
    //
		ResData.nType  = 1;//상세
		ResData.dNARotate = dInRota;
    ResData.dCb    = (1.0-dXn)*dDmax;
		ResData.dXn    = dXn;
    ResData.dPn    = dPn;
    ResData.dPncc  = dPncc;
    ResData.dPnsc  = dPnsc;
    ResData.dPnst  = dPnst;
    ResData.dMny   = dMny;
    ResData.dMnccy = dMnccy;
    ResData.dMnscy = dMnscy;
    ResData.dMnsty = dMnsty;
    ResData.dMnz   = dMnz;
    ResData.dMnccz = dMnccz;
    ResData.dMnscz = dMnscz;
    ResData.dMnstz = dMnstz; 
		ResData.dCForce= dCForce;
		ResData.dTForce= dTForce;
		ResData.desiMax= desiMax;    
		if(bRbarDetail)
		{ if(Get_Rebar_DatilData(dXn, dInRota*CMathFunc::m_trang, ResData.arRbarRes, dAsThk)) ResData.nType = 3; }		
	}
	return TRUE;
}

double CCalcPmcv::Calc_Deff_Gen(double dDmax, double dZmin, _UMD_RC_RBAR_LIST& arRbarUnit, BOOL bFixDgnRebarOnly)
{
  double dComp, dTens;
  //최상단에 중립축이 위치할때의 인장측 철근까지의 거리가 유효 깊이임
  if(!Cal_MaxDistRbar_Gen(1.0, dDmax, dZmin, arRbarUnit, dTens, dComp, bFixDgnRebarOnly)) return 0.0;
  return dTens; 
}

double CCalcPmcv::Get_Xn(int iCurDiv, int iDivNum, double dXb)
{
	double dXn=0.0;
	if(iDivNum%2!=0)	return dXn;	// Only even.
	if(iDivNum < 3)		return dXn;	// Minimum 2.
	
	int iDivHalf = iDivNum/2;
	double dX		 = (iCurDiv <= iDivHalf ? iCurDiv : iCurDiv-iDivHalf)/(double)iDivHalf;
	double dTag1 = (iCurDiv <= iDivHalf ? 0.0 : 1.0);	// Control Region(0=Comp, 1=Tens).
	double dTag2 = (iCurDiv <= iDivHalf ? 1.0 : 0.0);	// Control Region(1=Comp, 0=Tens).
	int iShp=m_iSectShp;
	
	double dUp=0.0;
	if     (iShp==DGN_SECT_SHAPE_INDEX_REG_B)			dUp = 0.6;
	else if(iShp==DGN_SECT_SHAPE_INDEX_REG_P)     dUp = 1.2;
	else if(iShp==DGN_SECT_SHAPE_INDEX_REG_SR)		dUp = 1.5;
	else if(iShp==DGN_SECT_SHAPE_INDEX_REG_SB)		dUp = 1.2;
	else if(iShp==DGN_SECT_SHAPE_INDEX_REG_OCT)	  dUp = 0.6;
	else if(iShp==DGN_SECT_SHAPE_INDEX_REG_SOCT)	dUp = 1.2;
	else if(iShp==DGN_SECT_SHAPE_INDEX_REG_TRK)	  dUp = 1.2;
	else if(iShp==DGN_SECT_SHAPE_INDEX_REG_STRK)	dUp = 1.5;
	else if(iShp==DGN_SECT_SHAPE_INDEX_REG_HTRK)	dUp = 1.2;
	else if(iShp==DGN_SECT_SHAPE_INDEX_REG_BSTF)	dUp = 0.6;
	else if(iShp==DGN_SECT_SHAPE_INDEX_REG_PSTF)	dUp = 1.2;
	//else if(iShp==DGN_SECT_SHAPE_INDEX_REG_PSC)	  dUp = 0.6;
  else if(iShp==DGN_SECT_SHAPE_INDEX_REG_GEN)   dUp = 1.2;// SHIN ('06.02.15) SB과 같이 만들었을 뿐 근거 없음
	else	ASSERT(0);
	// y = ax^dUp + b.
	dXn = (dTag1-dXb)*pow(fabs(dTag2-dX),dUp) + dXb;	// Distance Ratio from Tensile edge.
	return dXn;
}


BOOL CCalcPmcv::Calc_PnMn_Gen_SliceMethod(double dXn, double dInRota, double dAsThk, T_PM_SECT_RES& ResData, 
		                           _UMD_RC_CON_PART_LIST& arConcUnit, _UMD_RC_RBAR_LIST& arRbarUnit, BOOL bRbarDetail)
{
	//초기화
	ResData.Init();
	// c/s/r(Concrete/Rebar/철근면적 공재), c/t(Compression/Tension), _y/_z(with Y-Axis, with Z-Axis).		
	double dPn   = 0.0;// 공칭 축하중
  double dMny  = 0.0;// Y축에 대한 공칭 모멘트(전체좌표계)
  double dMnz  = 0.0;// Z축에 대한 공칭 모멘트(전체좌표계)
  double den   = 0.0;// 압축부의 최대변형률(인장시에는 콘크리트의 극한변형률) 
  double desi  = 0.0;// 인장부철근의 최대변형률(압축시 '-')
  double dEccen= 0.0;// 편심거리  
  double dCn   = 0.0;// 압축부에서부터 중립축까지의 거리
  double dCc   = 0.0;// 압축부 콘크리트가 부담하는 압축력(+)
  double dCs   = 0.0;// 압축부 철근이     부담하는 압축력(+)
	double dRs   = 0.0;// 압축부 철근면적에 공재되는 압축력(+)
  double dTs   = 0.0;// 인장부 철근이     부담하는 인장력(-)
  double dMccy = 0.0;// 압축부 콘크리트가 부담하는 Y축에 대한 모멘트(+)
  double dMcsy = 0.0;// 압축부 철근이     부담하는 Y축에 대한 모멘트(+)
  double dMtsy = 0.0;// 인장부 철근이     부담하는 Y축에 대한 모멘트(+)
	double dMrsy = 0.0;// 압축부 철근면적에 공재되는 Y축에 대한 모멘트(+)
  double dMccz = 0.0;// 압축부 콘크리트가 부담하는 Z축에 대한 모멘트(+)
  double dMcsz = 0.0;// 압축부 철근이     부담하는 Z축에 대한 모멘트(+)
  double dMtsz = 0.0;// 인장부 철근이     부담하는 Z축에 대한 모멘트(+)
	double dMrsz = 0.0;// 압축부 철근면적에 공재되는 Z축에 대한 모멘트(+)

  double dFc		= m_dFc;	// Concrete Design Strength.
	double dFyr		= Get_ReBarStrength();	// Rbar Strength.	
  double dSigcc = Get_Sigcc();
	double dDmax  = arConcUnit.dZmax- arConcUnit.dZmin;
	double dZmin  = arConcUnit.dZmin;
  double dRbarZmin = dDmax;
  double dEffectiveD = dDmax;

	double dRotCenter_z = m_dZbar;
	double dRotCenter_y = m_dYbar;

	_UMD_RC_CON_PART_UNIT ConUnit;
	_UMD_RC_RBAR_UNIT RbarUnit;

	BOOL bFixDgnRebarOnly = (m_bDesign && dAsThk<=0.0);
	
	int i=0, j=0;
  if(Is_PureComp(dXn))	// dXn = 0.
  {//순수압축
    den  = Get_ecu();  
		desi = -den;
		//////////////////////////
		// by Concrete.
		for(int i=0; i<arConcUnit.List.GetSize(); i++)
		{
			ConUnit.Init();
			ConUnit = arConcUnit.List.GetAt(i);
			double dArea = ConUnit.dArea;
      double dCeny = ConUnit.dyz[0];
			double dCenz = ConUnit.dyz[1];
      // 회전중심점에서의 거리
			double dCalz = dCenz - dRotCenter_z; 
      double dCaly = dCeny - dRotCenter_y;      
      double dFcc  = dSigcc*dArea;
			dCc += dFcc;	// Compression.
      // 순수 압축이라도 회전중심점에서 모멘트 발생가능
      dMccy += dFcc*dCalz;      
      dMccz += dFcc*dCaly;
		}
		//////////////////////////
		// by Rebar.
    for(int i=0; i<arRbarUnit.List.GetSize(); i++)
    {
      RbarUnit.Init();
      RbarUnit = arRbarUnit.List.GetAt(i);
      double dAs=0.0;
			if(m_bDesign)	dAs = (RbarUnit.iDgnType==0) ? RbarUnit.dLen * dAsThk : RbarUnit.dDgnArea;
			else					dAs = RbarUnit.dArea;
      double dCeny = RbarUnit.dyz[0];
			double dCenz = RbarUnit.dyz[1];

      if(dRbarZmin > dCenz) dRbarZmin = dCenz;
      // 회전중심점에서의 거리
			double dCalz = dCenz - dRotCenter_z; 
      double dCaly = dCeny - dRotCenter_y;
			double dFcs = dFyr*dAs;
			double dFrs = dSigcc*dAs;
      dCs += dFcs;	// Compression.
			dRs += dFrs;
      // 순수 압축이라도 회전중심점에서 모멘트 발생가능			
      dMcsy += dFcs*dCalz;
      dMcsz += dFcs*dCaly;
			dMrsy += dFrs*dCalz;
      dMrsz += dFrs*dCaly;
    }	
  }
	else if(Is_PureTens(dXn))	// dXn = 1.
  {// 순수인장
    den = Get_esu();  desi = den;

		//////////////////////////
		// by Rbar.
    for(int i=0; i<arRbarUnit.List.GetSize(); i++)
    {
      RbarUnit.Init();
      RbarUnit = arRbarUnit.List.GetAt(i);

      double dAs=0.0;
			if(m_bDesign)	dAs = (RbarUnit.iDgnType==0) ? RbarUnit.dLen * dAsThk : RbarUnit.dDgnArea;
			else					dAs = RbarUnit.dArea;
      double dCeny = RbarUnit.dyz[0];
			double dCenz = RbarUnit.dyz[1];

      if(dRbarZmin > dCenz) dRbarZmin = dCenz;
      // 회전중심점에서의 거리
			double dCalz = dCenz - dRotCenter_z; 
      double dCaly = dCeny - dRotCenter_y;
      double dFr   = (-1)*dFyr*dAs;
			dTs += dFr;	// Tension.
      // 순수 인장이라도 회전중심점에서 모멘트 발생가능
      dMtsy += dFr*dCalz;
      dMtsz += dFr*dCaly;
    }		
  }
  else	// 0 < dXn < 1.
  {		
		double dDcomp	= dDmax*(1.0-dXn);// 압축부 길이

		double dMaxDistTensRbar=0.0, dMaxDistCompRbar=0.0;
		// Calculate en (최외측에서의 변형률, 각 Part의 중심아님).
		double dec=0.0, det=0.0;
		double dDc=dDcomp;
		
		BOOL bMaxDistRbar = Cal_MaxDistRbar_Gen(dXn, dDmax, dZmin, arRbarUnit, dMaxDistTensRbar,dMaxDistCompRbar, bFixDgnRebarOnly);
		
    Calc_SigcU(dDcomp,dDc,dMaxDistCompRbar,dec);
    Calc_SigrU(dDcomp,dMaxDistTensRbar,dMaxDistCompRbar,det); 

		den  = dec;//압축부의 변위
    desi = det;//인장 철근의 변위
				
    double dDy1 = cUMDRC_Upon;
		double dDz1 = dDmax*dXn;
    double dDistEcc = dDcomp*(Get_ecc()/Get_ecu());
    double dHalfThick = dDmax/(2.0*double(m_iSliceNum));
    //////////////////////////
		// by Concrete.
		for(int i=0; i<arConcUnit.List.GetSize(); i++)
		{
			ConUnit.Init();
			ConUnit = arConcUnit.List.GetAt(i);    
			double dArea = ConUnit.dArea;
      double dCeny = ConUnit.dyz[0];
			double dCenz = ConUnit.dyz[1];
      double dDist = dCenz - dZmin - dDz1;	// Comp(+), Tens(-). 중립축에서 떨어진 거리
			double dEpsc = 0.0;
			double dSigc;
      if((m_iTypeStressStrain == 3 || m_iTypeStressStrain == 4) && fabs(dDist-dDistEcc) < dHalfThick)
      { // SHIN ('06.03.10)
				// 기본적으로 응력블럭일 경우 Calc_PnMn_Gen_SliceMethod가 아닌 Calc_PnMn_Gen_RectStress를 사용한다.
        //   ConUnit는 중심점과 면적만을 가지고 있으므로 실제 단면의 요소 형상을 알수는 없다... 하지만 슬라이스 수가 작을경우 
        //   응력블럭으로 하였을 경우 오차범위가 커지므로 2*dHalfThick크기의 직사각형으로 가정하고 보정해주기로 한다. 
        //   따라서 실제와는 다를수 있으나 슬라이스의 수가 작으므으로 생길수 있는 오차 보다는 작다고 볼수 있다.
        double dThickTop    = dDist+dHalfThick - dDistEcc;
        double dThickBottom = 2.0*dHalfThick - dThickTop;
        double dDistTop     = dDistEcc + dThickTop/2.0;
        double dDistBottom  = dDistEcc - dThickBottom/2.0;
        double dSigcTop     = Calc_SigcU(dDcomp,dDistTop,dMaxDistCompRbar,dEpsc);
        double dSigcBottom  = Calc_SigcU(dDcomp,dDistBottom,dMaxDistCompRbar,dEpsc);
        dSigc = (dThickTop*dSigcTop + dThickBottom*dSigcBottom) / (2.0*dHalfThick);
        dCenz = dCenz-dHalfThick + (dThickTop*dSigcTop*(2.0*dHalfThick-dThickTop/2.0) + dThickBottom*dSigcBottom*dThickBottom/2.0) / (dSigc*2.0*dHalfThick);
      }
      else
        dSigc = Calc_SigcU(dDcomp,dDist,dMaxDistCompRbar,dEpsc);
			double dFcc  = dSigc*dArea;
			// 회전중심점에서의 거리
			double dCalz = dCenz - dRotCenter_z; 
      double dCaly = dCeny - dRotCenter_y;
			// Comp(dDist > 0), Tens(dDist < 0).
			if(dDist > 0.0) {dCc += dFcc; dMccy += dFcc*dCalz; dMccz += dFcc*dCaly;}			
		}
		//////////////////////////
		// by Rebar.
		for(int i=0; i<arRbarUnit.List.GetSize(); i++)
		{
			RbarUnit.Init();
			RbarUnit = arRbarUnit.List.GetAt(i);
			double dArea=0.0;
			if(m_bDesign)	dArea = (RbarUnit.iDgnType==0) ? RbarUnit.dLen * dAsThk : RbarUnit.dDgnArea;
			else					dArea = RbarUnit.dArea;
      double dCeny = RbarUnit.dyz[0];
			double dCenz = RbarUnit.dyz[1];
			double dDist = dCenz - dZmin - dDz1;	// Comp(+), Tens(-).
			double dEpsr = 0.0;
			double dSigr = Calc_SigrU(dDcomp,dDist,dMaxDistCompRbar,dEpsr);
			double dFr	 = dSigr*dArea;

      if(dRbarZmin > dCenz) dRbarZmin = dCenz;
			// Change by ZINU.('05.05.07). m_dCzm -> m_dCzm or m_dCzp by iAxisDir.
			double dCalz = dCenz - dRotCenter_z;// 회전중심점에서의 거리
      double dCaly = dCeny - dRotCenter_y;
			// Comp(dDist > 0), Tens(dDist < 0).
			if(dDist > 0.0)	
			{// Comp.
				dCs += dFr; dMcsy += dFr*dCalz; dMcsz += dFr*dCaly;
				double dEpsc = 0.0;
				double dFrs  = dArea*Calc_SigcU(dDcomp,dDist,dMaxDistCompRbar,dEpsc);
				dRs += dFrs; dMrsy += dFrs*dCalz; dMrsz += dFrs*dCaly;
			}	
			else						
			{// Tens
				dTs += dFr; dMtsy += dFr*dCalz; dMtsz += dFr*dCaly;
			}	
		}		
   
	}

	// Calculate Pn, Mn.
	if(m_bReduceRbar || Is_PureComp(dXn))
	{// 합성단면 계산시 단면공제 여부(FALSE=공제안함, TRUE=공제함), 순수 압축시에는 무조건 공제
		dCc   -= dRs;
		dMccy -= dMrsy;
		dMccz -= dMrsz;
	}
	dPn = dCc + dCs + dTs;	
  dMny = dMccy + dMcsy + dMtsy;
  dMnz = dMccz + dMcsz + dMtsz;    
	double dTForce = fabs(dTs);
  double dCForce = fabs(dCc + dCs);

  //모멘트의 경우 회전중심을 중심으로 중립축의 직각좌표축에 대한 값이므로 전체 좌표계로 수정
  double dMnx=0.0;
	double dRotation = dInRota*CMathFunc::m_trang;

  CMathFunc::mathRotateX(dRotation, dMnx, dMnz,  dMny);  dMnx=0.0;
  CMathFunc::mathRotateX(dRotation, dMnx, dMccz, dMccy); dMnx=0.0;
  CMathFunc::mathRotateX(dRotation, dMnx, dMcsz, dMcsy); dMnx=0.0;
  CMathFunc::mathRotateX(dRotation, dMnx, dMtsz, dMtsy);
  
	ResData.nType  = 1;
	ResData.dCb    = (1.0-dXn)*dDmax;
	ResData.dXn    = dXn;
  ResData.dPn    = dPn;
  ResData.dPncc  = dCc;
  ResData.dPnsc  = dCs;
  ResData.dPnst  = dTs;
  ResData.dMny   = dMny;
  ResData.dMnccy = dMccy;
  ResData.dMnscy = dMcsy;
  ResData.dMnsty = dMtsy;
  ResData.dMnz   = dMnz;
  ResData.dMnccz = dMccz;
  ResData.dMnscz = dMcsz;
  ResData.dMnstz = dMtsz; 
	ResData.dCForce= dCForce;
	ResData.dTForce= dTForce;
	ResData.desiMax= desi;    
	if(bRbarDetail)
	{ if(Get_Rebar_DatilData(dXn, dRotation, ResData.arRbarRes, dAsThk)) ResData.nType = 3; }

	return TRUE;		
}

BOOL CCalcPmcv::Calc_PnMn_Gen_RectStress(double dXn, double dInRota, double dAsThk, T_PM_SECT_RES& ResData, 
		                          double dYmax, double dYmin, double dZmax, double dZmin, CPolyMaker* pPolyMaker, _UMD_RC_RBAR_LIST& arRbarUnit, BOOL bRbarDetail)
{
	//초기화
	ResData.Init();
	// c/s/r(Concrete/Rebar/철근면적 공재), c/t(Compression/Tension), _y/_z(with Y-Axis, with Z-Axis).		
	double dPn   = 0.0;// 공칭 축하중
  double dMny  = 0.0;// Y축에 대한 공칭 모멘트(전체좌표계)
  double dMnz  = 0.0;// Z축에 대한 공칭 모멘트(전체좌표계)
  double den   = 0.0;// 압축부의 최대변형률(인장시에는 콘크리트의 극한변형률) 
  double desi  = 0.0;// 인장부철근의 최대변형률(압축시 '-')
  double dEccen= 0.0;// 편심거리  
  double dCn   = 0.0;// 압축부에서부터 중립축까지의 거리
  double dCc   = 0.0;// 압축부 콘크리트가 부담하는 압축력(+)
  double dCs   = 0.0;// 압축부 철근이     부담하는 압축력(+)
	double dRs   = 0.0;// 압축부 철근면적에 공재되는 압축력(+)
  double dTs   = 0.0;// 인장부 철근이     부담하는 인장력(-)
  double dMccy = 0.0;// 압축부 콘크리트가 부담하는 Y축에 대한 모멘트(+)
  double dMcsy = 0.0;// 압축부 철근이     부담하는 Y축에 대한 모멘트(+)
  double dMtsy = 0.0;// 인장부 철근이     부담하는 Y축에 대한 모멘트(+)
	double dMrsy = 0.0;// 압축부 철근면적에 공재되는 Y축에 대한 모멘트(+)
  double dMccz = 0.0;// 압축부 콘크리트가 부담하는 Z축에 대한 모멘트(+)
  double dMcsz = 0.0;// 압축부 철근이     부담하는 Z축에 대한 모멘트(+)
  double dMtsz = 0.0;// 인장부 철근이     부담하는 Z축에 대한 모멘트(+)
	double dMrsz = 0.0;// 압축부 철근면적에 공재되는 Z축에 대한 모멘트(+)

  double dFc		= m_dFc;	// Concrete Design Strength.
	double dFyr		= Get_ReBarStrength();	// Rbar Strength.
  double dSigcc = Get_Sigcc();
	double dDmax  = dZmax- dZmin;
  double dRbarZmin = dDmax;
  double dEffectiveD = dDmax;
	
	double dRotCenter_z = m_dZbar;
	double dRotCenter_y = m_dYbar;

	BOOL bFixDgnRebarOnly = (m_bDesign && dAsThk<=0.0);

	_UMD_RC_CON_PART_UNIT ConUnit;
	_UMD_RC_RBAR_UNIT RbarUnit;
	
	CMdgnSectTool SectTool;
	int i=0, j=0;
  if(Is_PureComp(dXn))	// dXn = 0.
  {//순수압축
    den  = Get_ecu();  
		desi = -den;
		double dZ1 = dZmax+cUMDRC_Zero;
		double dZ2 = dZmin-cUMDRC_Zero;
		//////////////////////////
		// by Concrete.	
		if(pPolyMaker != NULL)
		{
			double dArea, dCeny, dCenz;
			dArea = dCeny = dCenz = 0.0;
			SectTool.Get_PolyDataCentroid(dZ1, dZ2, dYmax, dYmin, pPolyMaker, dCeny, dCenz, dArea);
			// 회전중심점에서의 거리
			double dCalz = dCenz - dRotCenter_z; 
			double dCaly = dCeny - dRotCenter_y;      
			double dFcc  = dSigcc*dArea;
			dCc = dFcc;	// Compression.
			// 순수 압축이라도 회전중심점에서 모멘트 발생가능
			dMccy = dFcc*dCalz;      
			dMccz = dFcc*dCaly;
		}
				
		//////////////////////////
		// by Rebar.
    for(int i=0; i<arRbarUnit.List.GetSize(); i++)
    {
      RbarUnit.Init();
      RbarUnit = arRbarUnit.List.GetAt(i);
      double dAs=0.0;
			if(m_bDesign)	dAs = (RbarUnit.iDgnType==0) ? RbarUnit.dLen * dAsThk : RbarUnit.dDgnArea;
			else					dAs = RbarUnit.dArea;
      double dCeny = RbarUnit.dyz[0];
			double dCenz = RbarUnit.dyz[1];

      if(dRbarZmin > dCenz) dRbarZmin = dCenz;
      // 회전중심점에서의 거리
			double dCalz = dCenz - dRotCenter_z; 
      double dCaly = dCeny - dRotCenter_y;
			double dFcs = dFyr*dAs;
			double dFrs = dSigcc*dAs;
      dCs += dFcs;	// Compression.
			dRs += dFrs;
      // 순수 압축이라도 회전중심점에서 모멘트 발생가능			
      dMcsy += dFcs*dCalz;
      dMcsz += dFcs*dCaly;
			dMrsy += dFrs*dCalz;
      dMrsz += dFrs*dCaly;
    }			
  }
	else if(Is_PureTens(dXn))	// dXn = 1.
  {// 순수인장
    den = Get_esu();  desi = den;

		//////////////////////////
		// by Rbar.
    for(int i=0; i<arRbarUnit.List.GetSize(); i++)
    {
      RbarUnit.Init();
      RbarUnit = arRbarUnit.List.GetAt(i);

      double dAs=0.0;
			if(m_bDesign)	dAs = (RbarUnit.iDgnType==0) ? RbarUnit.dLen * dAsThk : RbarUnit.dDgnArea;
			else					dAs = RbarUnit.dArea;
      double dCeny = RbarUnit.dyz[0];
			double dCenz = RbarUnit.dyz[1];

      if(dRbarZmin > dCenz) dRbarZmin = dCenz;
      // 회전중심점에서의 거리
			double dCalz = dCenz - dRotCenter_z; 
      double dCaly = dCeny - dRotCenter_y;
      double dFr   = (-1)*dFyr*dAs;
			dTs += dFr;	// Tension.
      // 순수 인장이라도 회전중심점에서 모멘트 발생가능
      dMtsy += dFr*dCalz;
      dMtsz += dFr*dCaly;
    }		
  }
  else	// 0 < dXn < 1.
  {		
		double dDcomp	= dDmax*(1.0-dXn);// 압축부 길이
		double dBeta = Calc_Beta();
		double dZ1 = dZmax+cUMDRC_Zero;
		double dZ2 = dZmax-dBeta*dDcomp;

		double dMaxDistTensRbar=0.0, dMaxDistCompRbar=0.0;
		// Calculate en (최외측에서의 변형률, 각 Part의 중심아님).
		double dec=0.0, det=0.0;
		double dDc=dDcomp;
		
		BOOL bMaxDistRbar = Cal_MaxDistRbar_Gen(dXn, dDmax, dZmin, arRbarUnit, dMaxDistTensRbar,dMaxDistCompRbar, bFixDgnRebarOnly);
		
    Calc_SigcU(dDcomp,dDc,dMaxDistCompRbar,dec);
    Calc_SigrU(dDcomp,dMaxDistTensRbar,dMaxDistCompRbar,det); 

		den  = dec;//압축부의 변위
    desi = det;//인장 철근의 변위
				
    double dDy1 = cUMDRC_Upon;
		double dDz1 = dDmax*dXn;
    //////////////////////////
		// by Concrete.		
		if(pPolyMaker != NULL)
		{
			double dArea, dCeny, dCenz;
			dArea = dCeny = dCenz = 0.0;
			SectTool.Get_PolyDataCentroid(dZ1, dZ2, dYmax, dYmin, pPolyMaker, dCeny, dCenz, dArea);
			// 회전중심점에서의 거리
			double dCalz = dCenz - dRotCenter_z; 
			double dCaly = dCeny - dRotCenter_y;      
			double dFcc  = dSigcc*dArea;
			dCc = dFcc;	// Compression.
			// 순수 압축이라도 회전중심점에서 모멘트 발생가능
			dMccy = dFcc*dCalz;      
			dMccz = dFcc*dCaly;
		}
				
		//////////////////////////
		// by Rebar.
		for(int i=0; i<arRbarUnit.List.GetSize(); i++)
		{
			RbarUnit.Init();
			RbarUnit = arRbarUnit.List.GetAt(i);
			double dArea=0.0;
			if(m_bDesign)	dArea = (RbarUnit.iDgnType==0) ? RbarUnit.dLen * dAsThk : RbarUnit.dDgnArea;
			else					dArea = RbarUnit.dArea;
      double dCeny = RbarUnit.dyz[0];
			double dCenz = RbarUnit.dyz[1];
			double dDist = dCenz - dZmin - dDz1;	// Comp(+), Tens(-).
			double dEpsr = 0.0;
			double dSigr = Calc_SigrU(dDcomp,dDist,dMaxDistCompRbar,dEpsr);
			double dFr	 = dSigr*dArea;

      if(dRbarZmin > dCenz) dRbarZmin = dCenz;
			// Change by ZINU.('05.05.07). m_dCzm -> m_dCzm or m_dCzp by iAxisDir.
			double dCalz = dCenz - dRotCenter_z;// 회전중심점에서의 거리
      double dCaly = dCeny - dRotCenter_y;
			// Comp(dDist > 0), Tens(dDist < 0).
			if(dDist > 0.0)	
			{// Comp.
				dCs += dFr; dMcsy += dFr*dCalz; dMcsz += dFr*dCaly;
				double dEpsc = 0.0;
				double dFrs  = dArea*Calc_SigcU(dDcomp,dDist,dMaxDistCompRbar,dEpsc);
				dRs += dFrs; dMrsy += dFrs*dCalz; dMrsz += dFrs*dCaly;
			}	
			else						
			{// Tens
				dTs += dFr; dMtsy += dFr*dCalz; dMtsz += dFr*dCaly;
			}	
		}		
   
	}

	// Calculate Pn, Mn.
	if(m_bReduceRbar || Is_PureComp(dXn))
	{// 합성단면 계산시 단면공제 여부(FALSE=공제안함, TRUE=공제함)
		dCc   -= dRs;
		dMccy -= dMrsy;
		dMccz -= dMrsz;
	}
	dPn = dCc + dCs + dTs;	
  dMny = dMccy + dMcsy + dMtsy;
  dMnz = dMccz + dMcsz + dMtsz;    
	double dTForce = fabs(dTs);
  double dCForce = fabs(dCc + dCs);

  //모멘트의 경우 회전중심을 중심으로 중립축의 직각좌표축에 대한 값이므로 전체 좌표계로 수정
  double dMnx=0.0;
	double dRotation = dInRota*CMathFunc::m_trang;
  CMathFunc::mathRotateX(dRotation, dMnx, dMnz,  dMny);  dMnx=0.0;
  CMathFunc::mathRotateX(dRotation, dMnx, dMccz, dMccy); dMnx=0.0;
  CMathFunc::mathRotateX(dRotation, dMnx, dMcsz, dMcsy); dMnx=0.0;
  CMathFunc::mathRotateX(dRotation, dMnx, dMtsz, dMtsy);
  
	ResData.nType  = 1;
	ResData.dCb    = (1.0-dXn)*dDmax;
	ResData.dXn    = dXn;
  ResData.dPn    = dPn;
  ResData.dPncc  = dCc;
  ResData.dPnsc  = dCs;
  ResData.dPnst  = dTs;
  ResData.dMny   = dMny;
  ResData.dMnccy = dMccy;
  ResData.dMnscy = dMcsy;
  ResData.dMnsty = dMtsy;
  ResData.dMnz   = dMnz;
  ResData.dMnccz = dMccz;
  ResData.dMnscz = dMcsz;
  ResData.dMnstz = dMtsz; 
	ResData.dCForce= dCForce;
	ResData.dTForce= dTForce;
	ResData.desiMax= desi;    
	if(bRbarDetail)
	{ if(Get_Rebar_DatilData(dXn, dRotation, ResData.arRbarRes, dAsThk)) ResData.nType = 3; }
	
	return TRUE;		
}

BOOL CCalcPmcv::Get_Rebar_DatilData(double dXn, double dRotation, CArray<T_PM_RBAR_RES, T_PM_RBAR_RES&>& arRbarRes, double dAsThk)
{
	arRbarRes.RemoveAll();
	
	double dCenter_y  = m_dYbar;
	double dCenter_z  = m_dZbar;
	double dYmax, dYmin, dZmax, dZmin;
	
	// 단면의 회전정보 생성
	// (회전계산은 CDgnCalcBase_Section_Tool::Get_RotateConPolyData, CDgnCalcBase_Section_Tool::Get_RotateRebarData과 일치시킬것) 
	
	// 최대 최소 초기화
  dYmax = dYmin = dCenter_y;
	dZmax = dZmin = dCenter_z;
	
  int nSize;
	double dx, dy, dz;
  // Concrete
	nSize = m_ConPolyData.OutPoly.aVertex.GetSize();	
	
	if(IsGeneralType())	
	{	
		_UMD_RC_GSEC_VERTEX PontUnit;
		for(int i=0; i<nSize; i++)
		{
			PontUnit = m_ConPolyData.OutPoly.aVertex.GetAt(i);
			dx = PontUnit.dX;
			dy = PontUnit.dY;
			dz = 0.0;
			
			// 중립축을 반시계 방향으로 돌릴려면 단면좌표를 시계방향으로 돌려야함
			CMathFunc::mathRotate(-dRotation, dCenter_y, dCenter_z, 0.0, 0.0, 0.0, 1.0, dx, dy, dz); 
			
			if(i==0)
			{// 첫번째 좌표값으로 최대 초소 초기화
				dYmax = dYmin = dx;
				dZmax = dZmin = dy;
			}
			else
			{// 두번째 좌표값부터 외각선의 최대 최소값을 찾음
				if(dx > dYmax)        dYmax = dx;
				if(dx < dYmin)        dYmin = dx;
				if(dy > dZmax)        dZmax = dy;
				if(dy < dZmin)        dZmin = dy;      
			}
		}
	}
	else 
	{
		if(!Calc_DistanceFromSectionToAngleLine(0, 0, dRotation, dZmax, dZmin)) return FALSE;
	}
	// 중립축 위치계산
	double dNAxisZ = dZmin + dXn*(dZmax-dZmin);
	double dDcomp  = dZmax - dNAxisZ;
	// Rebar
	nSize = m_RbarData.arRbarUnit.GetCount();
	arRbarRes.SetSize(nSize);
	
	_UMD_RC_RBAR_UNIT    RbarUnit;
	T_PM_RBAR_RES RbarRes;
	double dDist, desi;
	BOOL bComp;
	int iIndex=0;
	POSITION Pos = m_RbarData.arRbarUnit.GetStartPosition();
	int nCount = 0;
	while(Pos)
	{
		RbarUnit.Init();
		m_RbarData.arRbarUnit.GetNextAssoc(Pos, iIndex, RbarUnit);
		
		dx = RbarUnit.dyz[0];
		dy = RbarUnit.dyz[1];
    dz = 0.0;
		
    // 중립축을 반시계 방향으로 돌릴려면 단면좌표를 시계방향으로 돌려야함
    CMathFunc::mathRotate(-dRotation, dCenter_y, dCenter_z, 0.0, 0.0, 0.0, 1.0, dx, dy, dz); 
		
		dDist = fabs(dy-dNAxisZ);
		bComp = (dy >= dNAxisZ);
		desi  = Calc_esi(dDcomp, dDist);
		
    RbarRes.Init();
		RbarRes.dyz[0]  = RbarUnit.dyz[0];
		RbarRes.dyz[1]  = RbarUnit.dyz[1];
    RbarRes.dds     = dDist;   // Add Seungjun '20120203 for eGen
    RbarRes.desi    = desi; // Add Seungjun '20120203 for eGen
		if(m_bDesign) RbarRes.dAs = (RbarUnit.iDgnType==0) ? dAsThk*RbarUnit.dLen : RbarUnit.dDgnArea;
		else          RbarRes.dAs = RbarUnit.dArea;		
		RbarRes.dStress = Calc_Fsi(bComp, desi);
		RbarRes.dForce  = RbarRes.dAs * RbarRes.dStress;
		
		arRbarRes.SetAt(nCount, RbarRes);
		nCount++;
		if(nCount >= nSize)
			break;
	}	
	
	return TRUE;
}

BOOL CCalcPmcv::Is_ConcComp(int iIndexC, double dRota, _UMD_RC_CONC_UNIT& ConcUnit, double dDy2, double dDz2,
															 double& dCalcAc, double& dyInc, double& dzInc)
{	
	double dAc		= ConcUnit.dArea;
	double dyz[3] = {ConcUnit.dyz[0], ConcUnit.dyz[1], 0.0};
	double dyzPnt[CONST_UMDRC_iRECT][3]={0.0};
	for(int i=0; i<CONST_UMDRC_iRECT; i++)
	{
		dyzPnt[i][0] = ConcUnit.dyzPnt[i][0];
		dyzPnt[i][1] = ConcUnit.dyzPnt[i][1];
		dyzPnt[i][2] = 0.0;
	}
	// Check if Concrete Compression Region is OK by Offset Neutral axis.
	double dCoordZ[CONST_UMDRC_iRECT]={0.0};
  double dCoordY[CONST_UMDRC_iRECT]={0.0};
	BOOL bIsComp[CONST_UMDRC_iRECT]={FALSE};
	
	BOOL bIs00Deg = fabs(dRota) < cUMDRC_Zero;
	BOOL bIs90Deg = fabs(dRota-CMathFunc::m_pi/2.0) < cUMDRC_Zero;
	BOOL bAllComp = TRUE;
	if(bIs00Deg)	// 0 Deg.
	{
		for(int i=0; i<CONST_UMDRC_iRECT; i++)
		{
			dCoordZ[i] = dDz2;
			bIsComp[i] = (dyzPnt[i][1] > dCoordZ[i] ? TRUE : FALSE);			
		}
	}
	else if(bIs90Deg)	// 90 Deg.
	{
		for(int i=0; i<CONST_UMDRC_iRECT; i++)
		{
			dCoordY[i] = dDy2;
			bIsComp[i] = (dyzPnt[i][0] > dCoordY[i] ? TRUE : FALSE);
		}
	}
	else 
	{
		if(fabs(dDz2-1./(cUMDRC_Zero)) < cUMDRC_Zero)
		{
			for(int i=0; i<CONST_UMDRC_iRECT; i++)
			{
				dCoordY[i] = dDy2;
				bIsComp[i] = (dyzPnt[i][0] > dCoordY[i] ? TRUE : FALSE);
			}
		}
		else
		{
			for(int i=0; i<CONST_UMDRC_iRECT; i++)
			{
				dCoordZ[i] = (-1)*tan(dRota)*dyzPnt[i][0] + dDz2;
				bIsComp[i] = (dyzPnt[i][1] > dCoordZ[i] ? TRUE : FALSE);
			}
		}
	}
	
	dCalcAc	= 0.0;
	dyInc		= 0.0;
	dzInc		= 0.0;
	if(bIsComp[0] && bIsComp[1] && bIsComp[2] && bIsComp[3])	dCalcAc = dAc;
	else if(bIsComp[0] || bIsComp[1] || bIsComp[2] || bIsComp[3])
	{
    double dLineI[3]={dDy2,0.0,0.0};
    double dLineJ[3]={0.0,dDz2,0.0};
		if     (bIs00Deg) { dLineI[1] = dDz2; }
		else if(bIs90Deg) { dLineJ[0] = dDy2; }

		double dDistance=0.0;
		double dIntPnt[CONST_UMDRC_iRECT][3]={0.0};
		BOOL bExistPnt[CONST_UMDRC_iRECT]={FALSE};
		int iIntCount=0;
		for(int i=0; i<CONST_UMDRC_iRECT; i++)
		{
      for(int j=0; j<3; j++) dIntPnt[i][j] = 0.0;
      int iPrev = (i==0 ? CONST_UMDRC_iRECT-1 : i-1);
      int iNext = i;
			double dDistance=0.0;

			if(bIs00Deg && (dyzPnt[iPrev][1]==dDz2 || dyzPnt[iNext][1]==dDz2))
			{
				if(dyzPnt[iPrev][1]==dDz2 && dyzPnt[iNext][1]==dDz2) continue;
				bExistPnt[i] = TRUE;  iIntCount++;
				if(dyzPnt[iPrev][1]==dDz2) { for(int k=0; k<3; k++) dIntPnt[i][k] = dyzPnt[iPrev][k]; }
				else                       { for(int k=0; k<3; k++) dIntPnt[i][k] = dyzPnt[iNext][k]; }
			}
			else if(bIs90Deg && (dyzPnt[iPrev][0]==dDy2 || dyzPnt[iNext][0]==dDy2))
			{
				if(dyzPnt[iPrev][0]==dDy2 && dyzPnt[iNext][0]==dDy2) continue;
				bExistPnt[i] = TRUE;  iIntCount++;
				if(dyzPnt[iPrev][0]==dDy2) { for(int k=0; k<3; k++) dIntPnt[i][k] = dyzPnt[iPrev][k]; }
				else                       { for(int k=0; k<3; k++) dIntPnt[i][k] = dyzPnt[iNext][k]; }
			}
			else 
			{
				bExistPnt[i] = CMathFunc::mathIntersectLine2(dyzPnt[iPrev],dyzPnt[iNext],dLineI,dLineJ,cUMDRC_Zero,dDistance,dIntPnt[i]);
				if(bExistPnt[i] && dIntPnt[i][0]!=0.0 && dIntPnt[i][1]!=0.0)	iIntCount++;
			}
		}
		//if(iIntCount!=2)	ASSERT(0);	// Always 2.

		int iCompPntCount=0;
		for(int i=0; i<CONST_UMDRC_iRECT; i++)	{if(bIsComp[i])	iCompPntCount++;}

		int iCountPnt=0;
		while(iCountPnt < CONST_UMDRC_iRECT)
		{
			int iPrev = iCountPnt;
			int iNext = (iCountPnt < CONST_UMDRC_iRECT-1 ? iCountPnt+1 : 0);
			int iGoN1 = iPrev;
			int iGoN2=0;
			for(int k=0; k<CONST_UMDRC_iRECT; k++)	{if(k!=iPrev && bExistPnt[k])	iGoN2 = k;}
			if(iCompPntCount==1)	// Triangle.
			{
				if(bIsComp[iPrev])
				{
					double dyy[3] = {dyzPnt[iPrev][0], dIntPnt[iGoN1][0], dIntPnt[iGoN2][0]};
					double dzz[3] = {dyzPnt[iPrev][1], dIntPnt[iGoN1][1], dIntPnt[iGoN2][1]};
					double dyzCen[2]={0.0};
					CMathFunc::mathPolyCentroid(3,dyy,dzz,dyzCen[0],dyzCen[1],dCalcAc);
					dyInc	= dyzCen[0] - dyz[0];
					dzInc	= dyzCen[1] - dyz[1];
					break;
				}
			}
			else if(iCompPntCount==2)	// Rectangle.
			{
				if(bIsComp[iPrev] && bIsComp[iNext])
				{
					double dyy[4] = {dyzPnt[iPrev][0], dyzPnt[iNext][0], dIntPnt[iGoN2][0], dIntPnt[iGoN1][0]};
					double dzz[4] = {dyzPnt[iPrev][1], dyzPnt[iNext][1], dIntPnt[iGoN2][1], dIntPnt[iGoN1][1]};
					double dyzCen[2]={0.0};
					CMathFunc::mathPolyCentroid(4,dyy,dzz,dyzCen[0],dyzCen[1],dCalcAc);
					dyInc	= dyzCen[0] - dyz[0];
					dzInc	= dyzCen[1] - dyz[1];
					break;
				}
			}
			else if(iCompPntCount==3)	// Pentagon.
			{
			  int iPrev = iCountPnt;
        int iNext = (iCountPnt < CONST_UMDRC_iRECT-1 ? iCountPnt+1 : iPrev-3);
				int iMore = (iCountPnt < CONST_UMDRC_iRECT-2 ? iCountPnt+2 : iPrev-2);
				if(bIsComp[iPrev] && bIsComp[iNext] && bIsComp[iMore])
				{
					double dyy[5] = {dyzPnt[iPrev][0], dyzPnt[iNext][0], dyzPnt[iMore][0], dIntPnt[iGoN2][0], dIntPnt[iGoN1][0]};
					double dzz[5] = {dyzPnt[iPrev][1], dyzPnt[iNext][1], dyzPnt[iMore][1], dIntPnt[iGoN2][1], dIntPnt[iGoN1][1]};
					double dyzCen[2]={0.0};
					CMathFunc::mathPolyCentroid(5,dyy,dzz,dyzCen[0],dyzCen[1],dCalcAc);
					dyInc	= dyzCen[0] - dyz[0];
					dzInc	= dyzCen[1] - dyz[1];
					break;
				}
			}
			else	ASSERT(0);
			iCountPnt++;
		}
	}
//	if(fabs(dCalcAc) > dAc+cUMDRC_Zero)	ASSERT(0);
  dCalcAc = fabs(dCalcAc);
	return (fabs(dCalcAc) > 0.0);
}

double CCalcPmcv::Calc_esi(double dDcomp, double dDist) // 중립축에서 떨어진거리(dDist)에서의 단면의 변형률
{
  int iCode = m_iDgnCode;
  double decu = Get_ecu();//허용응력에 관련된 설계법은 지원하지 않습니다.(IRC_21_00지원하지 않음)
	// Get Required Data from MyDB.
  double desi = (dDcomp > 0.0 ? (dDist/dDcomp)*decu : 0.0);
	return desi;
}

double CCalcPmcv::Calc_Fsi(BOOL bComp, double desi)    // 중립축에서 떨어진거리(dDist)에서의 철근의 응력
{
	double dFsi=0.0;
	// Get Required Data from MyDB.
	double dEsr = m_dEsr;
	double dFyr = Get_ReBarStrength();
	//!!!if(CDgnBase_Tool::IsLSD(m_iDgnCode)) dFyr = Get_ft_Common_LSD(dFyr);
	//허용응력에 관련된 설계법은 지원하지 않습니다.(IRC_21_00지원하지 않음)
  if(bComp)	dFsi =      min(dFyr, dEsr*desi);	// Compression.
	else			dFsi = (-1)*min(dFyr, dEsr*desi);	// Tension.
  
	return dFsi;
}

double CCalcPmcv::Get_ecc()
{
	double decc=0.0;
	
  if(m_iTypeStressStrain==1) //Parabolic-Descending
		decc = 0.002;
	else if(m_iTypeStressStrain==2)	//Parabolic-Plateau
    decc = 2.0*Get_Sigcc()/m_dEc;
	else if(m_iTypeStressStrain==3 || m_iTypeStressStrain==4)	//Whitney(직사각형 응력블럭)
  {
    double decu = Get_ecu();
    decc = decu*(1.0-Calc_Beta());    
  }
  else if(m_iTypeStressStrain==5)	// user defined non-linear
  {
		ASSERT(0);
//     if(m_NonLinearProp.nConHysModelType==1)
//       decc = m_NonLinearProp.CON_CURVEC.dec1;
//     else if(m_NonLinearProp.nConHysModelType==2)
//       decc = m_NonLinearProp.CON_PAREEC.dec2;
//     else if(m_NonLinearProp.nConHysModelType==3)
//       decc = m_NonLinearProp.CON_BILNEC.dec3;
//     else if(m_NonLinearProp.nConHysModelType==4)      
//       decc = m_NonLinearProp.CON_KENTPK.dE0;
//     else if(m_NonLinearProp.nConHysModelType==5)
// 		{//Mander Model
// 			GSD_CALC_FIMP_CON_MANDER ConcD = m_NonLinearProp.CON_MANDER;
// 			if(ConcD.nConcType == 0)
// 				decc = ConcD.dConcDataeco;
// 			else 
// 				decc = ConcD.dConfinedConcStrainecc;
// 		}
//     else
//     {
//       double decu = Get_ecu();
//       decc = decu*(1.0-Calc_Beta());    
//     }
  }
	else	ASSERT(0);
	
	return decc;
}

BOOL CCalcPmcv::Calc_DistanceFromSectionToAngleLine(double dx, double dy, double dAngle, double& dDistMax, double& dDistMin)
{// SHIN ('06.03.15) 추가 
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
  case DGN_SECT_SHAPE_INDEX_REG_B: 
  case DGN_SECT_SHAPE_INDEX_REG_SB:
    dHc = m_dSize[0];
		dB1 = m_dSize[1];    
		// Point(4).
		pPointData[0].id = 1;	pPointData[0].uv[0] = 0.0;	pPointData[0].uv[1] = 0.0;
		pPointData[1].id = 2;	pPointData[1].uv[0] = dB1;	pPointData[1].uv[1] = 0.0;
		pPointData[2].id = 3;	pPointData[2].uv[0] = dB1;	pPointData[2].uv[1] = dHc;
		pPointData[3].id = 4;	pPointData[3].uv[0] = 0.0;	pPointData[3].uv[1] = dHc;    
    break;
  case DGN_SECT_SHAPE_INDEX_REG_P:
  case DGN_SECT_SHAPE_INDEX_REG_SR:
    dHc = m_dSize[0];
    // CirclePoint(1).
    pCirclePointData[0].id = 1;	pCirclePointData[0].uv[0] = dHc/2.0;	pCirclePointData[0].uv[1] = dHc/2.0;
    dCircleRadius[0] = dHc/2.0;
    break;
  case DGN_SECT_SHAPE_INDEX_REG_OCT:
  case DGN_SECT_SHAPE_INDEX_REG_SOCT:
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
  case DGN_SECT_SHAPE_INDEX_REG_TRK:
  case DGN_SECT_SHAPE_INDEX_REG_STRK:
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
  case DGN_SECT_SHAPE_INDEX_REG_HTRK: // Half Track        
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
  return TRUE;
}

BOOL CCalcPmcv::Cal_MaxDistRbar_Gen(double dXn, double dDmax, double dZmin, _UMD_RC_RBAR_LIST& arRbarUnit, double& dMaxTens, double& dMaxComp, BOOL bFixDgnRebarOnly)
{
  // Calculate maximum distance from Rebar to neutral axis.
	dMaxTens = 0.0;// 중립축에서 가장 멀리떨어진 인장철근의 거리
	dMaxComp = 0.0;// 중립축에서 가장 멀리떨어진 압축철근의 거리
	
	_UMD_RC_RBAR_UNIT RbarUnit;
	int i=0;
	
	for(int i=0; i<arRbarUnit.List.GetSize(); i++)
	{
		RbarUnit.Init();
		RbarUnit = arRbarUnit.List.GetAt(i);
		if(bFixDgnRebarOnly)
		{	if(m_bDesign && RbarUnit.iDgnType==0) continue; }
		// Check if compression or not.
		double dDz1	 = dDmax*dXn;// 
		double dDist = dZmin+dDz1  - RbarUnit.dyz[1];
		if(dDist > 0.0)	dMaxTens = max(dMaxTens, dDist);	// Tension.
		else						dMaxComp = min(dMaxComp, dDist);	// Compression.
	}
	// Absolute Distance.
	dMaxTens = fabs(dMaxTens);
	dMaxComp = fabs(dMaxComp);
	
	if(dMaxTens==0.0 && dMaxComp==0.0)	return FALSE;	// For Checking.
	return TRUE;
}
BOOL CCalcPmcv::Cal_MaxDistRbar_Gen(double dXn, int iAxisDir, double dDmax, double& dMaxTens, double& dMaxComp, BOOL bFixDgnRebarOnly)
{
	// Calculate maximum distance from Rebar to neutral axis.
	dMaxTens = 0.0;// 중립축에서 가장 멀리떨어진 인장철근의 거리
	dMaxComp = 0.0;// 중립축에서 가장 멀리떨어진 압축철근의 거리
	
	_UMD_RC_CON_PART_LIST arConcUnit;
	_UMD_RC_RBAR_LIST arRbarUnit;
	
	if(iAxisDir >= Get_HoldNumberToPmcvData())	{ASSERT(0); return FALSE;} 
	
	if(!Cng_AxisDir_Gen(iAxisDir, arConcUnit, arRbarUnit))	return FALSE;  	
  
  return Cal_MaxDistRbar_Gen(dXn, dDmax, arConcUnit.dZmin, arRbarUnit, dMaxTens, dMaxComp, bFixDgnRebarOnly);  
}

double CCalcPmcv::Get_Sigcc()
{   
  double dSigcc=0.0;
//   if(m_iTypeStressStrain==5 && m_NonLinearProp.nConHysModelType>0)  // user defined non-linear.
//   {
//     if(m_NonLinearProp.nConHysModelType==1)
//       dSigcc = m_NonLinearProp.CON_CURVEC.dfck;
//     else if(m_NonLinearProp.nConHysModelType==2)
//       dSigcc = m_NonLinearProp.CON_PAREEC.dfck;
//     else if(m_NonLinearProp.nConHysModelType==3)
//       dSigcc = m_NonLinearProp.CON_BILNEC.dfck;
//     else if(m_NonLinearProp.nConHysModelType==4)      
//       dSigcc = m_NonLinearProp.CON_KENTPK.dFc;
//     else if(m_NonLinearProp.nConHysModelType==5)
// 		{//Mander Model
// 			GSD_CALC_FIMP_CON_MANDER ConcD = m_NonLinearProp.CON_MANDER;
// 			if(ConcD.nConcType == 0)
// 				dSigcc = ConcD.dConcDataFco;
// 			else 
// 				dSigcc = ConcD.dConfinedConcStrengthFcc;
// 		}
//     else ASSERT(0);
// 		
//     return dSigcc;
//   }
// 	
  dSigcc = Calc_Alpha()*m_dFc;    	
	return dSigcc;
}

double CCalcPmcv::Get_ery()
{
	return (m_dEsr==0.0 ? 0.0 : m_dFyr/m_dEsr);
}

double CCalcPmcv::Calc_SigrU(double dDcomp, double dDist, double dMaxDistCompRbar, double& der)
{
	// Rbar Stress (Ultimate).
	double dSigr=0.0;
	double decu = Get_ecu();
	double desy	= Get_ery();
	double dSign= (dDist > 0.0 ? 1.0 : -1.0);	// Tens(-), Comp(+).
	double dDedge = dDcomp;
	der = (dDedge==0.0 ? dSign*decu : decu*dDist/dDedge);  
	
  double derTemp;
  double dSigcc = Calc_SigcU(dDcomp, dDist, dMaxDistCompRbar, derTemp);//double dSigcc= Cal_Sigcc();
  dSigr = min(m_dFyr, m_dEsr*fabs(der));

	dSigr *= (der > 0.0 ? +1.0 : -1.0);	// If der > 0, Comp.

	return dSigr;
}

double CCalcPmcv::Get_Db()
{
  double dDb=1.0;
	
  if(m_iTypeStressStrain==1)	//Parabolic-Descending
		dDb = 0.85;
	else if(m_iTypeStressStrain==2)	//Parabolic-Plateau
		dDb = 1.0;		
	else if(m_iTypeStressStrain==3 || m_iTypeStressStrain==4)	//Whitney(직사각형 응력블럭)
		dDb = 1.0;		
  else if(m_iTypeStressStrain==5)	// user defined non-linear
		dDb = 1.0;		
	else	ASSERT(0);
	
	return dDb;
}


double CCalcPmcv::Calc_SigcU(double dDcomp, double dDist, double dMaxDistCompRbar, double& dec)
{
	// Concrete Stress.
	double dSigc =0.0;
	double decc	 = Get_ecc();
	double decu  = Get_ecu();
	double dSigcc= Get_Sigcc();
  double dDb   = Get_Db(); 
	double dDedge= dDcomp;
	dec = (dDedge==0.0 ? -decu : decu*dDist/dDedge);//해당위치에서의 변형률

  if(m_iTypeStressStrain==1)	// Parabolic-Descending
	{
    if(dec <= 0.0 || dec > decu) dSigc = 0.0;
    else if(dec >= decc)         dSigc = dSigcc*(1.0 - dDb*(dec-decc)/(decu-decc));
    else                         dSigc = dSigcc*(2*dec/decc-(dec/decc)*(dec/decc));		
	}
	else if(m_iTypeStressStrain==2)	// Parabolic-Plateau
	{
    if(dec <= 0.0 || dec > decu) dSigc = 0.0;
    else if(dec >= decc)         dSigc = dSigcc;
    else                         dSigc = dSigcc*(2*dec/decc-(dec/decc)*(dec/decc));		
	}
	else if(m_iTypeStressStrain==3 || m_iTypeStressStrain==4)	// Whitney(직사각형 응력블럭)
	{
		if(dec < decc || dec > decu) dSigc = 0.0;
    else                         dSigc = dSigcc;
	}
  else if(m_iTypeStressStrain==5)  // user defined non-linear.
  {
		ASSERT(0);
//     if(m_NonLinearProp.nConHysModelType==1)
//     {
//       dSigc = m_NonLinearProp.CON_CURVEC.Get_Sigc(dec);
//       dSigc = Get_fc_Common_LSD(dSigc);	// Concrete Design Strength.
//     }
//     else if(m_NonLinearProp.nConHysModelType==2)
//       dSigc = m_NonLinearProp.CON_PAREEC.Get_Sigc(dec);
//     else if(m_NonLinearProp.nConHysModelType==3)
//       dSigc = m_NonLinearProp.CON_BILNEC.Get_Sigc(dec);
//     else if(m_NonLinearProp.nConHysModelType==4)      
//       dSigc = m_NonLinearProp.CON_KENTPK.Get_Sigc(dec);
//     else if(m_NonLinearProp.nConHysModelType==5)
// 		{//Mander Model
// 			GSD_CALC_FIMP_CON_MANDER ConcD = m_NonLinearProp.CON_MANDER;
// 			if(dec > 0.0)
// 			{
// 				double dfc = 0.0;
// 				double dx  = 0.0;
// 				double dEsec = 0.0;				 
// 				if(ConcD.nConcType==0)
// 				{// Unconfined
// 					if(ConcD.dConcDataeco>0.0)
// 					{					
// 						dfc = ConcD.dConcDataFco;
// 						dx = dec/ConcD.dConcDataeco;
// 						dEsec = ConcD.dConcDataFco/ConcD.dConcDataeco;
// 					}
// 				}
// 				else 
// 				{// Confined
// 					if(ConcD.dConfinedConcStrainecc>0.0)
// 					{					
// 						dfc = ConcD.dConfinedConcStrengthFcc;
// 						dx = dec/ConcD.dConfinedConcStrainecc;
// 						dEsec = ConcD.dConfinedConcStrengthFcc/ConcD.dConfinedConcStrainecc;
// 					}
// 				}
// 				double dr = (ConcD.dConcDataEc-dEsec)==0.0 ? 0.0 : ConcD.dConcDataEc/(ConcD.dConcDataEc-dEsec);
// 				dSigc = (dr-1.0+pow(dx, dr))==0.0 ? 0.0 : dfc*dx*dr/(dr-1.0+pow(dx, dr));
// 			}
// 			else
// 			{
// 				if(dec >= (-1.0)*ConcD.dTensConcDataet)
// 				{
// 					dSigc = (ConcD.dTensConcDataet==0.0) ? 0.0 : dec * ConcD.dTensConcDataft/ConcD.dTensConcDataet;					
// 				}
// 			}
// 		}
//     else  // Whitney(직사각형 응력블럭)
//     {
//       if(dec < decc || dec > decu) dSigc = 0.0;
//       else                         dSigc = dSigcc;
//     }
  }
	else	
		ASSERT(0);

	return dSigc;
}

BOOL CCalcPmcv::Cng_AxisDir_Gen(int iAxisDir, _UMD_RC_CON_PART_LIST& arConUnit, _UMD_RC_RBAR_LIST& arRbarUnit)
{
	arConUnit.Init();
	arRbarUnit.Init();
	
	double dyCenter = m_dYbar;
	double dzCenter = m_dZbar;
	
	// iAxisDir (0:0Deg,  1~:iAxisDir*90/m_divisionNumber_90deg)
	if(iAxisDir < 0 || iAxisDir >= Get_HoldNumberToPmcvData())	{ASSERT(0); return FALSE;}
	int iAreaKind = iAxisDir;
  if(m_iSymmetryType != 2 && iAxisDir >= m_iDivisionNumber_90deg*2)
    iAreaKind = iAxisDir-m_iDivisionNumber_90deg*2;// m_arConUnit및 m_arRbarUnit의 해당 행
	
	// SetSize for Speed Up.
  arConUnit  = m_RotSectData.arConUnitList.GetAt(iAreaKind);  
  arRbarUnit = m_RotSectData.arRbarUnitList.GetAt(iAreaKind);
		
	// Change Direction.
	BOOL byzCng = ((m_RotSectData.iSymmetryType != 2 && iAxisDir >= m_RotSectData.iDivisionNumber_90deg*2) ? TRUE : FALSE);	// Reverse Data about Axis.
	
	int iNumCon	 = arConUnit.List.GetSize();
	int iNumRbar = arRbarUnit.List.GetSize();
	if(byzCng)
	{		
		int i=0;
		// Concrete.
		for(int i=0; i<iNumCon; i++)
		{
			_UMD_RC_CON_PART_UNIT ConUnit = arConUnit.List.GetAt(i);
			if(byzCng)	ConUnit.dyz[0] = 2.0 * dyCenter - ConUnit.dyz[0];
			if(byzCng)	ConUnit.dyz[1] = 2.0 * dzCenter - ConUnit.dyz[1];
			arConUnit.List.SetAt(i,ConUnit);
		}
		// Rbar.
		for(int i=0; i<iNumRbar; i++)
		{
			_UMD_RC_RBAR_UNIT RbarUnit = arRbarUnit.List.GetAt(i);
			if(byzCng)	RbarUnit.dyz[0] = 2.0 * dyCenter - RbarUnit.dyz[0];
			if(byzCng)	RbarUnit.dyz[1] = 2.0 * dzCenter - RbarUnit.dyz[1];
			arRbarUnit.List.SetAt(i,RbarUnit);
		}	
		double dYmax = arConUnit.dYmax;
		double dZmax = arConUnit.dZmax;
		double dYmin = arConUnit.dYmin;
		double dZmin = arConUnit.dZmin;
		
		arConUnit.dYmin = 2.0 * dyCenter - dYmax;
		arConUnit.dZmin = 2.0 * dzCenter - dZmax;
		arConUnit.dYmax = 2.0 * dyCenter - dYmin;
		arConUnit.dZmax = 2.0 * dzCenter - dZmin;
	}
	
	if(m_iTypeStressStrain==3)
	{
		if(iNumRbar==0)
		{
			arConUnit.Init();
			arRbarUnit.Init();
			return FALSE;
		}
	}
	else
	{
		if(iNumCon*iNumRbar==0)
		{
			arConUnit.Init();
			arRbarUnit.Init();
			return FALSE;
		}
	}
	
	return TRUE;
}


double CCalcPmcv::Get_RotationAngle(int iAxisDir)
{
  // SHIN (06.01.27) 가로 세로가 차이 날경우 해석결과가 한쪽을 치우칠수 있으므로 보완방법 추가 요망
  if(m_iSymmetryType == 4)// 상하대칭일때에는 90도 부터 시작
    return double(iAxisDir)*90.0/double(m_iDivisionNumber_90deg) + 90.0;
  else                    // 나머지의 경우 0도 부터 시작
    return double(iAxisDir)*90.0/double(m_iDivisionNumber_90deg);
}


double CCalcPmcv::Calc_Phi(double dPn, double dPb, double desi)
{
	double dPhi = 0.0;
	//!!! Need to override.
	double dFyr = Get_ReBarStrength();
	double dEpsi_y = dFyr/m_dEsr;
	double dEpsi_t = dFyr<(400.+cUMDRC_Zero) ? 0.005 : 2.5*dEpsi_y; // Modify, Jaeoh. [8/11/2011] 0.004 -> 0.005
	int iHoopType = m_iHoopType;
	
	double dPhiUp = m_dPhi[1];		// Tension.
	double dPhiDn = (iHoopType==2 ? m_dPhi[2] : m_dPhi[3]);

	if(desi <= dEpsi_y)			dPhi = (dPn > 0.0 ? dPhiDn : dPhiUp);	// Distinguish Pure Tension.
	else if(desi >= dEpsi_t)	dPhi = dPhiUp;
	else
	{
		if(iHoopType==2)	dPhi = 0.70 + 0.15*(desi-dEpsi_y)/(dEpsi_t-dEpsi_y);	// Spiral
		else				      dPhi = 0.65 + 0.20*(desi-dEpsi_y)/(dEpsi_t-dEpsi_y);	// Tied
	}
	return dPhi;
}

double CCalcPmcv::Calc_Beta()
{
	//!!! need to override.
	double dFck = m_dFc;
	return (dFck<=  28.0 ? 0.85 : max(0.85-0.007*(dFck-28.0)/1.0, 0.65));
}

double CCalcPmcv::Get_BalancedXb(double dDmax, double dDeff)
{
	double decu   = Get_ecu();
	double dEs    = m_dEsr;
	double dFyr   = Get_ReBarStrength();
  return (dDmax==0.0 ? 0.0 : 1.0 - (dEs*decu/(dEs*decu+dFyr)*dDeff)/dDmax);// dXb값은 인장측에서 중립축 까지의 거리를 나타내는 비임;
}

double CCalcPmcv::Get_PhiPnmax(T_PMCV_2D& PmData)
{
	double dPhiPnmax=0.0;
	
	if(PmData.arPmUnit.GetSize() > 0)
	{
    for(int i=0; i<PmData.arPmUnit.GetSize(); i++)
    {
			T_PM_UNIT PmUnit;
			PmUnit = PmData.arPmUnit.GetAt(i);
      //(2012.02.28_PKS) DB단면과 Gen단면에서 최대 축력일때를 다르게 정의 : Is_PureComp()
      if(Is_PureComp(PmUnit.dXn))
      {
				double dFactor = Get_MaxPnReductionFactor();
				dPhiPnmax = dFactor * PmUnit.dPhi * PmUnit.dPn;	// Assume Tie.
        break;
      }
    }
	}
	
	return dPhiPnmax;
}

int CCalcPmcv::Get_HoldNumberToPmcvData()
{
  // ※참고문헌(기하 알고리즘 P1118)
  switch(m_iSymmetryType)
  {
  case 1:
    return m_iDivisionNumber_90deg*4;
  case 2:// 좌우상하대칭시에는 0-90까지만
    return m_iDivisionNumber_90deg + 1;// +1은 90부분을 포함하기 위해서임
  default:// 상하, 좌우, 원점대칭시에는 0-180또는 90-270까지만
    return m_iDivisionNumber_90deg*2 + 1;// +1은 180부분을 포함하기 위해서임
  }
}

CString CCalcPmcv::GetFileName()
{
	return _T("PMCurve.txt");
}