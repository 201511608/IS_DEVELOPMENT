// MdgnPmTool.cpp: implementation of the CMdgnPmTool class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Pmcv.h"
#include "MdgnPmTool.h"
#include "MathFunc.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMdgnPmTool::CMdgnPmTool()
{

}

CMdgnPmTool::~CMdgnPmTool()
{

}

void CMdgnPmTool::Set_Tolerance(double dTolAngle, double dTolM0, double dTolP0)
{
	m_dTolAngle = dTolAngle;
	m_dTolM0    = dTolM0;
	m_dTolP0    = dTolP0;
}

BOOL CMdgnPmTool::Get_ConvertToPMTool(int nType, T_PMCV_3D& InData, _UMD_3DPM_DATA& OutData)
{
	int nSize = InData.arPmcvData.GetSize();	
	//정보가 없으면 종료(현재(07.11.02) DB단면으로 계산시 P0값을 만들지 않고 있음)
	if(InData.bCbData && nSize != InData.arCbPmData.GetSize()) { ASSERT(0); return FALSE; }
	if(InData.bP0Data && nSize != InData.arP0PmData.GetSize()) { ASSERT(0); return FALSE; }
	if(!InData.b3DData || !InData.bCbData                    ) { ASSERT(0); return FALSE; }
	
	OutData.Init();
	OutData.nType = nType;
	OutData.iSymmetryType         = InData.iSymmetryType;
	OutData.iDivisionNumber_90deg = InData.iDivisionNumber_90deg;
	OutData.PmcvData.SetSize(nSize);
	
	_UMD_PMCV_DATA_POLY PmcvPoly;
	T_PMCV_2D          PmcvData;
	_UMD_PMCV_UNIT_POLY PmcvPolyUnit;
	T_PM_UNIT     PmcvDataUnit;
	T_PM_SECT_RES CbPmData, P0PmData;
	int nUnitSize;
	for(int i=0 ; i<nSize ; i++)
	{
		PmcvData = InData.arPmcvData.GetAt(i);
		CbPmData = InData.arCbPmData.GetAt(i);
		if(InData.bP0Data) P0PmData = InData.arP0PmData.GetAt(i);
		else               P0PmData.Init();
		
		PmcvPoly.Init();
		PmcvPoly.dRotate = PmcvData.dRotate * CMathFunc::m_trang;    // 소성중심축의 회전각 소성중심축이 y축방향이고 압축부가 Z축'+'방향일때 0으로 기준 방향은 반시계방향 (Deg)
		PmcvPoly.dXb     = CbPmData.dXn;  
		PmcvPoly.dPb     = CbPmData.dPn;  
		PmcvPoly.dMby    = CbPmData.dMny;
		PmcvPoly.dMbz    = CbPmData.dMnz;	
		PmcvPoly.desib   = CbPmData.desiMax;
		PmcvPoly.dX0     = P0PmData.dXn;  
		PmcvPoly.dM0y    = P0PmData.dMny;
		PmcvPoly.dM0z    = P0PmData.dMnz; 
		PmcvPoly.dDmax   = PmcvData.dDmax;
		PmcvPoly.dDeff   = PmcvData.dDeff;
		
		nUnitSize = PmcvData.arPmUnit.GetSize(); 	
		if(nType == 1)
		{
			for(int j=0 ; j<nUnitSize ; j++)
			{
				PmcvDataUnit = PmcvData.arPmUnit.GetAt(j); 		
				PmcvPolyUnit.Init();
				PmcvPolyUnit.dPn  = PmcvDataUnit.dPn;
				PmcvPolyUnit.dMny = PmcvDataUnit.dMny;
				PmcvPolyUnit.dMnz = PmcvDataUnit.dMnz;
				PmcvPolyUnit.dXn  = PmcvDataUnit.dXn;
				PmcvPolyUnit.desi = PmcvDataUnit.desi;
				PmcvPoly.arPmcvUnit.SetAt(j+1, PmcvPolyUnit);
			}		
		}
		else
		{
			for(int j=0 ; j<nUnitSize ; j++)
			{
				PmcvDataUnit = PmcvData.arPmUnit.GetAt(j); 		
				PmcvPolyUnit.Init();
				PmcvPolyUnit.dPn  = PmcvDataUnit.dPn  * PmcvDataUnit.dPhi;
				PmcvPolyUnit.dMny = PmcvDataUnit.dMny * PmcvDataUnit.dPhi;
				PmcvPolyUnit.dMnz = PmcvDataUnit.dMnz * PmcvDataUnit.dPhi;
				PmcvPolyUnit.dXn  = PmcvDataUnit.dXn;
				PmcvPolyUnit.desi = PmcvDataUnit.desi;
				PmcvPoly.arPmcvUnit.SetAt(j+1, PmcvPolyUnit);
			}		
		}
		
		OutData.PmcvData.SetAt(i, PmcvPoly);
	}
	return TRUE;
}

BOOL CMdgnPmTool::Get_PmcvUnitAt3DPmcvData(_UMD_3DPM_DATA& Umd3DPmData, double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, BOOL& bCheck, int& iLeftFitDirID, double& dModul_DirID)
{
	BOOL bResInStartPMM;
	return Get_PmcvUnitAt3DPmcvData(Umd3DPmData, 0.0, 0.0, 0.0, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, bCheck, iLeftFitDirID, dModul_DirID, FALSE, bResInStartPMM);
}


BOOL CMdgnPmTool::Get_PmcvUnitAt3DPmcvData(_UMD_3DPM_DATA& Umd3DPmData, double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, BOOL& bCheck, int& iLeftFitDirID, double& dModul_DirID, BOOL bChkInStartPMM, BOOL& bResInStartPMM)
{ 
	bResInStartPMM = TRUE;// 시작점이 PM상관도 내에 있는지 여부를 검토
	if(bChkInStartPMM)
	{
		if(dPs==0.0 && (dMs_y==0.0 && dMs_z==0.0))
			bResInStartPMM = TRUE;
		else
		{
			BOOL bResInStartPMM_Temp;
			int iLeftFitDirID_Temp;
			double dPn_Temp, dMn_y_Temp, dMn_z_Temp, desi_Temp, dXn_Temp, dModul_DirID_Temp;
			if(!Get_PmcvUnitAt3DPmcvData(Umd3DPmData, 0.0, 0.0, 0.0, dPs, dMs_y, dMs_z, dPn_Temp, dMn_y_Temp, dMn_z_Temp, desi_Temp, dXn_Temp, bResInStartPMM, iLeftFitDirID_Temp, dModul_DirID_Temp, FALSE, bResInStartPMM_Temp))
				return FALSE;
		}
	}
	int iDivisionNumber_90deg = Umd3DPmData.iDivisionNumber_90deg;
	if(Umd3DPmData.PmcvData.GetSize() < iDivisionNumber_90deg + 1)return FALSE;

	int iFitDirID, iFitUnitID;
  // 1단계 : 속도 증진을 위해 가장 해당 하중과 비슷한 방향의 점을 찾아 주변 평면을 검색함 ///////////////////////// 
	//         (시작점이 3D-PM내부에 있을 경우에만 적용)
	if(bResInStartPMM)
	{		
		// 해당 하중에 대하여 가장 가까운 방향의 PM좌표를 찾음
		Get_FitPmcvUnitID(Umd3DPmData, dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, iFitDirID, iFitUnitID);
		// 검색된 좌표주변 평면을 검색하여 교차점을 찾음
		CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&> arPmcvUnit;
		if(Get_PiercePointForContactPlane(Umd3DPmData, dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, iFitDirID, iFitUnitID, iLeftFitDirID, dModul_DirID))
		{    
			bCheck = (CMathFunc::mathLength(dPd-dPs, dMd_y-dMs_y, dMd_z-dMs_z) <= CMathFunc::mathLength(dPn-dPs, dMn_y-dMs_y, dMn_z-dMs_z));
			return TRUE;
		}
	}
  // 2단계 : 1단계 실패시에는 모든 평면을 검토하면서 교차점을 찾음 ////////////////////////////////////////////////
	if(bResInStartPMM)
  {// 시작점이 내부에 있는 경우
    int iPmcvDataCount = 4*iDivisionNumber_90deg;
    int iDirID;

    for(int i=0 ; i<iPmcvDataCount ; i++)
    {
      if(i%2==0)
        iDirID = iFitDirID + i/2;
      else
        iDirID = iFitDirID - 1 - (i-1)/2;

      if(Get_PiercePontForDirection(Umd3DPmData, dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, iDirID, iDirID+1, dModul_DirID))
      {    
        iLeftFitDirID = iDirID;
        bCheck = (CMathFunc::mathLength(dPd-dPs, dMd_y-dMs_y, dMd_z-dMs_z) <= CMathFunc::mathLength(dPn-dPs, dMn_y-dMs_y, dMn_z-dMs_z));
        return TRUE;
      }
    }
  }  
	else
	{// 시작점이 외부에 있는 경우
		int iPmcvDataCount = 4*iDivisionNumber_90deg;
    int iDirID;
		
		BOOL bCheck1, bCheck2;
		int  iLeftFitDirID1, iLeftFitDirID2;
		double dPn1, dMn_y1, dMn_z1, dModul_DirID1, dMsMn1;
		double dPn2, dMn_y2, dMn_z2, dModul_DirID2, dMsMn2;

		double dMsMd = CMathFunc::mathLength(dPd-dPs, dMd_y-dMs_y, dMd_z-dMs_z);


		int nCount = 0;
    for(int i=0 ; i<iPmcvDataCount ; i++)
    {
      if(i%2==0)
        iDirID = iFitDirID + i/2;
      else
        iDirID = iFitDirID - 1 - (i-1)/2;

      if(Get_PiercePontForDirection(Umd3DPmData, dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, iDirID, iDirID+1, dModul_DirID))
      {    
				if(nCount == 0)
				{
					dPn1   = dPn;
					dMn_y1 = dMn_y;
					dMn_z1 = dMn_z;
					dMsMn1 = CMathFunc::mathLength(dPn-dPs, dMn_y-dMs_y, dMn_z-dMs_z);
					dModul_DirID1  = dModul_DirID;
					iLeftFitDirID1 = iDirID;
					bCheck1 = (dMsMd <= dMsMn1);
					nCount++;
				}
				else if(nCount == 1)
				{
					dPn2   = dPn;
					dMn_y2 = dMn_y;
					dMn_z2 = dMn_z;
					dMsMn2 = CMathFunc::mathLength(dPn-dPs, dMn_y-dMs_y, dMn_z-dMs_z);
					dModul_DirID2  = dModul_DirID;
					iLeftFitDirID2 = iDirID;
					bCheck2 = (dMsMd <= dMsMn2);
					nCount++;
				}
				else
					break;  
      }
    }

		if(nCount == 1)
		{
			dPn   = dPn1;
			dMn_y = dMn_y1;
			dMn_z = dMn_z1;
			dModul_DirID  = dModul_DirID1;
			iLeftFitDirID = iLeftFitDirID1;
			bCheck = (dPn==dPd && (dMn_y==dMd_y && dMn_z==dMd_z));
			return TRUE;
		}
		else if(nCount == 2)
		{
			if(bCheck1 == bCheck1)
			{//2점 모두 dMsMd보다 크거나 작은 경우
				bCheck = FALSE;				
				if((bCheck1 && dMsMn1 >= dMsMn2) || (!bCheck1 && dMsMn1 <= dMsMn2))
				{//모두 큰 경우 둘 중 작은 값을, 모두 작은 경우  둘 중 큰 값을
					dPn   = dPn2;
					dMn_y = dMn_y2;
					dMn_z = dMn_z2;
					dModul_DirID  = dModul_DirID2;
					iLeftFitDirID = iLeftFitDirID2;						
				}
				else
				{
					dPn   = dPn1;
					dMn_y = dMn_y1;
					dMn_z = dMn_z1;
					dModul_DirID  = dModul_DirID1;
					iLeftFitDirID = iLeftFitDirID1;
				}
			}
			else
			{
				bCheck = TRUE;				
				if(dMsMn1 >= dMsMn2)
				{//dMsMd이 내부에 있으므로 큰값을 입력하여줌
					dPn   = dPn1;
					dMn_y = dMn_y1;
					dMn_z = dMn_z1;
					dModul_DirID  = dModul_DirID1;
					iLeftFitDirID = iLeftFitDirID1;				
				}
				else
				{
					dPn   = dPn2;
					dMn_y = dMn_y2;
					dMn_z = dMn_z2;
					dModul_DirID  = dModul_DirID2;
					iLeftFitDirID = iLeftFitDirID2;						
				}
			}
			return TRUE;
		}		
	}
  return FALSE;
}

BOOL CMdgnPmTool::Get_PmcvDataToDirID(_UMD_3DPM_DATA& Umd3DPmData, int iAxisDir_Global, double dModul_DirID, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcv, _UMD_PMCV_UNIT_BALANCE& PmBalance)
{
	int iDivisionNumber_90deg = Umd3DPmData.iDivisionNumber_90deg;
	int iSymmetryType = Umd3DPmData.iSymmetryType;
	if(Umd3DPmData.PmcvData.GetSize() < iDivisionNumber_90deg + 1)return FALSE;
  int iPmcvDataPrevID, iPmcvDataNextID, iPmcvSize;
  double dModul_PrevMy,   dModul_PrevMz,   dModul_NextMy,   dModul_NextMz;
  double dModul_PrevDmax, dModul_PrevDeff, dModul_NextDmax, dModul_NextDeff;
  double dModul_PrevAngle,dModul_PrevXb,   dModul_NextAngle,dModul_NextXb;
  
  _UMD_PMCV_DATA_POLY pmcvDataPrev, pmcvDataNext;
  _UMD_PMCV_UNIT_POLY pmcvUnit, pmcvUnitPrev, pmcvUnitNext;
	
  Get_ResultPmcvDataID(Umd3DPmData, iAxisDir_Global,   iPmcvDataPrevID, dModul_PrevMy, dModul_PrevMz);
  Get_ResultPmcvDataID(Umd3DPmData, iAxisDir_Global+1, iPmcvDataNextID, dModul_NextMy, dModul_NextMz);
  pmcvDataPrev = Umd3DPmData.PmcvData.GetAt(iPmcvDataPrevID);
  pmcvDataNext = Umd3DPmData.PmcvData.GetAt(iPmcvDataNextID);
  
  iPmcvSize = pmcvDataPrev.arPmcvUnit.GetCount();  
  
  arPmcv.RemoveAll();
	PmBalance.Init();
  
  arPmcv.SetSize(iPmcvSize);
	
  for(int i=1 ; i<=iPmcvSize ; i++)
  {
    pmcvUnitPrev = Get_CngSignPmcvUnit(pmcvDataPrev.arPmcvUnit[i], dModul_PrevMy, dModul_PrevMz);
    pmcvUnitNext = Get_CngSignPmcvUnit(pmcvDataNext.arPmcvUnit[i], dModul_NextMy, dModul_NextMz);
		
    pmcvUnit.Init();
    pmcvUnit.dPn  = pmcvUnitPrev.dPn +  (pmcvUnitNext.dPn  - pmcvUnitPrev.dPn )*dModul_DirID;
    pmcvUnit.dMny = pmcvUnitPrev.dMny + (pmcvUnitNext.dMny - pmcvUnitPrev.dMny)*dModul_DirID;
    pmcvUnit.dMnz = pmcvUnitPrev.dMnz + (pmcvUnitNext.dMnz - pmcvUnitPrev.dMnz)*dModul_DirID;
    pmcvUnit.dXn  = pmcvUnitPrev.dXn +  (pmcvUnitNext.dXn  - pmcvUnitPrev.dXn )*dModul_DirID;    
    pmcvUnit.desi = pmcvUnitPrev.desi +  (pmcvUnitNext.desi  - pmcvUnitPrev.desi )*dModul_DirID;
    arPmcv.SetAt(i-1, pmcvUnit);
  }
  dModul_PrevDmax = pmcvDataPrev.dDmax;
  dModul_NextDmax = pmcvDataNext.dDmax;
  dModul_PrevDeff = pmcvDataPrev.dDeff;
  dModul_NextDeff = pmcvDataNext.dDeff;
  dModul_PrevAngle= Get_RotationGlobalAngle(iDivisionNumber_90deg, iAxisDir_Global);
  dModul_NextAngle= Get_RotationGlobalAngle(iDivisionNumber_90deg, iAxisDir_Global+1);
	
  dModul_PrevXb   = pmcvDataPrev.dXb;
  dModul_NextXb   = pmcvDataNext.dXb;
  double dDmax = dModul_PrevDmax + (dModul_NextDmax - dModul_PrevDmax )*dModul_DirID;  
	
  PmBalance.dXb    = dModul_PrevXb   + (dModul_NextXb   - dModul_PrevXb   )*dModul_DirID;
  PmBalance.dDeffb = dModul_PrevDeff + (dModul_NextDeff - dModul_PrevDeff )*dModul_DirID;
	PmBalance.dDmax  = dDmax;
  PmBalance.dCb    = (1.0-PmBalance.dXb)*dDmax;//600.0/(600.0 + m_dFyr) * dDeffb;
  PmBalance.dAngleb = dModul_PrevAngle+ (90.0/double(iDivisionNumber_90deg))*dModul_DirID;
  if(PmBalance.dAngleb > 180.0)
    PmBalance.dAngleb = PmBalance.dAngleb - 360.0;
  if(PmBalance.dAngleb <= -180.0)
    PmBalance.dAngleb = PmBalance.dAngleb + 360.0;
  PmBalance.dAngleb = PmBalance.dAngleb/180.0*CMathFunc::m_pi;
	
  PmBalance.dPb    = pmcvDataPrev.dPb + (pmcvDataNext.dPb - pmcvDataPrev.dPb)*dModul_DirID;
  PmBalance.dMby   = pmcvDataPrev.dMby*dModul_PrevMy + (pmcvDataNext.dMby*dModul_NextMy - pmcvDataPrev.dMby*dModul_PrevMy)*dModul_DirID;
  PmBalance.dMbz   = pmcvDataPrev.dMbz*dModul_PrevMz + (pmcvDataNext.dMbz*dModul_NextMz - pmcvDataPrev.dMbz*dModul_PrevMz)*dModul_DirID;
	PmBalance.desib  = pmcvDataPrev.desib*dModul_PrevMz + (pmcvDataNext.desib*dModul_NextMz - pmcvDataPrev.desib*dModul_PrevMz)*dModul_DirID;
  
  return TRUE;
}

BOOL CMdgnPmTool::Get_PiercePointForContactPlane(_UMD_3DPM_DATA& Umd3DPmData, double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, int iFitDirID, int iFitUnitID, int& iLeftFitDirID, double& dModul_DirID)
{ 
	int iDivisionNumber_90deg = Umd3DPmData.iDivisionNumber_90deg;
  if(Umd3DPmData.PmcvData.GetSize() < iDivisionNumber_90deg + 1)return FALSE;

  dPn = dMn_y = dMn_z = desi = dXn = 0.0;
  CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&> arPmcvUnit;
  arPmcvUnit.SetSize(3);

  _UMD_PMCV_DATA_POLY pmcvData_Center,    pmcvData_Next,    pmcvData_Before;
  int                 iPmcvDataID_Center, iPmcvDataID_Next, iPmcvDataID_Before;
  double              dModul_My_Center,   dModul_My_Next,   dModul_My_Before;
  double              dModul_Mz_Center,   dModul_Mz_Next,   dModul_Mz_Before;  
  double              dLine_i[3], dLine_j[3], dPoint[3];
  double              dLine[2][2], dStartP[2], dEndP[2], dCross[2];
  
  Get_ResultPmcvDataID(Umd3DPmData, iFitDirID, iPmcvDataID_Center, dModul_My_Center, dModul_Mz_Center);
  pmcvData_Center = Umd3DPmData.PmcvData.GetAt(iPmcvDataID_Center);
  
  int iPmcvDataCount = pmcvData_Center.arPmcvUnit.GetCount();
  if(iFitUnitID > 1 && iFitUnitID < iPmcvDataCount)
  {    
    //before  center  next
    //    ┌──┬──┐iFitUnitID-1
    //상면│①/ │③/ │①③의 경우 iFitUnitID == 2 일때는 검토 안함
    //    │ /②│ /④│
    //    ├──┼──┤iFitUnitID
    //하면│⑤/ │⑦/ │⑥⑧의 경우 iFitUnitID == iPmcvDataCount-1 일때는 검토 안함
    //    │ /⑥│ /⑧│
    //    └──┴──┘iFitUnitID+1
    Get_ResultPmcvDataID(Umd3DPmData, iFitDirID+1, iPmcvDataID_Next, dModul_My_Next, dModul_Mz_Next);
    pmcvData_Next = Umd3DPmData.PmcvData.GetAt(iPmcvDataID_Next);
    Get_ResultPmcvDataID(Umd3DPmData, iFitDirID-1, iPmcvDataID_Before, dModul_My_Before, dModul_Mz_Before);
    pmcvData_Before = Umd3DPmData.PmcvData.GetAt(iPmcvDataID_Before);

    //상면 검토
    //①평면 검토
    arPmcvUnit[0] = Get_CngSignPmcvUnit(pmcvData_Before.arPmcvUnit[iFitUnitID-1], dModul_My_Before, dModul_Mz_Before);
    arPmcvUnit[1] = Get_CngSignPmcvUnit(pmcvData_Before.arPmcvUnit[iFitUnitID],   dModul_My_Before, dModul_Mz_Before);
    arPmcvUnit[2] = Get_CngSignPmcvUnit(pmcvData_Center.arPmcvUnit[iFitUnitID-1], dModul_My_Center, dModul_Mz_Center);
    if(iFitUnitID != 2 && Cal_PierceCheckToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, arPmcvUnit)) // 꼭지점과 접할경우 Pass
    {
      if(Cal_PiercePointToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit)) 
      { // 좌측에 있는 DirID와 해당 ID에서 떨어진 비율을 계산       
        dLine_i[0] = arPmcvUnit[0].dPn;  dLine_i[1] = arPmcvUnit[0].dMny;  dLine_i[2] = arPmcvUnit[0].dMnz; 
        dLine_j[0] = arPmcvUnit[1].dPn;  dLine_j[1] = arPmcvUnit[1].dMny;  dLine_j[2] = arPmcvUnit[1].dMnz;
        dPoint[0]  = arPmcvUnit[2].dPn;  dPoint[1]  = arPmcvUnit[2].dMny;  dPoint[2]  = arPmcvUnit[2].dMnz;
        double dLength = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint);
        dPoint[0]  = dPn;                dPoint[1]  = dMn_y;               dPoint[2]  = dMn_z;
        dModul_DirID = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint)/dLength;
        iLeftFitDirID= iFitDirID-1;

        arPmcvUnit.RemoveAll();
        return TRUE;
      }
    }
    //②평면 검토
    arPmcvUnit[0] = Get_CngSignPmcvUnit(pmcvData_Center.arPmcvUnit[iFitUnitID],   dModul_My_Center, dModul_Mz_Center);
    if(Cal_PierceCheckToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, arPmcvUnit)) 
    {
      if(Cal_PiercePointToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit)) 
      { // 좌측에 있는 DirID와 해당 ID에서 떨어진 비율을 계산  
        if(iFitUnitID != 2)
        {
          dLine_i[0] = arPmcvUnit[0].dPn;  dLine_i[1] = arPmcvUnit[0].dMny;  dLine_i[2] = arPmcvUnit[0].dMnz; 
          dLine_j[0] = arPmcvUnit[2].dPn;  dLine_j[1] = arPmcvUnit[2].dMny;  dLine_j[2] = arPmcvUnit[2].dMnz;
          dPoint[0]  = arPmcvUnit[1].dPn;  dPoint[1]  = arPmcvUnit[1].dMny;  dPoint[2]  = arPmcvUnit[1].dMnz;
          double dLength = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint);
          dPoint[0]  = dPn;                dPoint[1]  = dMn_y;               dPoint[2]  = dMn_z;
          dModul_DirID = 1.0 - CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint)/dLength;
          iLeftFitDirID= iFitDirID-1;
        }
        else
        {
          dLine[0][0] = arPmcvUnit[2].dMny;  dLine[0][1] = arPmcvUnit[2].dMnz;
          dLine[1][0] = dMn_y;  dLine[1][1] = dMn_z;
          dStartP[0] = arPmcvUnit[1].dMny;  dStartP[1] = arPmcvUnit[1].dMnz;
          dEndP[0]   = arPmcvUnit[0].dMny;  dEndP[1]   = arPmcvUnit[0].dMnz;
          if(CMathFunc::mathLineSegCross2D(dLine, dStartP, dEndP, dCross) == 0)
          {
            arPmcvUnit.RemoveAll();
            return FALSE;
          }
          double dLength = CMathFunc::mathLength(dStartP[0], dStartP[1], dEndP[0], dEndP[1]);
          dModul_DirID = CMathFunc::mathLength(arPmcvUnit[1].dMny, arPmcvUnit[1].dMnz, dCross[0], dCross[1])/dLength;
          iLeftFitDirID = iFitDirID-1;
        }

        arPmcvUnit.RemoveAll();
        return TRUE;
      }
    }
    //③평면 검토
    arPmcvUnit[1] = Get_CngSignPmcvUnit(pmcvData_Next.arPmcvUnit[iFitUnitID-1],   dModul_My_Next, dModul_Mz_Next);
    if(iFitUnitID != 2 && Cal_PierceCheckToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, arPmcvUnit)) // 꼭지점과 접할경우 Pass
    {
      if(Cal_PiercePointToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit)) 
      { // 좌측에 있는 DirID와 해당 ID에서 떨어진 비율을 계산       
        dLine_i[0] = arPmcvUnit[0].dPn;  dLine_i[1] = arPmcvUnit[0].dMny;  dLine_i[2] = arPmcvUnit[0].dMnz; 
        dLine_j[0] = arPmcvUnit[2].dPn;  dLine_j[1] = arPmcvUnit[2].dMny;  dLine_j[2] = arPmcvUnit[2].dMnz;
        dPoint[0]  = arPmcvUnit[1].dPn;  dPoint[1]  = arPmcvUnit[1].dMny;  dPoint[2]  = arPmcvUnit[1].dMnz;
        double dLength = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint);
        dPoint[0]  = dPn;                dPoint[1]  = dMn_y;               dPoint[2]  = dMn_z;
        dModul_DirID = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint)/dLength;
        iLeftFitDirID= iFitDirID;

        arPmcvUnit.RemoveAll();
        return TRUE;
      }
    }
    //④평면 검토
    arPmcvUnit[2] = Get_CngSignPmcvUnit(pmcvData_Next.arPmcvUnit[iFitUnitID],     dModul_My_Next, dModul_Mz_Next);
    if(Cal_PierceCheckToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, arPmcvUnit)) 
    {
      if(Cal_PiercePointToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit)) 
      { // 좌측에 있는 DirID와 해당 ID에서 떨어진 비율을 계산       
        if(iFitUnitID != 2)
        {
          dLine_i[0] = arPmcvUnit[1].dPn;  dLine_i[1] = arPmcvUnit[1].dMny;  dLine_i[2] = arPmcvUnit[1].dMnz; 
          dLine_j[0] = arPmcvUnit[2].dPn;  dLine_j[1] = arPmcvUnit[2].dMny;  dLine_j[2] = arPmcvUnit[2].dMnz;
          dPoint[0]  = arPmcvUnit[0].dPn;  dPoint[1]  = arPmcvUnit[0].dMny;  dPoint[2]  = arPmcvUnit[0].dMnz;
          double dLength = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint);
          dPoint[0]  = dPn;                dPoint[1]  = dMn_y;               dPoint[2]  = dMn_z;
          dModul_DirID = 1.0 - CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint)/dLength;
          iLeftFitDirID= iFitDirID;
        }
        else
        {
          dLine[0][0] = arPmcvUnit[1].dMny;  dLine[0][1] = arPmcvUnit[1].dMnz;
          dLine[1][0] = dMn_y;  dLine[1][1] = dMn_z;
          dStartP[0] = arPmcvUnit[0].dMny;  dStartP[1] = arPmcvUnit[0].dMnz;
          dEndP[0]   = arPmcvUnit[2].dMny;  dEndP[1]   = arPmcvUnit[2].dMnz;
          if(CMathFunc::mathLineSegCross2D(dLine, dStartP, dEndP, dCross) == 0)
          {
            arPmcvUnit.RemoveAll();
            return FALSE;
          }
          double dLength = CMathFunc::mathLength(dStartP[0], dStartP[1], dEndP[0], dEndP[1]);
          dModul_DirID = CMathFunc::mathLength(arPmcvUnit[0].dMny, arPmcvUnit[0].dMnz, dCross[0], dCross[1])/dLength;
          iLeftFitDirID = iFitDirID;
        }

        arPmcvUnit.RemoveAll();
        return TRUE;
      }
    }
    //하면검토
    //⑤평면 검토
    arPmcvUnit[0] = Get_CngSignPmcvUnit(pmcvData_Before.arPmcvUnit[iFitUnitID],   dModul_My_Before, dModul_Mz_Before);
    arPmcvUnit[1] = Get_CngSignPmcvUnit(pmcvData_Before.arPmcvUnit[iFitUnitID+1], dModul_My_Before, dModul_Mz_Before);
    arPmcvUnit[2] = Get_CngSignPmcvUnit(pmcvData_Center.arPmcvUnit[iFitUnitID],   dModul_My_Center, dModul_Mz_Center);
    if(Cal_PierceCheckToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, arPmcvUnit)) 
    {
      if(Cal_PiercePointToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit)) 
      { // 좌측에 있는 DirID와 해당 ID에서 떨어진 비율을 계산       
        if(iFitUnitID != iPmcvDataCount-1)
        {
          dLine_i[0] = arPmcvUnit[0].dPn;  dLine_i[1] = arPmcvUnit[0].dMny;  dLine_i[2] = arPmcvUnit[0].dMnz; 
          dLine_j[0] = arPmcvUnit[1].dPn;  dLine_j[1] = arPmcvUnit[1].dMny;  dLine_j[2] = arPmcvUnit[1].dMnz;
          dPoint[0]  = arPmcvUnit[2].dPn;  dPoint[1]  = arPmcvUnit[2].dMny;  dPoint[2]  = arPmcvUnit[2].dMnz;
          double dLength = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint);
          dPoint[0]  = dPn;                dPoint[1]  = dMn_y;               dPoint[2]  = dMn_z;
          dModul_DirID = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint)/dLength;
          iLeftFitDirID= iFitDirID-1;
        }
        else
        {
          dLine[0][0] = arPmcvUnit[1].dMny;  dLine[0][1] = arPmcvUnit[1].dMnz;
          dLine[1][0] = dMn_y;  dLine[1][1] = dMn_z;
          dStartP[0] = arPmcvUnit[0].dMny;  dStartP[1] = arPmcvUnit[0].dMnz;
          dEndP[0]   = arPmcvUnit[2].dMny;  dEndP[1]   = arPmcvUnit[2].dMnz;
          if(CMathFunc::mathLineSegCross2D(dLine, dStartP, dEndP, dCross) == 0)
          {
            arPmcvUnit.RemoveAll();
            return FALSE;
          }
          double dLength = CMathFunc::mathLength(dStartP[0], dStartP[1], dEndP[0], dEndP[1]);
          dModul_DirID = CMathFunc::mathLength(arPmcvUnit[0].dMny, arPmcvUnit[0].dMnz, dCross[0], dCross[1])/dLength;
          iLeftFitDirID = iFitDirID-1;
        }

        arPmcvUnit.RemoveAll();
        return TRUE;
      }
    }
    //⑥평면 검토
    arPmcvUnit[0] = Get_CngSignPmcvUnit(pmcvData_Center.arPmcvUnit[iFitUnitID+1], dModul_My_Center, dModul_Mz_Center);
    if(iFitUnitID != iPmcvDataCount-1 && Cal_PierceCheckToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, arPmcvUnit)) // 꼭지점과 접할경우 Pass
    {
      if(Cal_PiercePointToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit)) 
      { // 좌측에 있는 DirID와 해당 ID에서 떨어진 비율을 계산       
        dLine_i[0] = arPmcvUnit[0].dPn;  dLine_i[1] = arPmcvUnit[0].dMny;  dLine_i[2] = arPmcvUnit[0].dMnz; 
        dLine_j[0] = arPmcvUnit[2].dPn;  dLine_j[1] = arPmcvUnit[2].dMny;  dLine_j[2] = arPmcvUnit[2].dMnz;
        dPoint[0]  = arPmcvUnit[1].dPn;  dPoint[1]  = arPmcvUnit[1].dMny;  dPoint[2]  = arPmcvUnit[1].dMnz;
        double dLength = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint);
        dPoint[0]  = dPn;                dPoint[1]  = dMn_y;               dPoint[2]  = dMn_z;
        dModul_DirID = 1.0 - CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint)/dLength;
        iLeftFitDirID= iFitDirID-1;

        arPmcvUnit.RemoveAll();
        return TRUE;
      }
    }
    //⑦평면 검토
    arPmcvUnit[1] = Get_CngSignPmcvUnit(pmcvData_Next.arPmcvUnit[iFitUnitID],     dModul_My_Next, dModul_Mz_Next);
    if(Cal_PierceCheckToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, arPmcvUnit))
    {
      if(Cal_PiercePointToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit)) 
      { // 좌측에 있는 DirID와 해당 ID에서 떨어진 비율을 계산  
        if(iFitUnitID != iPmcvDataCount-1)
        {
          dLine_i[0] = arPmcvUnit[0].dPn;  dLine_i[1] = arPmcvUnit[0].dMny;  dLine_i[2] = arPmcvUnit[0].dMnz; 
          dLine_j[0] = arPmcvUnit[2].dPn;  dLine_j[1] = arPmcvUnit[2].dMny;  dLine_j[2] = arPmcvUnit[2].dMnz;
          dPoint[0]  = arPmcvUnit[1].dPn;  dPoint[1]  = arPmcvUnit[1].dMny;  dPoint[2]  = arPmcvUnit[1].dMnz;
          double dLength = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint);
          dPoint[0]  = dPn;                dPoint[1]  = dMn_y;               dPoint[2]  = dMn_z;
          dModul_DirID = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint)/dLength;
          iLeftFitDirID= iFitDirID;
        }
        else
        {
          dLine[0][0] = arPmcvUnit[0].dMny;  dLine[0][1] = arPmcvUnit[0].dMnz;
          dLine[1][0] = dMn_y;  dLine[1][1] = dMn_z;
          dStartP[0] = arPmcvUnit[2].dMny;  dStartP[1] = arPmcvUnit[2].dMnz;
          dEndP[0]   = arPmcvUnit[1].dMny;  dEndP[1]   = arPmcvUnit[1].dMnz;
          if(CMathFunc::mathLineSegCross2D(dLine, dStartP, dEndP, dCross) == 0)
          {
            arPmcvUnit.RemoveAll();
            return FALSE;
          }
          double dLength = CMathFunc::mathLength(dStartP[0], dStartP[1], dEndP[0], dEndP[1]);
          dModul_DirID = CMathFunc::mathLength(arPmcvUnit[2].dMny, arPmcvUnit[2].dMnz, dCross[0], dCross[1])/dLength;
          iLeftFitDirID = iFitDirID;
        }

        arPmcvUnit.RemoveAll();
        return TRUE;
      }
    }
    //⑧평면 검토
    arPmcvUnit[2] = Get_CngSignPmcvUnit(pmcvData_Next.arPmcvUnit[iFitUnitID+1],   dModul_My_Next, dModul_Mz_Next);
    if(iFitUnitID != iPmcvDataCount-1 && Cal_PierceCheckToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, arPmcvUnit)) // 꼭지점과 접할경우 Pass
    {
      if(Cal_PiercePointToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit)) 
      { // 좌측에 있는 DirID와 해당 ID에서 떨어진 비율을 계산       
        dLine_i[0] = arPmcvUnit[1].dPn;  dLine_i[1] = arPmcvUnit[1].dMny;  dLine_i[2] = arPmcvUnit[1].dMnz; 
        dLine_j[0] = arPmcvUnit[2].dPn;  dLine_j[1] = arPmcvUnit[2].dMny;  dLine_j[2] = arPmcvUnit[2].dMnz;
        dPoint[0]  = arPmcvUnit[0].dPn;  dPoint[1]  = arPmcvUnit[0].dMny;  dPoint[2]  = arPmcvUnit[0].dMnz;
        double dLength = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint);
        dPoint[0]  = dPn;                dPoint[1]  = dMn_y;               dPoint[2]  = dMn_z;
        dModul_DirID = 1.0 - CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint)/dLength;
        iLeftFitDirID= iFitDirID;

        arPmcvUnit.RemoveAll();
        return TRUE;
      }
    }
  }
  else if(iFitUnitID == 1 || iFitUnitID == iPmcvDataCount)
  {// 가장 가까운 방향이 꼭지점일 경우
    int iNexUnitID;
    if(iFitUnitID == 1)
      iNexUnitID = 2;
    else
      iNexUnitID = iFitUnitID-1;    

    arPmcvUnit[2] = Get_CngSignPmcvUnit(pmcvData_Center.arPmcvUnit[iFitUnitID], dModul_My_Center, dModul_Mz_Center);
    arPmcvUnit[0] = Get_CngSignPmcvUnit(pmcvData_Center.arPmcvUnit[iNexUnitID], dModul_My_Center, dModul_Mz_Center);
   
    //    
    double vectorD[3], vectorN[3];
    double dPmax = pmcvData_Center.arPmcvUnit[1].dPn;
    double dMmax = CMathFunc::mathLength(pmcvData_Center.dMby, pmcvData_Center.dMbz);
    vectorD[0] = (dPd-dPs)/dPmax;                vectorD[1] = (dMd_y-dMs_y)/dMmax;               vectorD[2] = (dMd_z-dMs_z)/dMmax;
    vectorN[0] = (arPmcvUnit[2].dPn-dPs)/dPmax;  vectorN[1] = (arPmcvUnit[2].dMny-dMs_y)/dMmax;  vectorN[2] = (arPmcvUnit[2].dMnz-dMs_z)/dMmax;
    CMathFunc::mathNormalize(vectorD, vectorD);
    CMathFunc::mathNormalize(vectorN, vectorN);
    if(fabs(CMathFunc::mathCrossAngle(vectorD, vectorN)) < m_dTolAngle)
    {
      dLine[0][0] = dMs_y;  dLine[0][1] = dMs_z;
      dLine[1][0] = dMd_y;  dLine[1][1] = dMd_z;
      if(CMathFunc::mathLength(dMd_y-dMs_y, dMd_z-dMs_z) < dMmax*m_dTolM0)
      {// 0일경우에는 Y축방향으로 설정함
        dLine[1][0] = dMd_y+0.0; dLine[1][1] = dMd_z+1.0;
      }

      for(int i=1 ; i<=4*iDivisionNumber_90deg ; i++)
      {
        Get_ResultPmcvDataID(Umd3DPmData, iFitDirID+i, iPmcvDataID_Next, dModul_My_Next, dModul_Mz_Next);
        pmcvData_Next = Umd3DPmData.PmcvData.GetAt(iPmcvDataID_Next);
        arPmcvUnit[i%2] = Get_CngSignPmcvUnit(pmcvData_Next.arPmcvUnit[iNexUnitID], dModul_My_Next, dModul_Mz_Next);
        
        dStartP[0] = arPmcvUnit[0].dMny;  dStartP[1] = arPmcvUnit[0].dMnz;
        dEndP[0]   = arPmcvUnit[1].dMny;  dEndP[1]   = arPmcvUnit[1].dMnz;

        if(CMathFunc::mathLineSegCross2D(dLine, dStartP, dEndP, dCross) == 1)
        {
          double dLength = CMathFunc::mathLength(dStartP[0], dStartP[1], dEndP[0], dEndP[1]);
          dModul_DirID = CMathFunc::mathLength(arPmcvUnit[(i-1)%2].dMny, arPmcvUnit[(i-1)%2].dMnz, dCross[0], dCross[1])/dLength;
          iLeftFitDirID = iFitDirID+i-1;

					dPn   = arPmcvUnit[2].dPn;
					dMn_y = arPmcvUnit[2].dMny;
					dMn_z = arPmcvUnit[2].dMnz;
          dXn   = arPmcvUnit[2].dXn;
					desi  = arPmcvUnit[2].desi;
          arPmcvUnit.RemoveAll();
          return TRUE;
        }
      }
      arPmcvUnit.RemoveAll();
      return FALSE;
    }          
    
    for(int i=1 ; i<=4*iDivisionNumber_90deg ; i++)
    {
      Get_ResultPmcvDataID(Umd3DPmData, iFitDirID+i, iPmcvDataID_Next, dModul_My_Next, dModul_Mz_Next);
      pmcvData_Next = Umd3DPmData.PmcvData.GetAt(iPmcvDataID_Next);
      arPmcvUnit[i%2] = Get_CngSignPmcvUnit(pmcvData_Next.arPmcvUnit[iNexUnitID], dModul_My_Next, dModul_Mz_Next);
      if(Cal_PierceCheckToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, arPmcvUnit)) 
      {
        if(Cal_PiercePointToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit)) 
        { //
          dLine[0][0] = arPmcvUnit[2].dMny;  dLine[0][1] = arPmcvUnit[2].dMnz;
          dLine[1][0] = dMn_y;  dLine[1][1] = dMn_z;
          dStartP[0] = arPmcvUnit[0].dMny;  dStartP[1] = arPmcvUnit[0].dMnz;
          dEndP[0]   = arPmcvUnit[1].dMny;  dEndP[1]   = arPmcvUnit[1].dMnz;

          if(CMathFunc::mathLineSegCross2D(dLine, dStartP, dEndP, dCross) == 0)
          {// 교차점을 찾을 수 없을 경우 dStartP, dEndP 두점이 근사하게 교차하는지 여부를 확인함
            double uLineX, uLineY, uStartPX, uStartPY, uEndPX, uEndPY;
            double dStartPCrossAngle, dEndPCrossAngle;

            CMathFunc::mathNormalize(dLine[1][0]-dLine[0][0], dLine[1][1]-dLine[0][1], uLineX,   uLineY);
            CMathFunc::mathNormalize(dStartP[0]-dLine[0][0],  dStartP[1]-dLine[0][1],  uStartPX, uStartPY);
            CMathFunc::mathNormalize(dEndP[0]-dLine[0][0],    dEndP[1]-dLine[0][1],    uEndPX,   uEndPY);

            dStartPCrossAngle = CMathFunc::mathCrossAngle2D(uLineX, uLineY, uStartPX, uStartPY);
            dEndPCrossAngle   = CMathFunc::mathCrossAngle2D(uLineX, uLineY, uEndPX,   uEndPY);

            if(fabs(dStartPCrossAngle) < m_dTolAngle/10.0)
            { dCross[0] = dStartP[0];  dCross[1] = dStartP[1]; }
            else if(fabs(dEndPCrossAngle) < m_dTolAngle/10.0)
            { dCross[0] = dEndP[0];  dCross[1] = dEndP[1]; }
            else
            {// 교차접이 없을 경우 종료
              arPmcvUnit.RemoveAll();
              return FALSE;
            }
          }
          double dLength = CMathFunc::mathLength(dStartP[0], dStartP[1], dEndP[0], dEndP[1]);
          dModul_DirID = CMathFunc::mathLength(arPmcvUnit[(i-1)%2].dMny, arPmcvUnit[(i-1)%2].dMnz, dCross[0], dCross[1])/dLength;
          iLeftFitDirID = iFitDirID+i-1;

          arPmcvUnit.RemoveAll();
          return TRUE;
        }
      }
    }
  }
  arPmcvUnit.RemoveAll();
  return FALSE;

}


BOOL CMdgnPmTool::Get_FitPmcvUnitID(_UMD_3DPM_DATA& Umd3DPmData, double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, int& iFitDirID, int& iFitUnitID)
{
	int iDivisionNumber_90deg = Umd3DPmData.iDivisionNumber_90deg;
  iFitDirID = iFitUnitID = 0;
  if(Umd3DPmData.PmcvData.GetSize() < iDivisionNumber_90deg + 1)return FALSE;

  double (*momentB_MzMy)[2];//(x:Mz  y:My)
  double vectorMdzMdyPd[3];// 외력에 대한 3D좌표에 대한 단위벡터
  double vectorMdPd[2];    // 외력에 대한 모멘트와 축하중에 대한 2D좌표에 대한 단위벡터 
  double vectorMdzMdy[2];  // 외력에 대한 z,y모멘트에 대한 2D좌표에 대한 단위벡터
  double vectorMnzMnyPn[3];// PM-Curve상의3D좌표에 대한 단위벡터
  double vectorMnPn[2];    // PM-Curve상의3D좌표에 대한 모멘트와 축하중에 대한 2D좌표에 대한 단위벡터 
  double vectorMnzMny[2];  // PM-Curve상의3D좌표에 대한 z,y모멘트에 대한 2D좌표에 대한 단위벡터

  double dPmax, dMmax, dMn;// Scale조절을 위한 최대 축하중 및 최대 모멘트
  dPmax = dMmax = 0.0;

  Get_UnitVectorToMzMyP(dPd-dPs, dMd_y-dMs_y, dMd_z-dMs_z, vectorMdzMdyPd, vectorMdPd, vectorMdzMdy);//여기서는 vectorMdPd와 vectorMdzMdy를 쓰기 위해서임
                                                                                     //이후 각도의 차이를 찾기 위한 vectorMdzMdyPd는 Pmax와 Mmax로 초기화 하여야함
  int iNumber = 4*iDivisionNumber_90deg;
  int iPmcvDataID;
  double dModul_My, dModul_Mz;
  
  // 최대 모멘트의 발생방향에 대한 단위 Vector배열및 외력의 방향과 가장 가까운 방향에 대한 0~분할수-1까지의 ID를 찾음
  momentB_MzMy = new double[iNumber][2];
  double dAngleMin = 180;
  double dAngle;
  for(int i=0 ; i<iNumber ; i++)
  {
    Get_ResultPmcvDataID(Umd3DPmData, i, iPmcvDataID, dModul_My, dModul_Mz);
    //pmcvData = Umd3DPmData.PmcvData.GetAt(iPmcvDataID);
    dMn = CMathFunc::mathLength(Umd3DPmData.PmcvData[iPmcvDataID].dMbz, Umd3DPmData.PmcvData[iPmcvDataID].dMby);
    if(dMmax < dMn) dMmax = dMn;//최대 모멘트를 찾음(Scale조절을 위해)
    CMathFunc::mathNormalize((Umd3DPmData.PmcvData[iPmcvDataID].dMbz-dMs_z)*dModul_Mz, (Umd3DPmData.PmcvData[iPmcvDataID].dMby-dMs_y)*dModul_My, momentB_MzMy[i][0], momentB_MzMy[i][1]);

    dAngle = fabs(CMathFunc::mathCrossAngle2D(vectorMdzMdy[0], vectorMdzMdy[1], momentB_MzMy[i][0], momentB_MzMy[i][1]));
    if(dAngle < dAngleMin)
    {
      iFitDirID = i;  dAngleMin = dAngle;
    }
  }
  Get_ResultPmcvDataID(Umd3DPmData, iFitDirID, iPmcvDataID, dModul_My, dModul_Mz);
  int iPmcvSize = Umd3DPmData.PmcvData[iPmcvDataID].arPmcvUnit.GetCount();
  _UMD_PMCV_UNIT_POLY pmcvUnit;
  
  // 순수 압축 및 순수 인장시에 모멘트가 0이 아닌경우에는 해당외력의 방향이 PM-Curve상에 없을 수가 있으므로, 이때에는 반대방향의 PM-Curve방향을 잡는다.  
  double dAngleTop, dAngleBottom;
  
  pmcvUnit = Umd3DPmData.PmcvData[iPmcvDataID].arPmcvUnit[1];
  dPmax = pmcvUnit.dPn; // PM-Curve의 첫번째 점은 순수 압축이므로 최대 축하중임
 
	double nMzMyP_Org[3];
	nMzMyP_Org[0] = dMs_z/dMmax;
	nMzMyP_Org[1] = dMs_y/dMmax;
	nMzMyP_Org[2] = dPs/dPmax;

  Get_UnitVectorToMzMyP((dPd-dPs)/dPmax, (dMd_y-dMs_y)/dMmax, (dMd_z-dMs_z)/dMmax, vectorMdzMdyPd, vectorMdPd, vectorMdzMdy);// Mmax와 Pmax로 일률화 시킨후 단위 Vector를 만듬
  
  Get_UnitVectorToMzMyP((pmcvUnit.dPn-dPs)/dPmax, (pmcvUnit.dMny-dMs_y)/dMmax, (pmcvUnit.dMnz-dMs_z)/dMmax, vectorMnzMnyPn, vectorMnPn, vectorMnzMny);
  dAngleTop = CMathFunc::mathCrossAngle2DSign(vectorMdPd[0], vectorMdPd[1], vectorMnPn[0], vectorMnPn[1]);
  pmcvUnit = Umd3DPmData.PmcvData[iPmcvDataID].arPmcvUnit[iPmcvSize];
  Get_UnitVectorToMzMyP((pmcvUnit.dPn-dPs)/dPmax, (pmcvUnit.dMny-dMs_y)/dMmax, (pmcvUnit.dMnz-dMs_z)/dMmax, vectorMnzMnyPn, vectorMnPn, vectorMnzMny);
  dAngleBottom = CMathFunc::mathCrossAngle2DSign(vectorMdPd[0], vectorMdPd[1], vectorMnPn[0], vectorMnPn[1]);
  if(dAngleTop < 0 || dAngleBottom > 0)
  {//해당외력의 방향이 PM-Curve상에 없을때
    for(int i=0 ; i<iNumber ; i++)
    {
      dAngle = fabs(CMathFunc::mathCrossAngle2D(-vectorMdzMdy[0], -vectorMdzMdy[1], momentB_MzMy[i][0], momentB_MzMy[i][1]));
      if(dAngle < dAngleMin)
      {
        iFitDirID = i;  dAngleMin = dAngle;
      }
    }    
  }

  //전체 3D 좌표값과 외력방향과 이루는 각이 가장 작은 점을 찾음
  int iStartID  = iFitDirID;
  int iNextID   = iStartID+1;
  int iBeforeID = iStartID-1;
  int iUnitStartID, iUnitNextID, iUnitBeforeID;
  double dStartAngleMin, dNextAngleMin, dBeforeAngleMin;
  BOOL bFindCheck = FALSE;

  Get_AngleMinIDToPmCurve_NonScale(Umd3DPmData, nMzMyP_Org, vectorMdzMdyPd, iStartID,  dPmax, dMmax, iUnitStartID,  dStartAngleMin); 
  Get_AngleMinIDToPmCurve_NonScale(Umd3DPmData, nMzMyP_Org, vectorMdzMdyPd, iNextID,   dPmax, dMmax, iUnitNextID,   dNextAngleMin); 
  Get_AngleMinIDToPmCurve_NonScale(Umd3DPmData, nMzMyP_Org, vectorMdzMdyPd, iBeforeID, dPmax, dMmax, iUnitBeforeID, dBeforeAngleMin); 
  for(int i=0 ; i<4*iDivisionNumber_90deg ; i++)
  {     
    if(dStartAngleMin <= dNextAngleMin && dStartAngleMin <= dBeforeAngleMin)
    {// dStartAngleMin 가장 작을때
      iFitDirID = iStartID;  iFitUnitID = iUnitStartID;
      bFindCheck = TRUE;  
      break;
    }
    //else if(dNextAngleMin < dStartAngleMin && dNextAngleMin < dBeforeAngleMin)
    else if(dNextAngleMin < dStartAngleMin && dNextAngleMin <= dBeforeAngleMin)// SHIN(06.06.29)좌우가 모두 작을 경우 한쪽이 큰것으로 취급함
    {// dNextAngleMin 가장 작을때
      iBeforeID       = iStartID;
      iUnitBeforeID   = iUnitStartID;
      dBeforeAngleMin = dStartAngleMin;
      iStartID        = iNextID;
      iUnitStartID    = iUnitNextID;
      dStartAngleMin  = dNextAngleMin;
      Get_AngleMinIDToPmCurve_NonScale(Umd3DPmData, nMzMyP_Org, vectorMdzMdyPd, iStartID+1, dPmax, dMmax, iUnitNextID,   dNextAngleMin); 
    }
    else if(dBeforeAngleMin < dStartAngleMin && dBeforeAngleMin < dNextAngleMin)
    {// dNextAngleMin 가장 작을때
      iNextID         = iStartID;
      iUnitNextID     = iUnitStartID;
      dNextAngleMin   = dStartAngleMin;
      iStartID        = iBeforeID;
      iUnitStartID    = iUnitBeforeID;
      dStartAngleMin  = dBeforeAngleMin;
      Get_AngleMinIDToPmCurve_NonScale(Umd3DPmData, nMzMyP_Org, vectorMdzMdyPd, iStartID-1, dPmax, dMmax, iUnitBeforeID,   dBeforeAngleMin); 
    }
    else if((iUnitStartID == 1 && iUnitNextID == 1 && iUnitBeforeID ==1) || (iUnitStartID==iPmcvSize && iUnitNextID==iPmcvSize && iUnitBeforeID==iPmcvSize))
    {// 모두 꼭지점에 가까울때
      iFitDirID = iStartID;  iFitUnitID = iUnitStartID;
      bFindCheck = TRUE;  
      break;
    }   
    else
      ASSERT(0);
  }

  delete[] momentB_MzMy;
  return TRUE;
}

BOOL CMdgnPmTool::Get_PiercePontForDirection(_UMD_3DPM_DATA& Umd3DPmData, double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, int iDirID_1, int iDirID_2, double& dModul_DirID)
{
	int iDivisionNumber_90deg = Umd3DPmData.iDivisionNumber_90deg;
  if(Umd3DPmData.PmcvData.GetSize() < iDivisionNumber_90deg + 1)return FALSE;

  dPn = dMn_y = dMn_z = 0.0;
  CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&> arPmcvUnit;
  arPmcvUnit.SetSize(3);

  _UMD_PMCV_DATA_POLY pmcvData_1,    pmcvData_2;
  int                 iPmcvDataID_1, iPmcvDataID_2;
  double              dModul_My_1,   dModul_My_2;
  double              dModul_Mz_1,   dModul_Mz_2;
  double              dLine_i[3], dLine_j[3], dPoint[3];
  double              dLine[2][2], dStartP[2], dEndP[2], dCross[2];
  
  Get_ResultPmcvDataID(Umd3DPmData, iDirID_1, iPmcvDataID_1, dModul_My_1, dModul_Mz_1);
  Get_ResultPmcvDataID(Umd3DPmData, iDirID_2, iPmcvDataID_2, dModul_My_2, dModul_Mz_2);
  pmcvData_1 = Umd3DPmData.PmcvData.GetAt(iPmcvDataID_1);  
  pmcvData_2 = Umd3DPmData.PmcvData.GetAt(iPmcvDataID_2);
  
  arPmcvUnit.RemoveAll();
  arPmcvUnit.SetSize(3);

  int iPmcvDataCount = pmcvData_1.arPmcvUnit.GetCount();

  //윗부분 꼭지점 3각형 계산
  arPmcvUnit[0] = Get_CngSignPmcvUnit(pmcvData_1.arPmcvUnit[1], dModul_My_1, dModul_Mz_1);
  arPmcvUnit[1] = Get_CngSignPmcvUnit(pmcvData_1.arPmcvUnit[2], dModul_My_1, dModul_Mz_1);
  arPmcvUnit[2] = Get_CngSignPmcvUnit(pmcvData_2.arPmcvUnit[2], dModul_My_2, dModul_Mz_2);
  if(Cal_PierceCheckToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, arPmcvUnit)) 
  {
    if(Cal_PiercePointToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit)) 
    {
      dLine[0][0] = arPmcvUnit[0].dMny;  dLine[0][1] = arPmcvUnit[0].dMnz;
      dLine[1][0] = dMn_y;  dLine[1][1] = dMn_z;

      dStartP[0] = arPmcvUnit[1].dMny;  dStartP[1] = arPmcvUnit[1].dMnz;
      dEndP[0]   = arPmcvUnit[2].dMny;  dEndP[1]   = arPmcvUnit[2].dMnz;

      if(CMathFunc::mathLineSegCross2D(dLine, dStartP, dEndP, dCross) == 0)
      {
        arPmcvUnit.RemoveAll();
        return FALSE;
      }
      double dLength = CMathFunc::mathLength(dStartP[0], dStartP[1], dEndP[0], dEndP[1]);
      dModul_DirID = CMathFunc::mathLength(arPmcvUnit[1].dMny, arPmcvUnit[1].dMnz, dCross[0], dCross[1])/dLength;     

      arPmcvUnit.RemoveAll();
      return TRUE;
    }
  }
  //중간 부분 검토
  for(int i=2 ; i<iPmcvDataCount-1 ; i++)
  {
    //①평면 검토
    arPmcvUnit[0] = Get_CngSignPmcvUnit(pmcvData_1.arPmcvUnit[i],  dModul_My_1, dModul_Mz_1);
    arPmcvUnit[1] = Get_CngSignPmcvUnit(pmcvData_1.arPmcvUnit[i+1],dModul_My_1, dModul_Mz_1);
    arPmcvUnit[2] = Get_CngSignPmcvUnit(pmcvData_2.arPmcvUnit[i],  dModul_My_2, dModul_Mz_2);
    if(Cal_PierceCheckToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, arPmcvUnit)) 
    {

      if(Cal_PiercePointToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit)) 
      {
        // 좌측에 있는 DirID와 해당 ID에서 떨어진 비율을 계산       
        dLine_i[0] = arPmcvUnit[0].dPn;  dLine_i[1] = arPmcvUnit[0].dMny;  dLine_i[2] = arPmcvUnit[0].dMnz; 
        dLine_j[0] = arPmcvUnit[1].dPn;  dLine_j[1] = arPmcvUnit[1].dMny;  dLine_j[2] = arPmcvUnit[1].dMnz;
        dPoint[0]  = arPmcvUnit[2].dPn;  dPoint[1]  = arPmcvUnit[2].dMny;  dPoint[2]  = arPmcvUnit[2].dMnz;
        double dLength = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint);
        dPoint[0]  = dPn;                dPoint[1]  = dMn_y;               dPoint[2]  = dMn_z;
        dModul_DirID = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint)/dLength;
       
        arPmcvUnit.RemoveAll();
        return TRUE;
      }
    }
    //②평면 검토
    arPmcvUnit[0] = Get_CngSignPmcvUnit(pmcvData_2.arPmcvUnit[i+1],dModul_My_2, dModul_Mz_2);    
    if(Cal_PierceCheckToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, arPmcvUnit)) 
    {

      if(Cal_PiercePointToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit)) 
      {
        dLine_i[0] = arPmcvUnit[0].dPn;  dLine_i[1] = arPmcvUnit[0].dMny;  dLine_i[2] = arPmcvUnit[0].dMnz; 
        dLine_j[0] = arPmcvUnit[2].dPn;  dLine_j[1] = arPmcvUnit[2].dMny;  dLine_j[2] = arPmcvUnit[2].dMnz;
        dPoint[0]  = arPmcvUnit[1].dPn;  dPoint[1]  = arPmcvUnit[1].dMny;  dPoint[2]  = arPmcvUnit[1].dMnz;
        double dLength = CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint);
        dPoint[0]  = dPn;                dPoint[1]  = dMn_y;               dPoint[2]  = dMn_z;
        dModul_DirID = 1.0 - CMathFunc::mathDistanceFromIntersectPointToLine(dLine_i, dLine_j, dPoint)/dLength;

        arPmcvUnit.RemoveAll();
        return TRUE;
      }
    }
  }
  //아랫부분 꼭지점 3각형 계산
  arPmcvUnit[0] = Get_CngSignPmcvUnit(pmcvData_1.arPmcvUnit[iPmcvDataCount],   dModul_My_1, dModul_Mz_1);
  arPmcvUnit[1] = Get_CngSignPmcvUnit(pmcvData_1.arPmcvUnit[iPmcvDataCount-1], dModul_My_1, dModul_Mz_1);
  arPmcvUnit[2] = Get_CngSignPmcvUnit(pmcvData_2.arPmcvUnit[iPmcvDataCount-1], dModul_My_2, dModul_Mz_2);
  if(Cal_PierceCheckToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, arPmcvUnit)) 
  {

    if(Cal_PiercePointToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit)) 
    {
      dLine[0][0] = arPmcvUnit[0].dMny;  dLine[0][1] = arPmcvUnit[0].dMnz;
      dLine[1][0] = dMn_y;  dLine[1][1] = dMn_z;

      dStartP[0] = arPmcvUnit[1].dMny;  dStartP[1] = arPmcvUnit[1].dMnz;
      dEndP[0]   = arPmcvUnit[2].dMny;  dEndP[1]   = arPmcvUnit[2].dMnz;

      if(CMathFunc::mathLineSegCross2D(dLine, dStartP, dEndP, dCross) == 0)
      {
        arPmcvUnit.RemoveAll();
        return FALSE;
      }
      double dLength = CMathFunc::mathLength(dStartP[0], dStartP[1], dEndP[0], dEndP[1]);
      dModul_DirID = CMathFunc::mathLength(arPmcvUnit[1].dMny, arPmcvUnit[1].dMnz, dCross[0], dCross[1])/dLength;     

      arPmcvUnit.RemoveAll();
      return TRUE;
    }
  }
  return FALSE;
}

void CMdgnPmTool::Get_ResultPmcvDataID(_UMD_3DPM_DATA& Umd3DPmData, int iAxisDir_Global, int& iPmcvDataID, double& dModul_My, double& dModul_Mz)
{

	int iDivisionNumber_90deg = Umd3DPmData.iDivisionNumber_90deg;
	int iSymmetryType         = Umd3DPmData.iSymmetryType; 
	if(Umd3DPmData.PmcvData.GetSize() < iDivisionNumber_90deg + 1)return;
  //iAxisDir_Global을 0 ~ iDivisionNumber_90deg*4-1 의 범위로 조절
  if(iAxisDir_Global >= iDivisionNumber_90deg*4)
    iAxisDir_Global = iAxisDir_Global%(iDivisionNumber_90deg*4);
  else if(iAxisDir_Global < 0)
  {
    while(iAxisDir_Global < 0)
      iAxisDir_Global += iDivisionNumber_90deg*4;
  }
	
  switch(iSymmetryType)
  {
  case 1:// 비대칭
    iPmcvDataID = iAxisDir_Global;
    dModul_My = 1.0;  dModul_Mz = 1.0;
    return;
  case 2:// 좌우상하대칭시에는 0-90까지만
    if(iAxisDir_Global < iDivisionNumber_90deg)
    {
      iPmcvDataID = iAxisDir_Global;
      dModul_My = 1.0;  dModul_Mz = 1.0;
    }
    else if(iAxisDir_Global < iDivisionNumber_90deg*2)
    {
      iPmcvDataID = iDivisionNumber_90deg*2 - iAxisDir_Global;
      dModul_My = -1.0;  dModul_Mz = 1.0;
    }
    else if(iAxisDir_Global < iDivisionNumber_90deg*3)
    {
      iPmcvDataID = iAxisDir_Global - iDivisionNumber_90deg*2;
      dModul_My = -1.0;  dModul_Mz = -1.0;
    }
    else 
    {
      iPmcvDataID = iDivisionNumber_90deg*4 - iAxisDir_Global;
      dModul_My = 1.0;  dModul_Mz = -1.0;
    }
    return;
  case 3:// z축에 대칭
    if(iAxisDir_Global < iDivisionNumber_90deg*2)
    {
      iPmcvDataID = iAxisDir_Global;
      dModul_My = 1.0;  dModul_Mz = 1.0;
    }
    else 
    {
      iPmcvDataID = iDivisionNumber_90deg*4 - iAxisDir_Global;
      dModul_My = 1.0;  dModul_Mz = -1.0;
    }
    return;
  case 4:// y축에 대칭
    if(iAxisDir_Global < iDivisionNumber_90deg)
    {
      iPmcvDataID = iDivisionNumber_90deg - iAxisDir_Global;
      dModul_My = -1.0;  dModul_Mz = 1.0;
    }
    else if(iAxisDir_Global < iDivisionNumber_90deg*3)
    {
      iPmcvDataID = iAxisDir_Global - iDivisionNumber_90deg;
      dModul_My = 1.0;  dModul_Mz = 1.0;
    }
    else 
    {
      iPmcvDataID = iDivisionNumber_90deg*5 - iAxisDir_Global;
      dModul_My = -1.0;  dModul_Mz = 1.0;
    }
    return;  
  default:// 원점에 대칭(case 5)
    if(iAxisDir_Global < iDivisionNumber_90deg*2)
    {
      iPmcvDataID = iAxisDir_Global;
      dModul_My = 1.0;  dModul_Mz = 1.0;
    }
    else 
    {
      iPmcvDataID = iAxisDir_Global - iDivisionNumber_90deg*2;
      dModul_My = -1.0;  dModul_Mz = -1.0;
    }
    return;
  }
	
}

_UMD_PMCV_UNIT_POLY CMdgnPmTool::Get_CngSignPmcvUnit(_UMD_PMCV_UNIT_POLY& pmcvUnit_In, double dModul_My, double dModul_Mz)
{
  _UMD_PMCV_UNIT_POLY pmcvUnit_Out = pmcvUnit_In;
  pmcvUnit_Out.dMny *= dModul_My;
  pmcvUnit_Out.dMnz *= dModul_Mz;
  return pmcvUnit_Out;
}

double CMdgnPmTool::Get_RotationGlobalAngle(int iDivisionNumber_90deg, int iAxisDir_Global)
{
  return double(iAxisDir_Global)*90.0/double(iDivisionNumber_90deg);
}
double CMdgnPmTool::Get_RotationGlobalAngle(int iDivisionNumber_90deg, int iAxisDir_Global, double dModul_DirID)
{
	return (double(iAxisDir_Global)+dModul_DirID)*90.0/double(iDivisionNumber_90deg);
}

BOOL CMdgnPmTool::Cal_PierceCheckToPolyLine(double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit)
{
		
		double p0[3], p1[3];
		p0[0] = dPs;  p0[1] = dMs_z;  p0[2] = dMs_y;
		p1[0] = dPd;  p1[1] = dMd_z;  p1[2] = dMd_y;
		
		double (*polyLine)[3];  
		int nData = arPmcvUnit.GetSize();
		
		if(nData < 3) return FALSE;
		
		polyLine = new double[nData][3];
		
		_UMD_PMCV_UNIT_POLY pmcvUnit;
		
		for(int i=0 ; i<nData ; i++)
		{
			pmcvUnit = arPmcvUnit.GetAt(i);
			polyLine[i][0] = pmcvUnit.dPn;
			polyLine[i][1] = pmcvUnit.dMnz;
			polyLine[i][2] = pmcvUnit.dMny;
		}
		
		BOOL bCheck = Cal_MathPierceCheckToPolyLine(p0, p1, nData, polyLine);
		delete[] polyLine;
		return bCheck;
}

BOOL CMdgnPmTool::Cal_MathPierceCheckToPolyLine(double p0[3], double p1[3], const int nData, double polyLine[][3])
{
	
  if(nData < 3)return FALSE;
  //원점을 기준으로 하여 새로운 좌표축 생성
  double ux[3], uy[3], uz[3];
  ux[0] = p1[0]-p0[0];  ux[1] = p1[1]-p0[1];  ux[2] = p1[2]-p0[2];
  CMathFunc::mathNormalize(ux, ux);
  uz[0] = 0.0;          uz[1] = -ux[2];       uz[2] = ux[1];
  if(uz[1] == 0.0 && uz[2] == 0.0)
    uz[2] = 1.0;
  CMathFunc::mathNormalize(uz, uz);
  CMathFunc::mathCross(uz, ux, uy);// ux와 uz평면에 직각이 백터는 외적을 이용하여 구함 또한 직각이 단위 벡터이므로 외적의 크기도 단위벡터임
	
  //새로우 좌표축으로 좌표변환
  double (*coor3D)[3];  
  coor3D = new double[nData][3];// 원본 보존을 위해서 복사본 사용
  for(int i=0 ; i<nData ; i++)
  {
    for(int j=0 ; j<3 ; j++)
      coor3D[i][j] = polyLine[i][j]-p0[j];
  }
  CMathFunc::mathTranUCS(ux, uy, uz, nData, coor3D);
	
  // 2D자료로 변환
  double (*coor2D)[2];
  coor2D = new double[nData][2];  
  for(int i=0 ; i<nData ; i++)
  {
    for(int j=0 ; j<2 ; j++)
      coor2D[i][j] = coor3D[i][j+1];
  }  
	
  // 2D 평면상에서 포함여부 체크
  double origin[2];
  origin[0] = origin[1] = 0.0;  
  BOOL bCheck = CMathFunc::mathIsInsidePoint2D(origin, nData, coor2D, TRUE);
	
  delete [] coor3D;
  delete [] coor2D;
	
  return bCheck;
}
BOOL CMdgnPmTool::Cal_MathPierceCheckToPolyLine(double p1[3], const int nData, double polyLine[][3])
{
	double p0[3];
	p0[0] = p0[1] = p0[2] = 0.0;
	return Cal_MathPierceCheckToPolyLine(p0, p1, nData, polyLine);
}

BOOL CMdgnPmTool::Cal_PierceCheckToPolyLine(double dPd, double dMd_y, double dMd_z, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit)
{
		return Cal_PierceCheckToPolyLine(0.0, 0.0, 0.0, dPd, dMd_y, dMd_z, arPmcvUnit);
}

BOOL CMdgnPmTool::Cal_PiercePointToPolyLine(double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit)
{
  if(!Cal_CrossPointToPolyLine(dPs, dMs_y, dMs_z, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit)) return FALSE;
  double vector1[3];
  double vector2[3];
  vector1[0] = dPd-dPs;  vector1[1] = dMd_y-dMs_y;  vector1[2] = dMd_z-dMs_z;
  vector2[0] = dPn-dPs;  vector2[1] = dMn_y-dMs_y;  vector2[2] = dMn_z-dMs_z;
	
  double angle = CMathFunc::mathCrossAngle(vector1, vector2);
  if(angle < 0.0001) return TRUE;
  return FALSE;
}
BOOL CMdgnPmTool::Cal_PiercePointToPolyLine(double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit)
{
	return Cal_PiercePointToPolyLine(0.0, 0.0, 0.0, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit);
}

BOOL CMdgnPmTool::Cal_CrossPointToPolyLine(double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit)
{
  dPn = dMn_y = dMn_z = 0.0;
  if(arPmcvUnit.GetSize() != 3) return FALSE;
  // [0]:Mz [1]:My [2]:P 
  double lineVector[3];
  double linePoint[3];
  double planeNormal[3];
  double planePoint[3];
  double pInts[3];
  double vector1[3];
  double vector2[3];
  lineVector[0] = dMd_z-dMs_z;  lineVector[1] = dMd_y-dMs_y;  lineVector[2] = dPd-dPs;
  linePoint[0]  = dMs_z;        linePoint[1]  = dMs_y;        linePoint[2]  = dPs;
  planePoint[0] = arPmcvUnit[0].dMnz;  planePoint[1] = arPmcvUnit[0].dMny;  planePoint[2] = arPmcvUnit[0].dPn;
  vector1[0] = arPmcvUnit[1].dMnz - arPmcvUnit[0].dMnz;
  vector1[1] = arPmcvUnit[1].dMny - arPmcvUnit[0].dMny;
  vector1[2] = arPmcvUnit[1].dPn  - arPmcvUnit[0].dPn;
  vector2[0] = arPmcvUnit[2].dMnz - arPmcvUnit[0].dMnz;
  vector2[1] = arPmcvUnit[2].dMny - arPmcvUnit[0].dMny;
  vector2[2] = arPmcvUnit[2].dPn  - arPmcvUnit[0].dPn;
  CMathFunc::mathCross(vector1, vector2, planeNormal);
	
  if(!CMathFunc::mathIntersectPointToPlane(lineVector, linePoint, planeNormal, planePoint, pInts)) return FALSE;
  dPn = pInts[2]; dMn_y = pInts[1]; dMn_z = pInts[0];
	
	double dLang, dLang0, dLang1, dLang2, dModul0, dModul1, dModul2;
	dLang0 = CMathFunc::mathLength(dPn, dMn_y, dMn_z, arPmcvUnit[0].dPn, arPmcvUnit[0].dMny, arPmcvUnit[0].dMnz);
	dLang1 = CMathFunc::mathLength(dPn, dMn_y, dMn_z, arPmcvUnit[1].dPn, arPmcvUnit[1].dMny, arPmcvUnit[1].dMnz);
	dLang2 = CMathFunc::mathLength(dPn, dMn_y, dMn_z, arPmcvUnit[2].dPn, arPmcvUnit[2].dMny, arPmcvUnit[2].dMnz);
	dModul0 = dModul1 = dModul2 = 1.0/3.0;
	dLang = dLang0+dLang1+dLang2;
	if(dLang != 0.0)
	{
		dModul0 = (dLang1+dLang2)/(2.0*dLang);
		dModul1 = (dLang0+dLang2)/(2.0*dLang);
		dModul2 = (dLang0+dLang1)/(2.0*dLang);
	}
	desi = dModul0*arPmcvUnit[0].desi + dModul1*arPmcvUnit[1].desi + dModul2*arPmcvUnit[2].desi;
	dXn  = dModul0*arPmcvUnit[0].dXn  + dModul1*arPmcvUnit[1].dXn  + dModul2*arPmcvUnit[2].dXn;
	
  return TRUE;
}
BOOL CMdgnPmTool::Cal_CrossPointToPolyLine(double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit)
{// 원점을 기준으로 계산
	return Cal_CrossPointToPolyLine(0.0, 0.0, 0.0, dPd, dMd_y, dMd_z, dPn, dMn_y, dMn_z, desi, dXn, arPmcvUnit);
}

BOOL CMdgnPmTool::Get_UnitVectorToMzMyP(double dP, double dMy, double dMz, double vectorMzMyP[3], double vectorMP[2], double vectorMzMy[2])
{
  vectorMzMyP[0] = vectorMzMy[0] = dMz;
  vectorMzMyP[1] = vectorMzMy[1] = dMy;
  vectorMzMyP[2] = vectorMP[1] = dP;
  vectorMP[0] = CMathFunc::mathLength(dMz, dMy);
	
  BOOL bCheck = TRUE;
  if(!CMathFunc::mathNormalize(vectorMzMyP, vectorMzMyP)) bCheck = FALSE;
  if(!CMathFunc::mathNormalize2D(vectorMzMy, vectorMzMy)) bCheck = FALSE;
  if(!CMathFunc::mathNormalize2D(vectorMP, vectorMP)) bCheck = FALSE;
  return bCheck;
}

void CMdgnPmTool::Get_AngleMinIDToPmCurve_NonScale(_UMD_3DPM_DATA& Umd3DPmData, double nMzMyP_Org[3], double nVectorMzMyP[3], int iAxisDir_Global, double dPmax, double dMmax, int& iUnitID, double& dAngleMin)
{
	int iDivisionNumber_90deg = Umd3DPmData.iDivisionNumber_90deg;
  if(Umd3DPmData.PmcvData.GetSize() < iDivisionNumber_90deg + 1)return;
	
  iUnitID = 0 ; dAngleMin = 180.0;
  int iPmcvDataID, iPmcvSize;
  double dModul_My, dModul_Mz;
  double dAngleN;
  double vectorN[3], vectorD[3];
  _UMD_PMCV_DATA_POLY pmcvData;
  _UMD_PMCV_UNIT_POLY pmcvUnit;
	
  Get_ResultPmcvDataID(Umd3DPmData, iAxisDir_Global, iPmcvDataID, dModul_My, dModul_Mz);
  pmcvData = Umd3DPmData.PmcvData.GetAt(iPmcvDataID);
  iPmcvSize = pmcvData.arPmcvUnit.GetCount();
	
  vectorD[0] = nVectorMzMyP[0];
  vectorD[1] = nVectorMzMyP[1];
  vectorD[2] = nVectorMzMyP[2];
  
  for(int i=1 ; i<=iPmcvSize ; i++)
  {
    pmcvUnit = pmcvData.arPmcvUnit[i];
    vectorN[0] = pmcvUnit.dMnz*dModul_Mz/dMmax - nMzMyP_Org[0];  
		vectorN[1] = pmcvUnit.dMny*dModul_My/dMmax - nMzMyP_Org[1];  
		vectorN[2] = pmcvUnit.dPn/dPmax            - nMzMyP_Org[2];
		
    CMathFunc::mathNormalize(vectorN, vectorN);
    dAngleN = CMathFunc::mathCrossAngle(vectorD, vectorN);
    if(fabs(dAngleN) < fabs(dAngleMin))
    {
      iUnitID = i;  dAngleMin = dAngleN;
    }
    //else
    //  break;
  }
  return;
}