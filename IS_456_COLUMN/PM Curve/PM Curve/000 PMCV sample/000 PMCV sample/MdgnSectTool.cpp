// MdgnSectTool.cpp: implementation of the CMdgnSectTool class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Pmcv.h"
#include "MdgnSectTool.h"
#include "MathFunc.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CMdgnSectTool::CMdgnSectTool()
{

}

CMdgnSectTool::~CMdgnSectTool()
{

}

BOOL CMdgnSectTool::Get_RotateSectData(int nType, double center_y, double center_z, int nSliceNum, int nDivisionNumber_90deg, 
																			 _UMD_RC_CON_POLY_DATA& PolyData, _UMD_RC_RBAR& RbarData, _UMD_RC_COL_ROTATE_SECT& RotSectData)
{
	RotSectData.Init();
	
	int iSymmetryType = Get_CheckedSymmetryType(PolyData, RbarData, center_y, center_z);
	int iHoldNumberToUnitList = Get_HoldNumberToUnitList(iSymmetryType, nDivisionNumber_90deg);
  double rotation = 0;
	
	RotSectData.nDataType             = nType;
	RotSectData.iDivisionNumber_90deg = nDivisionNumber_90deg;
	RotSectData.iSymmetryType         = iSymmetryType;
	RotSectData.arConUnitList.SetSize(iHoldNumberToUnitList);
	RotSectData.arRbarUnitList.SetSize(iHoldNumberToUnitList);
  
	for(int i=0; i<iHoldNumberToUnitList; i++)
	{
    rotation = Get_RotationAngle(i, iSymmetryType, nDivisionNumber_90deg);// 단위 : deg
    
		_UMD_RC_CON_PART_LIST RotConData;
		_UMD_RC_RBAR_LIST RotRbarData;
		
    if(!Get_RotateConUnitPartList(nType, center_y, center_z, rotation, nSliceNum, PolyData, RotConData)) return FALSE;
		if(!Get_RotateRebarData(center_y, center_z, rotation, RbarData, RotRbarData)) return FALSE;
    
    RotSectData.arConUnitList.SetAt(i,RotConData);
		RotSectData.arRbarUnitList.SetAt(i,RotRbarData);
	}	
	
	return TRUE;
}

int CMdgnSectTool::Get_CheckedSymmetryType(_UMD_RC_CON_POLY_DATA& ConPolyData, _UMD_RC_RBAR& RbarUnitData, double dCenter_y, double dCenter_z)
{
	double dSymmetry_Zero = min(1.0E-05, cUMDRC_Zero);//*100.0);
  // 단면 대칭 검토 //////////////////////////////////////////////////////////////
  CPolyMaker  PolyMaker;	    
	_UMD_RC_RBAR_LIST arRbarData; //CMap형식의 철근 정보로 CArray형시으로 바꾸어서 진행
  
  if(ConPolyData.OutPoly.aVertex.GetSize() < 1)
    return FALSE;

  // 최대 최소 초기화
  double max_y, min_y, max_z, min_z;
  max_y = min_y = dCenter_y;  max_z = min_z = dCenter_z;

  // 단면정보를 I_PolyMaker형식의 pPolyMaker로 변경
	// 비회전된 단면 정보로 PolyMaker정보 구성
	if(!Get_RotateConPolyData(dCenter_y, dCenter_z, 0.0, ConPolyData, max_y, min_y, max_z, min_z, &PolyMaker))
		return FALSE;
	// 비회전된 철근 정보로 arRbarData를 채움
	if(!Get_RotateRebarData(dCenter_y, dCenter_z, 0.0, RbarUnitData, arRbarData))
		return FALSE;
	

  double dTop_My=0.0,   dTop_Mz=0.0,      dBottom_My=0.0, dBottom_Mz=0.0;// 단면을 상하로 나누었을때 분리된 Polygon들의 각각의 면적과 도심과의 거리에 대한 각 축별 단면 1차모멘트 
  double dLeft_My=0.0,  dLeft_Mz=0.0,     dRight_My =0.0, dRight_Mz=0.0;// 단면을 좌우로 나누었을때 분리된 Polygon들의 각각의 면적과 도심과의 거리에 대한 각 축별 단면 1차모멘트
  double dTop_Area=0.0, dBottom_Area=0.0, dLeft_Area=0.0, dRight_Area=0.0;// 분리된 단면의 각각의 단면적  
  double dTotal_Area=0.0;
  // 각 회전 중심점을 기준으로 상하 조각내기	

  CArray<DGL_Line3d,DGL_Line3d&> arSlicingLines;
  arSlicingLines.SetSize(3);
	DGL_Line3d LineGL;
  LineGL.m_P1.Set(max_y+cUMDRC_Zero, max_z+cUMDRC_Zero, 0.0);	// Sta.
	LineGL.m_P2.Set(min_y-cUMDRC_Zero, max_z+cUMDRC_Zero, 0.0);	// End.
	arSlicingLines.SetAt(0,LineGL);
	LineGL.m_P1.Set(max_y+cUMDRC_Zero, dCenter_z, 0.0);	// Sta.
	LineGL.m_P2.Set(min_y-cUMDRC_Zero, dCenter_z, 0.0);	// End.
	arSlicingLines.SetAt(1,LineGL);
  LineGL.m_P1.Set(max_y+cUMDRC_Zero, min_z-cUMDRC_Zero, 0.0);	// Sta.
	LineGL.m_P2.Set(min_y-cUMDRC_Zero, min_z-cUMDRC_Zero, 0.0);	// End.
	arSlicingLines.SetAt(2,LineGL);
	
	PolyMaker.SetSlicingLines(arSlicingLines);

	if(!PolyMaker.MakeSlicedPolygons())	ASSERT(0);

	int iResPolySize = PolyMaker.GetResultPolyCount();
  for(int i=0 ; i<iResPolySize ; i++)
  {
    CArray<DGL_3dp,DGL_3dp&> arResPoly;
		PolyMaker.GetResultPoly(i,arResPoly);
		int iResNum = arResPoly.GetSize();
		// Calculate Center(x,y), Area.
		//---------------------------------
		double* dpx = new double[iResNum];
  	double* dpy = new double[iResNum];
		//---------------------------------
		for(int k=0; k<iResNum; k++)
		{
			dpx[k] = arResPoly[k].x();
			dpy[k] = arResPoly[k].y();
		}
		double dxCen=0.0, dyCen=0.0, dArea=0.0;
		if(!CMathFunc::mathPolyCentroid(iResNum,dpx,dpy,dxCen,dyCen,dArea))	ASSERT(0);
		//---------------------------------
		delete []dpx;
		delete []dpy;
		//---------------------------------    
    if(dyCen > dCenter_z)
    {
      dTop_Area += dArea;
      dTop_My   += dArea*(dyCen-dCenter_z);
      dTop_Mz   += dArea*(dxCen-dCenter_y);
    }
    else
    {
      dBottom_Area += dArea;
      dBottom_My   += dArea*(dyCen-dCenter_z);
      dBottom_Mz   += dArea*(dxCen-dCenter_y);
    }
  }
  // 각 회전 중심점을 기준으로 좌우 조각내기		
  LineGL.m_P1.Set(max_y+cUMDRC_Zero, max_z+cUMDRC_Zero, 0.0);	// Sta.
	LineGL.m_P2.Set(max_y+cUMDRC_Zero, min_z-cUMDRC_Zero, 0.0);	// End.
	arSlicingLines.SetAt(0,LineGL);
	LineGL.m_P1.Set(dCenter_y, max_z+cUMDRC_Zero, 0.0);	// Sta.
	LineGL.m_P2.Set(dCenter_y, min_z-cUMDRC_Zero, 0.0);	// End.
	arSlicingLines.SetAt(1,LineGL);
  LineGL.m_P1.Set(min_y-cUMDRC_Zero, max_z+cUMDRC_Zero, 0.0);	// Sta.
	LineGL.m_P2.Set(min_y-cUMDRC_Zero, min_z-cUMDRC_Zero, 0.0);	// End.
	arSlicingLines.SetAt(2,LineGL);
	
	PolyMaker.SetSlicingLines(arSlicingLines);

	if(!PolyMaker.MakeSlicedPolygons())	ASSERT(0);
  
	iResPolySize = PolyMaker.GetResultPolyCount();
  for(int i=0 ; i<iResPolySize ; i++)
  {
    CArray<DGL_3dp,DGL_3dp&> arResPoly;
		PolyMaker.GetResultPoly(i,arResPoly);
		int iResNum = arResPoly.GetSize();
		// Calculate Center(x,y), Area.
		//---------------------------------
		double* dpx = new double[iResNum];
		double* dpy = new double[iResNum];
		//---------------------------------
		for(int k=0; k<iResNum; k++)
		{
			dpx[k] = arResPoly[k].x();
			dpy[k] = arResPoly[k].y();
		}
		double dxCen=0.0, dyCen=0.0, dArea=0.0;
		if(!CMathFunc::mathPolyCentroid(iResNum,dpx,dpy,dxCen,dyCen,dArea))	ASSERT(0);
		//---------------------------------
		delete []dpx;
		delete []dpy;
		//---------------------------------
    if(dxCen > dCenter_y)
    {
      dRight_Area += dArea;
      dRight_My   += dArea*(dyCen-dCenter_z);
      dRight_Mz   += dArea*(dxCen-dCenter_y);
    }
    else
    {
      dLeft_Area += dArea;
      dLeft_My   += dArea*(dyCen-dCenter_z);
      dLeft_Mz   += dArea*(dxCen-dCenter_y);
    }
  }
  dTotal_Area = dTop_Area + dBottom_Area;
  BOOL bSection_Symmetry_YAxis=FALSE, bSection_Symmetry_ZAxis=FALSE, bSection_Symmetry_Org=FALSE;// 각각 단면의 Y축 대칭여부, Z축 대칭여부, 원점 대칭여부 
  // 여기서 비교항목은 분리된 Polygon의 단면 1차모멘트와 단면적을 비교함(기본적으로 도심을 기준으로 하므로 단면1차모멘트항목은 대칭으로 나타남)
  if(fabs((dTop_My+dBottom_My)/dTotal_Area/(max_z-min_z)) < dSymmetry_Zero && fabs((dTop_Mz-dBottom_Mz)/dTotal_Area/(max_y-min_y)) < dSymmetry_Zero && fabs((dTop_Area-dBottom_Area)/dTotal_Area) < dSymmetry_Zero)
    bSection_Symmetry_YAxis = TRUE;
  if(fabs((dRight_My-dLeft_My)/dTotal_Area/(max_z-min_z)) < dSymmetry_Zero && fabs((dRight_Mz+dLeft_Mz)/dTotal_Area/(max_y-min_y)) < dSymmetry_Zero && fabs((dLeft_Area-dRight_Area)/dTotal_Area) < dSymmetry_Zero)
    bSection_Symmetry_ZAxis = TRUE;
  if(fabs((dTop_My+dBottom_My)/dTotal_Area/(max_z-min_z)) < dSymmetry_Zero && fabs((dTop_Mz+dBottom_Mz)/dTotal_Area/(max_y-min_y)) < dSymmetry_Zero 
    && fabs((dRight_My+dLeft_My)/dTotal_Area/(max_z-min_z)) < dSymmetry_Zero && fabs((dRight_Mz+dLeft_Mz)/dTotal_Area/(max_y-min_y)) < dSymmetry_Zero 
    && fabs((dTop_Area-dBottom_Area)/dTotal_Area) < dSymmetry_Zero && fabs((dLeft_Area-dRight_Area)/dTotal_Area) < dSymmetry_Zero)
    bSection_Symmetry_Org   = TRUE;

  if(bSection_Symmetry_YAxis==FALSE && bSection_Symmetry_ZAxis==FALSE && bSection_Symmetry_Org==FALSE)
	{ PolyMaker.ResetAll();  return 1; } //비대칭임 철근을 검토할 필요없음

  
  // 철근 대칭 검토 //////////////////////////////////////////////////////////////
  //초기화
  _UMD_RC_RBAR_UNIT RbarUnit;
  dTop_My=0.0;   dTop_Mz=0.0;      dBottom_My=0.0; dBottom_Mz=0.0;
  dLeft_My=0.0;  dLeft_Mz=0.0;     dRight_My =0.0; dRight_Mz=0.0;
  dTop_Area=0.0; dBottom_Area=0.0; dLeft_Area=0.0; dRight_Area=0.0;  
  dTotal_Area=0.0;
  int iRebarUnitSize = arRbarData.List.GetSize();
  for(int i=0 ; i<iRebarUnitSize ; i++)
  {
    RbarUnit.Init();
		RbarUnit = arRbarData.List.GetAt(i);

    if(fabs((RbarUnit.dyz[1]-dCenter_z)/(max_z-min_z)) < cUMDRC_Zero)
    {
      dTop_Area    += RbarUnit.dArea/2.0;
      dBottom_Area += RbarUnit.dArea/2.0;
      //모멘트는 0이므로 더할필요 없음
    }
    else if(RbarUnit.dyz[1] > dCenter_z)
    {
      dTop_Area    += RbarUnit.dArea;
      dTop_My      += RbarUnit.dArea*(RbarUnit.dyz[1]-dCenter_z);
      dTop_Mz      += RbarUnit.dArea*(RbarUnit.dyz[0]-dCenter_y);
    }
    else
    {
      dBottom_Area += RbarUnit.dArea;
      dBottom_My   += RbarUnit.dArea*(RbarUnit.dyz[1]-dCenter_z);
      dBottom_Mz   += RbarUnit.dArea*(RbarUnit.dyz[0]-dCenter_y);
    }

    if(fabs((RbarUnit.dyz[0]-dCenter_y)/(max_y-min_y)) < cUMDRC_Zero)
    {
      dRight_Area  += RbarUnit.dArea/2.0;
      dLeft_Area   += RbarUnit.dArea/2.0;
      //모멘트는 0이므로 더할필요 없음
    }
    else if(RbarUnit.dyz[0] > dCenter_y)
    {
      dRight_Area  += RbarUnit.dArea;
      dRight_My    += RbarUnit.dArea*(RbarUnit.dyz[1]-dCenter_z);
      dRight_Mz    += RbarUnit.dArea*(RbarUnit.dyz[0]-dCenter_y);
    }
    else
    {
      dLeft_Area   += RbarUnit.dArea;
      dLeft_My     += RbarUnit.dArea*(RbarUnit.dyz[1]-dCenter_z);
      dLeft_Mz     += RbarUnit.dArea*(RbarUnit.dyz[0]-dCenter_y);
    }
  }
  dTotal_Area = dTop_Area + dBottom_Area;
  BOOL bRebar_Symmetry_YAxis=FALSE, bRebar_Symmetry_ZAxis=FALSE, bRebar_Symmetry_Org=FALSE;// 각각 철근의 Y축 대칭여부, Z축 대칭여부, 원점 대칭여부 
  if(fabs(dTotal_Area) < cUMDRC_Zero)
  {// 철근이 없으면 비교대상 제외
    bRebar_Symmetry_YAxis=TRUE; bRebar_Symmetry_ZAxis=TRUE; bRebar_Symmetry_Org=TRUE;
  }
  else
  {
    // 여기서 비교항목은 분리된 철근들의 단면의 단면 1차모멘트와 단면적을 비교함
    if(fabs((dTop_My+dBottom_My)/dTotal_Area/(max_z-min_z)) < dSymmetry_Zero && fabs((dTop_Mz-dBottom_Mz)/dTotal_Area/(max_y-min_y)) < dSymmetry_Zero && fabs((dTop_Area-dBottom_Area)/dTotal_Area) < dSymmetry_Zero)
      bRebar_Symmetry_YAxis = TRUE;
    if(fabs((dRight_My-dLeft_My)/dTotal_Area/(max_z-min_z)) < dSymmetry_Zero && fabs((dRight_Mz+dLeft_Mz)/dTotal_Area/(max_y-min_y)) < dSymmetry_Zero && fabs((dLeft_Area-dRight_Area)/dTotal_Area) < dSymmetry_Zero)
      bRebar_Symmetry_ZAxis = TRUE;
    if(fabs((dTop_My+dBottom_My)/dTotal_Area/(max_z-min_z)) < dSymmetry_Zero && fabs((dTop_Mz+dBottom_Mz)/dTotal_Area/(max_y-min_y)) < dSymmetry_Zero 
    && fabs((dRight_My+dLeft_My)/dTotal_Area/(max_z-min_z)) < dSymmetry_Zero && fabs((dRight_Mz+dLeft_Mz)/dTotal_Area/(max_y-min_y)) < dSymmetry_Zero 
    && fabs((dTop_Area-dBottom_Area)/dTotal_Area) < dSymmetry_Zero && fabs((dLeft_Area-dRight_Area)/dTotal_Area) < dSymmetry_Zero)
      bRebar_Symmetry_Org   = TRUE;
  }

	PolyMaker.ResetAll();
  if(bSection_Symmetry_YAxis==TRUE && bSection_Symmetry_ZAxis==TRUE && bSection_Symmetry_Org==TRUE && bRebar_Symmetry_YAxis==TRUE && bRebar_Symmetry_ZAxis==TRUE && bRebar_Symmetry_Org==TRUE)
    return 2; //상하 좌우 대칭
  if(bSection_Symmetry_ZAxis==TRUE && bRebar_Symmetry_ZAxis==TRUE)
    return 3; //좌우 대칭
  if(bSection_Symmetry_YAxis==TRUE && bRebar_Symmetry_YAxis==TRUE)
    return 4; //상하 대칭
  if(bSection_Symmetry_Org==TRUE && bRebar_Symmetry_Org==TRUE)
    return 5; //원점 대칭
  
  return 1;   //비대칭
}

BOOL CMdgnSectTool::Get_RotateRebarData(double center_y, double center_z, double rotation, _UMD_RC_RBAR& RbarData, _UMD_RC_RBAR_LIST& RotateRbarData)
{
	RotateRbarData.Init();
	
	int nRbarSize = RbarData.arRbarUnit.GetCount();
	if(nRbarSize <= 0) return FALSE;
  RotateRbarData.List.SetSize(nRbarSize);
	
	_UMD_RC_RBAR_UNIT RbarUnit, RbarUnit_Rot;
	int iIndex=0;
	POSITION Pos = RbarData.arRbarUnit.GetStartPosition();
	int nCount = 0;
	while(Pos)
	{
		RbarUnit.Init();
		RbarData.arRbarUnit.GetNextAssoc(Pos, iIndex, RbarUnit);
		
		double dx = RbarUnit.dyz[0];
		double dy = RbarUnit.dyz[1];
    double dz = 0.0;
		
    // 중립축을 반시계 방향으로 돌릴려면 단면좌표를 시계방향으로 돌려야함
    CMathFunc::mathRotate(-rotation, center_y, center_z, 0.0, 0.0, 0.0, 1.0, dx, dy, dz); 
		
    RbarUnit_Rot.Init();
		RbarUnit_Rot = RbarUnit;
    RbarUnit_Rot.dyz[0] = dx;
    RbarUnit_Rot.dyz[1] = dy;
		
		RotateRbarData.List.SetAt(nCount, RbarUnit_Rot);
		nCount++;
		if(nCount >= nRbarSize)
			break;
	}		  
  return TRUE;
}

BOOL CMdgnSectTool::Get_RotateConPolyData(double center_y, double center_z, double rotation, _UMD_RC_CON_POLY_DATA& PolyData, double& max_y, double& min_y, double& max_z, double& min_z, CPolyMaker* pPolyMaker)
{	
	// 최대 최소 초기화
  max_y = min_y = center_y;
  max_z = min_z = center_z;
	
  // Outer.
  CArray<DGL_3dp,DGL_3dp&> arOutCut;
	
	_UMD_RC_GSEC_VERTEX PontUnit;
  for(int k=0; k<PolyData.OutPoly.aVertex.GetSize(); k++)
  {
		PontUnit = PolyData.OutPoly.aVertex.GetAt(k);
    double dx = PontUnit.dX;
		double dy = PontUnit.dY;
    double dz = 0.0;
		
    // 중립축을 반시계 방향으로 돌릴려면 단면좌표를 시계방향으로 돌려야함
    CMathFunc::mathRotate(-rotation, center_y, center_z, 0.0, 0.0, 0.0, 1.0, dx, dy, dz); 
		
    if(k==0)
    {// 첫번째 좌표값으로 최대 초소 초기화
      max_y = min_y = dx;
      max_z = min_z = dy;
    }
    else
    {// 두번째 좌표값부터 외각선의 최대 최소값을 찾음
      if(dx > max_y)        max_y = dx;
      if(dx < min_y)        min_y = dx;
      if(dy > max_z)        max_z = dy;
      if(dy < min_z)        min_z = dy;      
    }
		
    DGL_3dp OutXGL;
		OutXGL.Set(dx,dy,0.0);
		arOutCut.Add(OutXGL);	
  }
  pPolyMaker->SetOutterPoly(arOutCut);
	
  // Inner.
  for(int k=0; k<PolyData.arInnPoly.GetSize(); k++)
  {
		CArray<DGL_3dp,DGL_3dp&> arInnCut;
		int nVertexSize = PolyData.arInnPoly[k].aVertex.GetSize();
		for(int j=0; j<nVertexSize; j++)
		{
			PontUnit = PolyData.arInnPoly[k].aVertex.GetAt(j);
			double dx = PontUnit.dX;
			double dy = PontUnit.dY;
      double dz = 0.0;
			
      // 중립축을 반시계 방향으로 돌릴려면 단면좌표를 시계방향으로 돌려야함
      CMathFunc::mathRotate(-rotation, center_y, center_z, 0.0, 0.0, 0.0, 1.0, dx, dy, dz); 
			
			DGL_3dp InnXGL;
			InnXGL.Set(dx,dy,0.0);
			arInnCut.Add(InnXGL);			
		}
		pPolyMaker->AddInnerPoly(arInnCut);    
  }  
	return TRUE;
}

BOOL CMdgnSectTool::Get_RotateConUnitPartList(int nType, double center_y, double center_z, double rotation, int nSliceNum, _UMD_RC_CON_POLY_DATA& PolyData, _UMD_RC_CON_PART_LIST& arConUnit)
{
	arConUnit.Init();

	CPolyMaker PolyMaker;

	_UMD_RC_CON_PART_UNIT PartUnit;
  
  // 회전된 단면 정보로 PolyMaker정보 구성
	double max_y, min_y, max_z, min_z;
	if(!Get_RotateConPolyData(center_y, center_z, rotation, PolyData, max_y, min_y, max_z, min_z, &PolyMaker))
		return FALSE;
  if(nType == 2)
	{// nType 1: 슬라이스 정보 보관 하지 않음  2:슬라이스 정보 보관
		double dyDim, dzDim, dPartH;
		dyDim = max_y - min_y;
		dzDim = max_z - min_z;
		dPartH = dzDim/(double)nSliceNum;

		CArray<DGL_Line3d,DGL_Line3d&> arSlicingLines;
		for(int j=0; j<=nSliceNum; j++)
		{
			// For Tolerance.
			double dRat	 = (j==0 ? cUMDRC_Zero : j/(double)nSliceNum);
			double dxSta = min_y - cUMDRC_Zero;
			double dxEnd = max_y + cUMDRC_Zero;        
			double dySta = min_z + dRat*dzDim;
			double dyEnd = min_z + dRat*dzDim;
			DGL_Line3d LineGL;
			LineGL.m_P1.Set(dxSta,dySta,0.0);	// Sta.
			LineGL.m_P2.Set(dxEnd,dyEnd,0.0);	// End.
			arSlicingLines.Add(LineGL);
		}
		PolyMaker.SetSlicingLines(arSlicingLines);

		if(!PolyMaker.MakeSlicedPolygons())	ASSERT(0);

		// Set Size.
		int iResPolySize = PolyMaker.GetResultPolyCount();

		arConUnit.List.SetSize(iResPolySize);

		int iResPolyCount=0;

		// For Debugging.
		int iDebugCount=0;
		// Result Polygon.
		for(int j=0; j<iResPolySize; j++)
		{				
			CArray<DGL_3dp,DGL_3dp&> arResPoly;
			PolyMaker.GetResultPoly(j,arResPoly);
			int iResNum = arResPoly.GetSize();
			// Calculate Center(x,y), Area.
			//---------------------------------
			double* dpx = new double[iResNum];
			double* dpy = new double[iResNum];
			//---------------------------------
			for(int k=0; k<iResNum; k++)
			{
				dpx[k] = arResPoly[k].x();
				dpy[k] = arResPoly[k].y();
			}
			double dxCen=0.0, dyCen=0.0, dArea=0.0;
			if(!CMathFunc::mathPolyCentroid(iResNum,dpx,dpy,dxCen,dyCen,dArea))
			{// SHIN ('06.02.26) 모서리에 걸릴경우 3개의 좌표값이 동일한 경우가 생기는데 이때 Debug시 ASSERT(0)의 호출을 피하기 위해
				BOOL bXCheck, bYCheck;
				bXCheck = bYCheck = TRUE;
				for(int k=0; k<iResNum-1; k++)
				{
					if(fabs(dpx[k]-dpx[k+1]) > cUMDRC_Zero)
					{
						bXCheck = FALSE; break;
					}
					if(fabs(dpy[k]-dpy[k+1]) > cUMDRC_Zero)
					{
						bYCheck = FALSE; break;
					}
				}
				if((fabs(dArea) > cUMDRC_Zero) || bXCheck == FALSE || bYCheck == FALSE)  ASSERT(0);
			}
			//---------------------------------
			delete []dpx;
			delete []dpy;
			//---------------------------------
			PartUnit.Init();
			PartUnit.dyz[0] = dxCen;
			PartUnit.dyz[1] = dyCen;
			PartUnit.dArea	= dArea;
			arConUnit.List.SetAt(iResPolyCount, PartUnit);
			iResPolyCount++;		    
		}
	}
	arConUnit.dYmax = max_y;
	arConUnit.dYmin = min_y;
	arConUnit.dZmax = max_z;
	arConUnit.dZmin = min_z;

	PolyMaker.ResetAll();
  return TRUE;
}


int CMdgnSectTool::Get_HoldNumberToUnitList(int iSymmetryType, int nDivisionNumber_90deg)
{
  if(iSymmetryType == 2)// 좌우상하대칭일 경에는 0-90까지만
    return nDivisionNumber_90deg+1;// +1은 90부분을 포함하기 위해서임
  else                    // 기타의 경우 0-180까지(나머지 부분은 좌표만 변환 하면 되므로
    return nDivisionNumber_90deg*2;
}

double CMdgnSectTool::Get_RotationAngle(int iRot, int iSymmetryType, int nDivisionNumber_90deg)
{
  // SHIN (06.01.27) 가로 세로가 차이 날경우 해석결과가 한쪽을 치우칠수 있으므로 보완방법 추가 요망
  if(iSymmetryType == 4)// 상하대칭일때에는 90도 부터 시작
    return double(iRot)*90.0/double(nDivisionNumber_90deg) + 90.0;
  else                    // 나머지의 경우 0도 부터 시작
    return double(iRot)*90.0/double(nDivisionNumber_90deg);
}

BOOL CMdgnSectTool::Get_PolyDataCentroid(double dZ1, double dZ2, double dYmax, double dYmin, CPolyMaker* pPolyMaker, double& center_y, double& center_z, double& dArea_Inner)
{
	center_y = center_z = dArea_Inner = 0.0;

	double dZ_Top = max(dZ1, dZ2);
	double dZ_Bot = min(dZ1, dZ2);
	CArray<DGL_Line3d,DGL_Line3d&> arSlicingLines;
	arSlicingLines.SetSize(2);
	DGL_Line3d LineGL;
	LineGL.m_P1.Set(dYmax + cUMDRC_Zero*cUMDRC_Zero,dZ_Top,0.0);	// Sta.
	LineGL.m_P2.Set(dYmin - cUMDRC_Zero,dZ_Top,0.0);	// End.
	arSlicingLines.SetAt(0, LineGL);
	LineGL.m_P1.Set(dYmax + cUMDRC_Zero*cUMDRC_Zero,dZ_Bot,0.0);	// Sta.
	LineGL.m_P2.Set(dYmin - cUMDRC_Zero,dZ_Bot,0.0);	// End.
	arSlicingLines.SetAt(1, LineGL);


	CPolyMaker PolyMakerUnit;
	CArray<DGL_3dp,DGL_3dp&> arOutCut;
	PolyMakerUnit.RemoveAllCutPolys();
	PolyMakerUnit.RemoveAllResult();
	
	//외곽 
	pPolyMaker->GetOutterPoly(arOutCut);
	PolyMakerUnit.SetOutterPoly(arOutCut);

	PolyMakerUnit.SetSlicingLines(arSlicingLines);

	if(!PolyMakerUnit.MakeSlicedPolygons())	ASSERT(0);

	// Set Size.
	int iResPolySize = PolyMakerUnit.GetResultPolyCount();

	// For Debugging.
	int iDebugCount=0;
	BOOL bInChk;
	// Result Polygon.
	for(int i=0; i<iResPolySize; i++)
	{				
		CArray<DGL_3dp,DGL_3dp&> arResPoly;
		PolyMakerUnit.GetResultPoly(i,arResPoly);
		int iResNum = arResPoly.GetSize();
		// Calculate Center(x,y), Area.
		//---------------------------------
		double* dpx = new double[iResNum];
		double* dpy = new double[iResNum];
		//---------------------------------
		bInChk = TRUE;
		for(int j=0; j<iResNum; j++)
		{
			dpx[j] = arResPoly[j].x();
			dpy[j] = arResPoly[j].y();
			if(dpy[j] > dZ_Top+cUMDRC_Zero || dpy[j] < dZ_Bot-cUMDRC_Zero)
				bInChk = FALSE;
		}
		if(!bInChk) continue;

		double dxCen=0.0, dyCen=0.0, dArea=0.0;
		if(!CMathFunc::mathPolyCentroid(iResNum,dpx,dpy,dxCen,dyCen,dArea))
		{// SHIN ('06.02.26) 모서리에 걸릴경우 3개의 좌표값이 동일한 경우가 생기는데 이때 Debug시 ASSERT(0)의 호출을 피하기 위해
			BOOL bXCheck, bYCheck;
			bXCheck = bYCheck = TRUE;
			for(int j=0; j<iResNum-1; j++)
			{
				if(fabs(dpx[j]-dpx[j+1]) > cUMDRC_Zero)
				{
					bXCheck = FALSE; break;
				}
				if(fabs(dpy[j]-dpy[j+1]) > cUMDRC_Zero)
				{
					bYCheck = FALSE; break;
				}
			}
			if((fabs(dArea) > cUMDRC_Zero) || bXCheck == FALSE || bYCheck == FALSE)  ASSERT(0);
		}
		//---------------------------------
		delete []dpx;
		delete []dpy;
		//---------------------------------
		
		center_y += dArea*dxCen;
		center_z += dArea*dyCen;
		dArea_Inner += dArea;   
	}
	PolyMakerUnit.RemoveAllCutPolys();
	PolyMakerUnit.RemoveAllResult();

	
	int nInnPolySize = pPolyMaker->GetInnerPolyCount();
	for(int i=0 ; i<nInnPolySize ; i++)
	{
		pPolyMaker->GetInnerPoly(i, arOutCut);
		PolyMakerUnit.SetOutterPoly(arOutCut);

		PolyMakerUnit.SetSlicingLines(arSlicingLines);

		if(!PolyMakerUnit.MakeSlicedPolygons())	ASSERT(0);

		// Set Size.
		iResPolySize = PolyMakerUnit.GetResultPolyCount();

		// For Debugging.
		iDebugCount=0;		
		// Result Polygon.		
		for(int j=0; j<iResPolySize; j++)
		{				
			CArray<DGL_3dp,DGL_3dp&> arResPoly;
			PolyMakerUnit.GetResultPoly(j,arResPoly);
			int iResNum = arResPoly.GetSize();
			// Calculate Center(x,y), Area.
			//---------------------------------
			double* dpx = new double[iResNum];
			double* dpy = new double[iResNum];
			//---------------------------------
			bInChk = TRUE;
			for(int k=0; k<iResNum; k++)
			{
				dpx[k] = arResPoly[k].x();
				dpy[k] = arResPoly[k].y();
				if(dpy[k] > dZ_Top+cUMDRC_Zero || dpy[k] < dZ_Bot-cUMDRC_Zero)
					bInChk = FALSE;
			}
			if(!bInChk) continue;

			double dxCen=0.0, dyCen=0.0, dArea=0.0;
			if(!CMathFunc::mathPolyCentroid(iResNum,dpx,dpy,dxCen,dyCen,dArea))
			{// SHIN ('06.02.26) 모서리에 걸릴경우 3개의 좌표값이 동일한 경우가 생기는데 이때 Debug시 ASSERT(0)의 호출을 피하기 위해
				BOOL bXCheck, bYCheck;
				bXCheck = bYCheck = TRUE;
				for(int k=0; k<iResNum-1; k++)
				{
					if(fabs(dpx[k]-dpx[k+1]) > cUMDRC_Zero)
					{
						bXCheck = FALSE; break;
					}
					if(fabs(dpy[k]-dpy[k+1]) > cUMDRC_Zero)
					{
						bYCheck = FALSE; break;
					}
				}
				if((fabs(dArea) > cUMDRC_Zero) || bXCheck == FALSE || bYCheck == FALSE)  ASSERT(0);
			}
			//---------------------------------
			delete []dpx;
			delete []dpy;
			//---------------------------------
			
			center_y -= dArea*dxCen;
			center_z -= dArea*dyCen;
			dArea_Inner -= dArea;   
		}
		PolyMakerUnit.RemoveAllCutPolys();
		PolyMakerUnit.RemoveAllResult();
	}


	if(dArea_Inner > 0.0)
	{
		center_y = center_y/dArea_Inner;
		center_z = center_z/dArea_Inner;
	}
	else
	{
		center_y = center_z = dArea_Inner = 0.0;
	}
	
	return TRUE;
}

BOOL CMdgnSectTool::MakeConDataReg(int iShape, double dSize[8], _UMD_RC_CONC& UmdRcConc, int iDivNum)
{	
	UmdRcConc.Init();

	double dPi = 4.*atan(1.0);
  double dDivNum = (double)iDivNum;
  
  _UMD_RC_CONC_UNIT ConcUnit;
	BOOL bMeshOK=TRUE;
  
  if(iShape==DGN_SECT_SHAPE_INDEX_REG_T)  //T
  {
		double dHc = dSize[0];
		double dBc = dSize[1];
		double dTw = dSize[2];
		double dTf = dSize[3];		

    int iConcKey = 1;
    double dArea1=0.0, duv1[2]={0.0, 0.0}, duvPnt1[4][2]={0.0};
		// Flange.
		double dxx1[4] = {0.0, dBc, dBc, 0.0};
		double dyy1[4] = {dHc-dTf, dHc-dTf, dHc, dHc};
		// Get area, centroid, points.
		CMathFunc::mathPolyCentroid(4,dxx1,dyy1,duv1[0],duv1[1],dArea1);
    dArea1 = fabs(dArea1);
		duvPnt1[0][0] = dxx1[0]; duvPnt1[0][1] = dyy1[0];
		duvPnt1[1][0] = dxx1[1]; duvPnt1[1][1] = dyy1[1];
		duvPnt1[2][0] = dxx1[2]; duvPnt1[2][1] = dyy1[2];
		duvPnt1[3][0] = dxx1[3]; duvPnt1[3][1] = dyy1[3];
		ConcUnit.Init();
		ConcUnit.dArea	= dArea1;
		ConcUnit.dyz[0] = duv1[0];	// y.
		ConcUnit.dyz[1] = duv1[1];	// z.
		for(int k=0; k<4; k++)
		{
			ConcUnit.dyzPnt[k][0] = duvPnt1[k][0];
			ConcUnit.dyzPnt[k][1] = duvPnt1[k][1];
		}
		UmdRcConc.arConcUnit.SetAt(iConcKey, ConcUnit);

    // Web.
    iConcKey = 2;
    double dArea2=0.0, duv2[2]={0.0, 0.0}, duvPnt2[4][2]={0.0};
		double dxx2[4] = {(dBc-dTw)/2., (dBc+dTw)/2., (dBc+dTw)/2., (dBc-dTw)/2.};
		double dyy2[4] = {0.0, 0.0, dHc-dTf, dHc-dTf};
		// Get area, centroid, points.
		CMathFunc::mathPolyCentroid(4,dxx2,dyy2,duv2[0],duv2[1],dArea2);
    dArea2 = fabs(dArea2);
		duvPnt2[0][0] = dxx2[0]; duvPnt2[0][1] = dyy2[0];
		duvPnt2[1][0] = dxx2[1]; duvPnt2[1][1] = dyy2[1];
		duvPnt2[2][0] = dxx2[2]; duvPnt2[2][1] = dyy2[2];
		duvPnt2[3][0] = dxx2[3]; duvPnt2[3][1] = dyy2[3];

		ConcUnit.Init();
		ConcUnit.dArea	= dArea2;
		ConcUnit.dyz[0] = duv2[0];	// y.
		ConcUnit.dyz[1] = duv2[1];	// z.
		for(int k=0; k<4; k++)
		{
			ConcUnit.dyzPnt[k][0] = duvPnt2[k][0];
			ConcUnit.dyzPnt[k][1] = duvPnt2[k][1];
		}
		UmdRcConc.arConcUnit.SetAt(iConcKey, ConcUnit);
  }
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_B)	// Box.
	{
		double dH  = dSize[0];
		double dB  = dSize[1];
		double dTw = dSize[2];
		double dTf = dSize[3];

		double dxx[4][4];
		double dyy[4][4];
    for(int i=0; i<4; i++)
    {
      for(int j=0; j<4; j++)
      {
        dxx[i][j] = 0.0; dyy[i][j] = 0.0;
      }
    }
    //Bottom Rect
		dxx[0][0] =    0.0; dxx[0][1] =     dB; dxx[0][2] =     dB; dxx[0][3] =    0.0;
		dyy[0][0] =    0.0; dyy[0][1] =    0.0; dyy[0][2] =    dTf; dyy[0][3] =    dTf;
    //Top Rect
		dxx[1][0] =    0.0; dxx[1][1] =     dB; dxx[1][2] =     dB; dxx[1][3] =    0.0;
		dyy[1][0] = dH-dTf; dyy[1][1] = dH-dTf; dyy[1][2] =     dH; dyy[1][3] =     dH;
    //Left Rect
		dxx[2][0] =    0.0; dxx[2][1] =    dTw; dxx[2][2] =    dTw; dxx[2][3] =    0.0;
		dyy[2][0] =    dTf; dyy[2][1] =    dTf; dyy[2][2] = dH-dTf; dyy[2][3] = dH-dTf;
    //Right Rect
		dxx[3][0] = dB-dTw; dxx[3][1] =     dB; dxx[3][2] =     dB; dxx[3][3] = dB-dTw;
		dyy[3][0] =    dTf; dyy[3][1] =    dTf; dyy[3][2] = dH-dTf; dyy[3][3] = dH-dTf;

		/////////////////////////////////////////////
		int iConcKey=0, iCountElem=0;
		double dTArea = 0.0;
    for(int i=0; i<4; i++)
    {
		  double dArea=0.0, duv[2]={0.0, 0.0}, duvPnt[4][2]={0.0};
		  // Bttm Rectangle.
		  // Get area, centroid, points.
		  CMathFunc::mathPolyCentroid(4,dxx[i],dyy[i],duv[0],duv[1],dArea);
      dArea = fabs(dArea);
      dTArea += dArea;
		  duvPnt[0][0] = dxx[i][0]; duvPnt[0][1] = dyy[i][0];
		  duvPnt[1][0] = dxx[i][1]; duvPnt[1][1] = dyy[i][1];
		  duvPnt[2][0] = dxx[i][2]; duvPnt[2][1] = dyy[i][2];
		  duvPnt[3][0] = dxx[i][3]; duvPnt[3][1] = dyy[i][3];

      iCountElem++;
		  ConcUnit.Init();
		  iConcKey = iCountElem;
		  ConcUnit.dArea	= dArea;
		  ConcUnit.dyz[0] = duv[0];	// y.
		  ConcUnit.dyz[1] = duv[1];	// z.
		  for(int k=0; k<4; k++)
		  {
			  ConcUnit.dyzPnt[k][0] = duvPnt[k][0];
			  ConcUnit.dyzPnt[k][1] = duvPnt[k][1];
		  }
		  UmdRcConc.arConcUnit.SetAt(iConcKey, ConcUnit);
    }
  }
  else if(iShape==DGN_SECT_SHAPE_INDEX_REG_P)	// P.
  {
    double dR = dSize[0]/2.;
    double dr = dR - dSize[1];

		int iCountElem=0;		
    double dPreAngle = 0.0, dNextAngle = 0.0;
    double dTArea = 0.0;
	  for(int i=0; i<dDivNum; i++)
		{
			double dArea=0.0, duv[2]={0.0, 0.0}, duvPnt[4][2]={0.0};
      dPreAngle = 2.*dPi*i/dDivNum;
      dNextAngle = 2.*dPi*(i+1)/dDivNum;
			// Get area, centroid, points.
			double dxx[4] = {dR+dr*cos(dPreAngle), dR+dr*cos(dNextAngle), dR+dR*cos(dNextAngle), dR+dR*cos(dPreAngle)};
			double dyy[4] = {dR-dr*sin(dPreAngle), dR-dr*sin(dNextAngle), dR-dR*sin(dNextAngle), dR-dR*sin(dPreAngle)};
			CMathFunc::mathPolyCentroid(4,dxx,dyy,duv[0],duv[1],dArea);
      dArea = fabs(dArea);
			duvPnt[0][0] = dxx[0]; duvPnt[0][1] = dyy[0];
			duvPnt[1][0] = dxx[1]; duvPnt[1][1] = dyy[1];
			duvPnt[2][0] = dxx[2]; duvPnt[2][1] = dyy[2];
			duvPnt[3][0] = dxx[3]; duvPnt[3][1] = dyy[3];

      iCountElem++;
			ConcUnit.Init();
			int iConcKey = iCountElem;
			ConcUnit.dArea	= dArea;
			ConcUnit.dyz[0] = duv[0];	// y.
			ConcUnit.dyz[1] = duv[1];	// z.
			for(int k=0; k<4; k++)
			{
				ConcUnit.dyzPnt[k][0] = duvPnt[k][0];
				ConcUnit.dyzPnt[k][1] = duvPnt[k][1];
			}
			UmdRcConc.arConcUnit.SetAt(iConcKey, ConcUnit);
    }
  }
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_SR)	// SR.
	{
    double dR = dSize[0]/2.;

		int iCountElem=0;		
    double dPreAngle = 0.0, dNextAngle = 0.0;
    double dTArea = 0.0;
	  for(int i=0; i<dDivNum; i++)
		{
			double dArea=0.0, duv[2]={0.0, 0.0}, duvPnt[4][2]={0.0};
      dPreAngle = 2.*dPi*i/dDivNum;
      dNextAngle = 2.*dPi*(i+1)/dDivNum;
			// Triangle.
			// Get area, centroid, points.
  		double dxx[4] = {dR+cUMDRC_Zero*cos(dPreAngle), dR+dR*cos(dPreAngle), dR+dR*cos(dNextAngle), dR+cUMDRC_Zero*cos(dNextAngle)};
			double dyy[4] = {dR-cUMDRC_Zero*sin(dPreAngle), dR-dR*sin(dPreAngle), dR-dR*sin(dNextAngle), dR-cUMDRC_Zero*sin(dNextAngle)};
			CMathFunc::mathPolyCentroid(4,dxx,dyy,duv[0],duv[1],dArea);
      dArea = fabs(dArea);
      dTArea += dArea;
			duvPnt[0][0] = dxx[0]; duvPnt[0][1] = dyy[0];
			duvPnt[1][0] = dxx[1]; duvPnt[1][1] = dyy[1];
			duvPnt[2][0] = dxx[2]; duvPnt[2][1] = dyy[2];
			duvPnt[3][0] = dxx[3]; duvPnt[3][1] = dyy[3];

      iCountElem++;
			ConcUnit.Init();
			int iConcKey = iCountElem;
			ConcUnit.dArea	= dArea;
			ConcUnit.dyz[0] = duv[0];	// y.
			ConcUnit.dyz[1] = duv[1];	// z.
			for(int j=0; j<4; j++)
			{
				ConcUnit.dyzPnt[j][0] = duvPnt[j][0];
				ConcUnit.dyzPnt[j][1] = duvPnt[j][1];
			}
			UmdRcConc.arConcUnit.SetAt(iConcKey, ConcUnit);
		}

	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_SB)	// SB.
	{
    double dDivBNum = 1;
    double dDivHNum = 1;
		double dH = dSize[0];
		double dB = dSize[1];

		/////////////////////////////////////////////
		int iCountElem=0;		
	  for(int i=0; i<dDivHNum; i++)
		{
      for(int j=0; j<dDivBNum; j++)
      {
			  double dArea=0.0, duv[2]={0.0, 0.0}, duvPnt[4][2]={0.0};
			  // Rectangle.
			  double dxx[4] = {dB*j/dDivBNum, dB*(j+1)/dDivBNum, dB*(j+1)/dDivBNum, dB*j/dDivBNum};
			  double dyy[4] = {dH*i/dDivHNum,     dH*i/dDivHNum, dH*(i+1)/dDivHNum, dH*(i+1)/dDivHNum};
			  // Get area, centroid, points.
			  CMathFunc::mathPolyCentroid(4,dxx,dyy,duv[0],duv[1],dArea);
        dArea = fabs(dArea);
			  duvPnt[0][0] = dxx[0]; duvPnt[0][1] = dyy[0];
			  duvPnt[1][0] = dxx[1]; duvPnt[1][1] = dyy[1];
			  duvPnt[2][0] = dxx[2]; duvPnt[2][1] = dyy[2];
			  duvPnt[3][0] = dxx[3]; duvPnt[3][1] = dyy[3];

        iCountElem++;
			  ConcUnit.Init();
			  int iConcKey = iCountElem;
			  ConcUnit.dArea	= dArea;
			  ConcUnit.dyz[0] = duv[0];	// y.
			  ConcUnit.dyz[1] = duv[1];	// z.
			  for(int k=0; k<4; k++)
			  {
				  ConcUnit.dyzPnt[k][0] = duvPnt[k][0];
				  ConcUnit.dyzPnt[k][1] = duvPnt[k][1];
			  }
			  UmdRcConc.arConcUnit.SetAt(iConcKey, ConcUnit);
      }
		}
  }
  else if(iShape==DGN_SECT_SHAPE_INDEX_REG_OCT) //Oct
  {
		double dH   = dSize[0];
		double dB   = dSize[1];
		double dDa  = dSize[2];
    double dDb  = dSize[3];
    double dThk = dSize[4];

    double dH1=0.0, dB1=0.0, dDa1=0.0, dDb1=0.0, dDist0=0.0, dTheta1=0.0;    
    CalcSectRegular_InOCT(dSize,dH1,dB1,dDa1,dDb1,dDist0,dTheta1);

		double dxx[8][4];
		double dyy[8][4];
    for(int i=0; i<8; i++)
    {
      for(int j=0; j<4; j++)
      {
        dxx[i][j] = 0.0; dyy[i][j] = 0.0;
      }
    }
    //Rect
		dxx[0][0] =     dDa; dxx[0][1] =  dB-dDa; dxx[0][2] =  dB1-dDa1+dThk; dxx[0][3] =      dDa1+dThk;
		dyy[0][0] =     0.0; dyy[0][1] =     0.0; dyy[0][2] =           dThk; dyy[0][3] =           dThk;

    dxx[1][0] =  dB-dDa; dxx[1][1] =      dB; dxx[1][2] =       dB1+dThk; dxx[1][3] =  dB1-dDa1+dThk;
		dyy[1][0] =     0.0; dyy[1][1] =     dDb; dyy[1][2] =      dDb1+dThk; dyy[1][3] =           dThk;

    dxx[2][0] =      dB; dxx[2][1] =      dB; dxx[2][2] =       dB1+dThk; dxx[2][3] =       dB1+dThk;
		dyy[2][0] =     dDb; dyy[2][1] =  dH-dDb; dyy[2][2] =  dH1-dDb1+dThk; dyy[2][3] =      dDb1+dThk;

    dxx[3][0] =      dB; dxx[3][1] =  dB-dDa; dxx[3][2] =  dB1-dDa1+dThk; dxx[3][3] =       dB1+dThk;
		dyy[3][0] =  dH-dDb; dyy[3][1] =      dH; dyy[3][2] =       dH1+dThk; dyy[3][3] =  dH1-dDb1+dThk;

		dxx[4][0] =  dB-dDa; dxx[4][1] =     dDa; dxx[4][2] =      dDa1+dThk; dxx[4][3] =  dB1-dDa1+dThk;
		dyy[4][0] =      dH; dyy[4][1] =      dH; dyy[4][2] =       dH1+dThk; dyy[4][3] =       dH1+dThk;

    dxx[5][0] =     dDa; dxx[5][1] =     0.0; dxx[5][2] =           dThk; dxx[5][3] =      dDa1+dThk;
		dyy[5][0] =      dH; dyy[5][1] =  dH-dDb; dyy[5][2] =  dH1-dDb1+dThk; dyy[5][3] =       dH1+dThk;

    dxx[6][0] =     0.0; dxx[6][1] =     0.0; dxx[6][2] =           dThk; dxx[6][3] =           dThk;
		dyy[6][0] =  dH-dDb; dyy[6][1] =     dDb; dyy[6][2] =      dDb1+dThk; dyy[6][3] =  dH1-dDb1+dThk;

    dxx[7][0] =     0.0; dxx[7][1] =     dDb; dxx[7][2] =      dDb1+dThk; dxx[7][3] =           dThk;
		dyy[7][0] =     dDb; dyy[7][1] =     0.0; dyy[7][2] =           dThk; dyy[7][3] =      dDb1+dThk;

    /////////////////////////////////////////////
		int iConcKey=0, iCountElem=0;
		double dTArea = 0.0;
    for(int i=0; i<8; i++)
    {
		  double dArea=0.0, duv[2]={0.0, 0.0}, duvPnt[4][2]={0.0};
		  // Bttm Rectangle.
		  // Get area, centroid, points.
		  CMathFunc::mathPolyCentroid(4,dxx[i],dyy[i],duv[0],duv[1],dArea);
      dArea = fabs(dArea);
      dTArea += dArea;
		  duvPnt[0][0] = dxx[i][0]; duvPnt[0][1] = dyy[i][0];
		  duvPnt[1][0] = dxx[i][1]; duvPnt[1][1] = dyy[i][1];
		  duvPnt[2][0] = dxx[i][2]; duvPnt[2][1] = dyy[i][2];
		  duvPnt[3][0] = dxx[i][3]; duvPnt[3][1] = dyy[i][3];

      iCountElem++;
		  ConcUnit.Init();
		  iConcKey = iCountElem;
		  ConcUnit.dArea	= dArea;
		  ConcUnit.dyz[0] = duv[0];	// y.
		  ConcUnit.dyz[1] = duv[1];	// z.
		  for(int k=0; k<4; k++)
		  {
			  ConcUnit.dyzPnt[k][0] = duvPnt[k][0];
			  ConcUnit.dyzPnt[k][1] = duvPnt[k][1];
		  }
		  UmdRcConc.arConcUnit.SetAt(iConcKey, ConcUnit);
    }
  }
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_SOCT)	// SOct.
	{
		double dH   = dSize[0];
		double dB   = dSize[1];
		double dDa  = dSize[2];
    double dDb  = dSize[3];

		double dxx[3][4];
		double dyy[3][4];
    for(int i=0; i<3; i++)
    {
      for(int j=0; j<4; j++)
      {
        dxx[i][j] = 0.0; dyy[i][j] = 0.0;
      }
    }

    //Rect
		dxx[0][0] =    dDa; dxx[0][1] = dB-dDa; dxx[0][2] = dB-dDa; dxx[0][3] =    dDa;
		dyy[0][0] =    0.0; dyy[0][1] =    0.0; dyy[0][2] =     dH; dyy[0][3] =     dH;
    //Left Rect
		dxx[1][0] =    0.0; dxx[1][1] =    dDa; dxx[1][2] =    dDa; dxx[1][3] =    0.0;
		dyy[1][0] =    dDb; dyy[1][1] =    0.0; dyy[1][2] =     dH; dyy[1][3] = dH-dDb;
    //Right Rect
		dxx[2][0] = dB-dDa; dxx[2][1] =     dB; dxx[2][2] =     dB; dxx[2][3] = dB-dDa;
		dyy[2][0] =    0.0; dyy[2][1] =    dDb; dyy[2][2] = dH-dDb; dyy[2][3] =     dH;

		/////////////////////////////////////////////
		int iConcKey=0, iCountElem=0;
		double dTArea = 0.0;
    for(int i=0; i<3; i++)
    {
		  double dArea=0.0, duv[2]={0.0, 0.0}, duvPnt[4][2]={0.0};
		  // Bttm Rectangle.
		  // Get area, centroid, points.
		  CMathFunc::mathPolyCentroid(4,dxx[i],dyy[i],duv[0],duv[1],dArea);
      dArea = fabs(dArea);
      dTArea += dArea;
		  duvPnt[0][0] = dxx[i][0]; duvPnt[0][1] = dyy[i][0];
		  duvPnt[1][0] = dxx[i][1]; duvPnt[1][1] = dyy[i][1];
		  duvPnt[2][0] = dxx[i][2]; duvPnt[2][1] = dyy[i][2];
		  duvPnt[3][0] = dxx[i][3]; duvPnt[3][1] = dyy[i][3];

      iCountElem++;
		  ConcUnit.Init();
		  iConcKey = iCountElem;
		  ConcUnit.dArea	= dArea;
		  ConcUnit.dyz[0] = duv[0];	// y.
		  ConcUnit.dyz[1] = duv[1];	// z.
		  for(int k=0; k<4; k++)
		  {
			  ConcUnit.dyzPnt[k][0] = duvPnt[k][0];
			  ConcUnit.dyzPnt[k][1] = duvPnt[k][1];
		  }
		  UmdRcConc.arConcUnit.SetAt(iConcKey, ConcUnit);
    }
  }
  else if(iShape==DGN_SECT_SHAPE_INDEX_REG_TRK)	// Trk.
  {
    double dH  = dSize[0];
    double dB  = dSize[1];
    double dThk= dSize[2];
    double dR = dH/2.;
    double dr = dR - dThk;
    double dIncDeg = 2.*dPi/dDivNum;

		int iCountElem=0;
		double dPreAngle = 0.0, dNextAngle = 0.0;
    double dTArea = 0.0;
	  for(int i=0; i<dDivNum; i++)
		{
			double dArea=0.0, duv[2]={0.0, 0.0}, duvPnt[4][2]={0.0};
      dPreAngle = 2.*dPi*(i-dDivNum/4.)/dDivNum;
      dNextAngle = dPreAngle + dIncDeg;
			// Get area, centroid, points.
      double dxx[4], dyy[4];
      if(i >= dDivNum/2.)
      {
			  dxx[0] = dR+dr*cos(dPreAngle); dxx[1] = dR+dr*cos(dNextAngle); dxx[2] = dR+dR*cos(dNextAngle); dxx[3] = dR+dR*cos(dPreAngle);
			  dyy[0] = dR-dr*sin(dPreAngle); dyy[1] = dR-dr*sin(dNextAngle); dyy[2] = dR-dR*sin(dNextAngle); dyy[3] = dR-dR*sin(dPreAngle);
      }
      else
      {
			  dxx[0] = dB-dR+dr*cos(dPreAngle); dxx[1] = dB-dR+dr*cos(dNextAngle); dxx[2] = dB-dR+dR*cos(dNextAngle); dxx[3] = dB-dR+dR*cos(dPreAngle);
			  dyy[0] = dR-dr*sin(dPreAngle);    dyy[1] = dR-dr*sin(dNextAngle);    dyy[2] = dR-dR*sin(dNextAngle);    dyy[3] = dR-dR*sin(dPreAngle);
      }

			CMathFunc::mathPolyCentroid(4,dxx,dyy,duv[0],duv[1],dArea);
      dArea = fabs(dArea);
      dTArea += dArea;
			duvPnt[0][0] = dxx[0]; duvPnt[0][1] = dyy[0];
			duvPnt[1][0] = dxx[1]; duvPnt[1][1] = dyy[1];
			duvPnt[2][0] = dxx[2]; duvPnt[2][1] = dyy[2];
			duvPnt[3][0] = dxx[3]; duvPnt[3][1] = dyy[3];

      iCountElem++;
			ConcUnit.Init();
			int iConcKey = iCountElem;
			ConcUnit.dArea	= dArea;
			ConcUnit.dyz[0] = duv[0];	// y.
			ConcUnit.dyz[1] = duv[1];	// z.
			for(int k=0; k<4; k++)
			{
				ConcUnit.dyzPnt[k][0] = duvPnt[k][0];
				ConcUnit.dyzPnt[k][1] = duvPnt[k][1];
			}
			UmdRcConc.arConcUnit.SetAt(iConcKey, ConcUnit);
    }
		// Rectangle.
    double dxx[2][4], dyy[2][4];
		dxx[0][0] = dR     ; dxx[0][1] = dB-dR  ; dxx[0][2] = dB-dR; dxx[0][3] =   dR;
		dyy[0][0] = 0.0    ; dyy[0][1] =   0.0  ; dyy[0][2] =  dThk; dyy[0][3] = dThk;
		dxx[1][0] = dR     ; dxx[1][1] = dB-dR  ; dxx[1][2] = dB-dR; dxx[1][3] =   dR;
		dyy[1][0] = dH-dThk; dyy[1][1] = dH-dThk; dyy[1][2] =  dH  ; dyy[1][3] =   dH;
    for(int i=0; i<2; i++)
    {
      double dArea=0.0, duv[2]={0.0, 0.0}, duvPnt[4][2]={0.0};
		  // Get area, centroid, points.
		  CMathFunc::mathPolyCentroid(4,dxx[i],dyy[i],duv[0],duv[1],dArea);
      dArea = fabs(dArea);
      dTArea += dArea;
		  duvPnt[0][0] = dxx[i][0]; duvPnt[0][1] = dyy[i][0];
		  duvPnt[1][0] = dxx[i][1]; duvPnt[1][1] = dyy[i][1];
		  duvPnt[2][0] = dxx[i][2]; duvPnt[2][1] = dyy[i][2];
		  duvPnt[3][0] = dxx[i][3]; duvPnt[3][1] = dyy[i][3];

      iCountElem++;
		  ConcUnit.Init();
		  int iConcKey = iCountElem;
		  ConcUnit.dArea	= dArea;
		  ConcUnit.dyz[0] = duv[0];	// y.
		  ConcUnit.dyz[1] = duv[1];	// z.
		  for(int k=0; k<4; k++)
		  {
			  ConcUnit.dyzPnt[k][0] = duvPnt[k][0];
			  ConcUnit.dyzPnt[k][1] = duvPnt[k][1];
		  }
		  UmdRcConc.arConcUnit.SetAt(iConcKey, ConcUnit);
    }
  }
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_STRK)	// STRK.
	{
    double dH = dSize[0];
    double dB = dSize[1];
    double dR = dH/2.;
    double dTArea = 0.0;
		int iCountElem=0;
		double dPreAngle = 0.0, dNextAngle = 0.0;
    double dIncDeg = 2.*dPi/dDivNum;
    //Circle Part
	  for(int i=0; i<dDivNum; i++)
		{
			double dArea=0.0, duv[2]={0.0, 0.0}, duvPnt[4][2]={0.0};
      dPreAngle = 2.*dPi*(i-dDivNum/4.)/dDivNum;
      dNextAngle = dPreAngle + dIncDeg;
			// Triangle.
			// Get area, centroid, points.
      double dxx[4], dyy[4];
      if(i >= dDivNum/2.)
      {
			  dxx[0] = dR+cUMDRC_Zero*cos(dPreAngle); dxx[1] = dR+dR*cos(dPreAngle); dxx[2] = dR+dR*cos(dNextAngle); dxx[3] = dR+cUMDRC_Zero*cos(dNextAngle);
			  dyy[0] = dR-cUMDRC_Zero*sin(dPreAngle); dyy[1] = dR-dR*sin(dPreAngle); dyy[2] = dR-dR*sin(dNextAngle); dyy[3] = dR-cUMDRC_Zero*sin(dNextAngle);
      }
      else
      {
			  dxx[0] = dB-dR+cUMDRC_Zero*cos(dPreAngle); dxx[1] = dB-dR+dR*cos(dPreAngle); dxx[2] = dB-dR+dR*cos(dNextAngle);  dxx[3] = dB-dR+cUMDRC_Zero*cos(dNextAngle);
			  dyy[0] = dR-cUMDRC_Zero*sin(dPreAngle);    dyy[1] = dR-dR*sin(dPreAngle);    dyy[2] = dR-dR*sin(dNextAngle);     dyy[3] = dR-cUMDRC_Zero*sin(dNextAngle);
      }
			CMathFunc::mathPolyCentroid(4,dxx,dyy,duv[0],duv[1],dArea);
      dArea = fabs(dArea);
      dTArea += dArea;
			duvPnt[0][0] = dxx[0]; duvPnt[0][1] = dyy[0];
			duvPnt[1][0] = dxx[1]; duvPnt[1][1] = dyy[1];
			duvPnt[2][0] = dxx[2]; duvPnt[2][1] = dyy[2];
			duvPnt[3][0] = dxx[3]; duvPnt[3][1] = dyy[3];

      iCountElem++;
			ConcUnit.Init();
			int iConcKey = iCountElem;
			ConcUnit.dArea	= dArea;
			ConcUnit.dyz[0] = duv[0];	// y.
			ConcUnit.dyz[1] = duv[1];	// z.
			for(int j=0; j<4; j++)
			{
				ConcUnit.dyzPnt[j][0] = duvPnt[j][0];
				ConcUnit.dyzPnt[j][1] = duvPnt[j][1];
			}
			UmdRcConc.arConcUnit.SetAt(iConcKey, ConcUnit);
		}

		// Rectangle.
		double dxx[4] = {dR , dB-dR, dB-dR, dR};
		double dyy[4] = {0.0,   0.0,    dH, dH};
    double dArea=0.0, duv[2]={0.0, 0.0}, duvPnt[4][2]={0.0};
		// Get area, centroid, points.
		CMathFunc::mathPolyCentroid(4,dxx,dyy,duv[0],duv[1],dArea);
    dArea = fabs(dArea);
    dTArea += dArea;
		duvPnt[0][0] = dxx[0]; duvPnt[0][1] = dyy[0];
		duvPnt[1][0] = dxx[1]; duvPnt[1][1] = dyy[1];
		duvPnt[2][0] = dxx[2]; duvPnt[2][1] = dyy[2];
		duvPnt[3][0] = dxx[3]; duvPnt[3][1] = dyy[3];

    iCountElem++;
		ConcUnit.Init();
		int iConcKey = iCountElem;
		ConcUnit.dArea	= dArea;
		ConcUnit.dyz[0] = duv[0];	// y.
		ConcUnit.dyz[1] = duv[1];	// z.
		for(int k=0; k<4; k++)
		{
			ConcUnit.dyzPnt[k][0] = duvPnt[k][0];
			ConcUnit.dyzPnt[k][1] = duvPnt[k][1];
		}
		UmdRcConc.arConcUnit.SetAt(iConcKey, ConcUnit);    
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_HTRK)	// HTRK.
	{
    dDivNum = dDivNum/2.;
    double dH = dSize[0];
    double dB = dSize[1];
    double dR = dH/2.;
    double dTArea = 0.0;
		int iCountElem=0;		
    double dPreAngle = 0.0, dNextAngle = 0.0;
    double dIncDeg = dPi/dDivNum;
    //Circle Part
	  for(int i=0; i<dDivNum; i++)
		{
			double dArea=0.0, duv[2]={0.0, 0.0}, duvPnt[4][2]={0.0};
      dPreAngle = dPi*(i-dDivNum/2.)/dDivNum;
      dNextAngle = dPreAngle + dIncDeg;
			// Triangle.
			// Get area, centroid, points.
      double dxx[4], dyy[4];
			dxx[0] = dB-dR+cUMDRC_Zero*cos(dPreAngle); dxx[1] = dB-dR+dR*cos(dPreAngle); dxx[2] = dB-dR+dR*cos(dNextAngle);  dxx[3] = dB-dR+cUMDRC_Zero*cos(dNextAngle);
			dyy[0] = dR-cUMDRC_Zero*sin(dPreAngle);    dyy[1] = dR-dR*sin(dPreAngle);    dyy[2] = dR-dR*sin(dNextAngle);     dyy[3] = dR-cUMDRC_Zero*sin(dNextAngle);
			CMathFunc::mathPolyCentroid(4,dxx,dyy,duv[0],duv[1],dArea);
      dArea = fabs(dArea);
      dTArea += dArea;
			duvPnt[0][0] = dxx[0]; duvPnt[0][1] = dyy[0];
			duvPnt[1][0] = dxx[1]; duvPnt[1][1] = dyy[1];
			duvPnt[2][0] = dxx[2]; duvPnt[2][1] = dyy[2];
			duvPnt[3][0] = dxx[3]; duvPnt[3][1] = dyy[3];

      iCountElem++;
			ConcUnit.Init();
			int iConcKey = iCountElem;
			ConcUnit.dArea	= dArea;
			ConcUnit.dyz[0] = duv[0];	// y.
			ConcUnit.dyz[1] = duv[1];	// z.
			for(int j=0; j<4; j++)
			{
				ConcUnit.dyzPnt[j][0] = duvPnt[j][0];
				ConcUnit.dyzPnt[j][1] = duvPnt[j][1];
			}
			UmdRcConc.arConcUnit.SetAt(iConcKey, ConcUnit);
		}

		// Rectangle.
		double dxx[4] = {0.0, dB-dR, dB-dR, 0.0};
		double dyy[4] = {0.0,   0.0,    dH, dH};
    double dArea=0.0, duv[2]={0.0, 0.0}, duvPnt[4][2]={0.0};
		// Get area, centroid, points.
		CMathFunc::mathPolyCentroid(4,dxx,dyy,duv[0],duv[1],dArea);
    dArea = fabs(dArea);
    dTArea += dArea;
		duvPnt[0][0] = dxx[0]; duvPnt[0][1] = dyy[0];
		duvPnt[1][0] = dxx[1]; duvPnt[1][1] = dyy[1];
		duvPnt[2][0] = dxx[2]; duvPnt[2][1] = dyy[2];
		duvPnt[3][0] = dxx[3]; duvPnt[3][1] = dyy[3];

    iCountElem++;
		ConcUnit.Init();
		int iConcKey = iCountElem;
		ConcUnit.dArea	= dArea;
		ConcUnit.dyz[0] = duv[0];	// y.
		ConcUnit.dyz[1] = duv[1];	// z.
		for(int k=0; k<4; k++)
		{
			ConcUnit.dyzPnt[k][0] = duvPnt[k][0];
			ConcUnit.dyzPnt[k][1] = duvPnt[k][1];
		}
		UmdRcConc.arConcUnit.SetAt(iConcKey, ConcUnit);
	}
  else 
	{
		ASSERT(0);
		bMeshOK = FALSE;
	}

	return bMeshOK;
}

BOOL CMdgnSectTool::CalcSectRegular_InOCT(double dSize[8], 
																												double& dH1, double& dB1, double& da1, double& db1, double& dDist0, double& dTheta1)
{
	double dPi = 4.*atan(1.0);
	
	double dH	= dSize[0];
	double dB	= dSize[1];
	double da	= dSize[2];
	double db	= dSize[3];
	double dt	= dSize[4];
	
	
	if(da == 0.0 || db == 0.0) 
	{
		dH1 =  max(0.0, dH - 2.0*dt);
		dB1 =  max(0.0, dB - 2.0*dt);
		da1 = db1 = 0.0;
		dTheta1 = dPi/2.;
		dDist0  = sqrt(2.0*dt*dt);
	}
	else
	{
		double dLineI[3]={-da,0.0,0.0};
		double dLineJ[3]={0.0, db,0.0};
		double dPoint[3]={-dt, dt,0.0};
		double dDist	= CMathFunc::mathDistanceFromIntersectPointToLine(dLineI,dLineJ,dPoint);
		double dPontO[3]={0.0,0.0,0.0};
		dDist0	= CMathFunc::mathDistanceFromIntersectPointToLine(dLineI,dLineJ,dPontO);
		dTheta1 = dPi/2. - atan(db/da);
		
		double dDist1	 = dDist0 + dt;
		double dx1 = (-1.)*dDist1*cos(dTheta1);
		double dy1 =       dDist1*sin(dTheta1);
		double dc1 = dy1 - (db/da)*dx1;
		
		double dPosiY = (-1.)*dt*(db/da) + db;
		dH1 = dH - 2.*dt;
		dB1 = dB - 2.*dt;
		da1 = 0.0;
		db1 = 0.0;
		if(dPosiY <= dt && dDist >= dt)	return FALSE;
		else
		{
			// Formula by ZINU. When dCompY1,dCompY2 = dt. 
			double dCompX1 =       (dt-dc1)*(da/db);
			double dCompX2 = (-1.)*(dt-dc1)*(da/db) - dB;
			if(dCompX1 <= dCompX2)	return FALSE;
			// Formula by ZINU. When dCompX1,dCompX2 = -dt. 
			double dCompY1 =       (db/da)*(-dt) + dc1;
			double dCompY2 = (-1.)*(db/da)*(-dt) + (dH-dc1);
			if(dCompY1 >= dCompY2)	return FALSE;
			// Calculate a1, b1.
			da1 = (-1.)*(dt-dc1)*(da/db) - dt;
			db1 = (-1.)*dt*(db/da) + dc1 - dt;
			if(da1 <= 0.0 || db1 <= 0.0)	return FALSE;
		}
	}
	return TRUE;
}

BOOL CMdgnSectTool::MakeRbarDataReg(int iShape, double dSize[8], const _UMD_RC_COL_MAINRBAR &RbarData, _UMD_RC_RBAR& UmdRcRbar)
{
	UmdRcRbar.Init();
	double dPi = 4.*atan(1.0);
	// To calculate Covering.
	// Only Design.
	BOOL   bDesign = FALSE;
	double dDcDgn = 0.0; ///!!!
	double dUserDcDgn=Get_CvrThkOnColmDgn(iShape, dSize, dDcDgn);		
	
	int iGrup;
	
	if(iShape==DGN_SECT_SHAPE_INDEX_REG_B)	// B.
	{
		// Column, Brace.
		double dHc=0., dB1=0., dTw=0., dTf1=0., dB2=0., dTf2=0.;
      
		dHc = dSize[0];
		dB1 = dSize[1];
		dTw = dSize[2];
		dTf1= dSize[3];
		dB2 = dSize[4];
		dTf2= dSize[5];

		if(bDesign)	// Design.
		{
			// Covering.
			double dDc = dUserDcDgn;
      
			// Top, Left, Bot, Right (StaY, EndY, StaZ, EndZ).
			Create_LineRbarDgn(1, dB1-dDc,     dDc, dHc-dDc, dHc-dDc, UmdRcRbar);
			Create_LineRbarDgn(2,     dDc,     dDc, dHc-dDc,     dDc, UmdRcRbar);
			Create_LineRbarDgn(1, dDc,     dB1-dDc,     dDc,     dDc, UmdRcRbar);
			Create_LineRbarDgn(2, dB1-dDc, dB1-dDc,     dDc, dHc-dDc, UmdRcRbar);
		}
		else
		{
			for(int i=0; i<CONST_UMDRC_iBAR_LAY; i++)
			{
				// Covering.
				double dDc = RbarData.dDc[i];
				if(dDc <= 0.0)	break;
				// Top, Left, Bot, Right (StaY, EndY, StaZ, EndZ).
				iGrup = 1;
				Create_LineRbarChk(iGrup, dB1-dDc,     dDc, dHc-dDc, dHc-dDc, TRUE, 
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);
				iGrup = 2;
				Create_LineRbarChk(iGrup,     dDc,     dDc, dHc-dDc,     dDc, FALSE,   
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 1;
				Create_LineRbarChk(iGrup, dDc,     dB1-dDc,     dDc,     dDc, TRUE,   
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);
				iGrup = 2;
				Create_LineRbarChk(iGrup, dB1-dDc, dB1-dDc,     dDc, dHc-dDc, FALSE,   
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);
			}
		}
		
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_P)	// P.
	{
		// Column, Brace.		
    double dHc=0., dB1=0.;      
		dHc = dSize[0];
		dB1 = dSize[1];
    
		if(bDesign)	// Design.
		{
			// Covering.
			double dDc = dUserDcDgn;
			// From Top (StaY, CenY, StaZ, CenZ).
			Create_CircRbarDgn(1, dHc/2.0, dHc/2.0, dHc-dDc, dHc/2.0, UmdRcRbar);
		}
		else	// Checking.
		{
			for(int i=0; i<CONST_UMDRC_iBAR_LAY; i++)
			{
				// Covering.
				double dDc = RbarData.dDc[i];
				if(dDc <= 0.0)	break;
				// From Top (StaY, CenY, StaZ, CenZ).					
				iGrup = 1;
				Create_CircRbarChk(iGrup, dHc/2.0, dHc/2.0, dHc-dDc, dHc/2.0, 
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);
			}
		}
		
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_SR)	// SR.
	{
		// Column, Brace.
		double dHc=0.;
		dHc = dSize[0];

		if(bDesign)	// Design.
		{
			// Covering.
			double dDc = dUserDcDgn;
			// From Top (StaY, CenY, StaZ, CenZ).
			Create_CircRbarDgn(1, dHc/2.0, dHc/2.0, dHc-dDc, dHc/2.0, UmdRcRbar);
		}
		else // Checking.
		{
			for(int i=0; i<CONST_UMDRC_iBAR_LAY; i++)
			{
				// Covering.
				double dDc = RbarData.dDc[i];
				if(dDc <= 0.0)	break;
				// From Top (StaY, CenY, StaZ, CenZ).
				iGrup = 1;
				Create_CircRbarChk(iGrup, dHc/2.0, dHc/2.0, dHc-dDc, dHc/2.0, 
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);
			}
		}		
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_SB)	// SB.
	{
		// Beam, Column, Brace.
		double dHc=0., dB1=0.;
    dHc = dSize[0];
		dB1 = dSize[1];
				
		if(bDesign)	// Design.
		{
			// Covering.
			double dDc = dUserDcDgn;
			// Top, Left, Bot, Right (StaY, EndY, StaZ, EndZ).
			Create_LineRbarDgn(1, dB1-dDc,     dDc, dHc-dDc, dHc-dDc, UmdRcRbar);
			Create_LineRbarDgn(2,     dDc,     dDc, dHc-dDc,     dDc, UmdRcRbar);
			Create_LineRbarDgn(1, dDc,     dB1-dDc,     dDc,     dDc, UmdRcRbar);
			Create_LineRbarDgn(2, dB1-dDc, dB1-dDc,     dDc, dHc-dDc, UmdRcRbar);
		}
		else	// Checking.
		{
			for(int i=0; i<CONST_UMDRC_iBAR_LAY; i++)
			{
				// Covering.
				double dDc = RbarData.dDc[i];
				if(dDc <= 0.0)	break;
				// Top, Left, Bot, Right (StaY, EndY, StaZ, EndZ).
				iGrup = 1;
				Create_LineRbarChk(iGrup, dB1-dDc,     dDc, dHc-dDc, dHc-dDc, TRUE,
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);				
				iGrup = 2;
				Create_LineRbarChk(iGrup,     dDc,     dDc, dHc-dDc,     dDc, FALSE,
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);		
				iGrup = 1;
				Create_LineRbarChk(iGrup, dDc,     dB1-dDc,     dDc,     dDc, TRUE,
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);		
				iGrup = 2;
				Create_LineRbarChk(iGrup, dB1-dDc, dB1-dDc,     dDc, dHc-dDc, FALSE,
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);		
			}
		}
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_OCT)	// OCT.
	{
		// Column, Brace.		
    double dH=0., dB=0., da=0., db=0., dt=0.;
    dH = dSize[0];
		dB = dSize[1];
		da = dSize[2];
		db = dSize[3];
    dt = dSize[4];
    
		// Calculate Inner Section Dimension.		  
		double dH1=0.0, dB1=0.0, da1=0.0, db1=0.0, dDist0=0.0, dTheta1=0.0;
		if(!CalcSectRegular_InOCT(dSize,dH1,dB1,da1,db1,dDist0,dTheta1))	return FALSE;

		double dCenY = dB/2.0;
		double dCenZ = dH/2.0;
		if(bDesign)	// Design.
		{
			// Covering.
			double dDc = dUserDcDgn;
			// Top, Diag, Right, Diag, Bot, Diag, Left, Diag (StaY, EndY, StaZ, EndZ).
			double dY1 = (dB/2.0-da) - ((dB/2.0-da)-(dB1/2.0-da1))*(dDc/dt);
			double dZ1 = dH/2.0 - dDc;
			double dY2 = dB/2.0 - dDc;
			double dZ2 = (dH/2.0-db) - ((dH/2.0-db)-(dH1/2.0-db1))*(dDc/dt);
			Create_LineRbarDgn(1,  dY1+dCenY, -dY1+dCenY,  dZ1+dCenZ,  dZ1+dCenZ, UmdRcRbar);
			Create_LineRbarDgn(2, -dY1+dCenY, -dY2+dCenY,  dZ1+dCenZ,  dZ2+dCenZ, UmdRcRbar);
			Create_LineRbarDgn(3, -dY2+dCenY, -dY2+dCenY,  dZ2+dCenZ, -dZ2+dCenZ, UmdRcRbar);
			Create_LineRbarDgn(2, -dY2+dCenY, -dY1+dCenY, -dZ2+dCenZ, -dZ1+dCenZ, UmdRcRbar);
			Create_LineRbarDgn(1, -dY1+dCenY,  dY1+dCenY, -dZ1+dCenZ, -dZ1+dCenZ, UmdRcRbar);
			Create_LineRbarDgn(2,  dY1+dCenY,  dY2+dCenY, -dZ1+dCenZ, -dZ2+dCenZ, UmdRcRbar);
			Create_LineRbarDgn(3,  dY2+dCenY,  dY2+dCenY, -dZ2+dCenZ,  dZ2+dCenZ, UmdRcRbar);
			Create_LineRbarDgn(2,  dY2+dCenY,  dY1+dCenY,  dZ2+dCenZ,  dZ1+dCenZ, UmdRcRbar);
		}
		else	// Checking.
		{
			for(int i=0; i<CONST_UMDRC_iBAR_LAY; i++)
			{
				// Covering.
				double dDc = RbarData.dDc[i];
				if(dDc <= 0.0)	break;
				// Top, Diag, Right, Diag, Bot, Diag, Left, Diag (StaY, EndY, StaZ, EndZ).
				double dY1, dZ1, dY2, dZ2;
				if(dt==0.0)
				{
					dY1 = (dB/2.0-da);
					dZ1 = dH/2.0 - dDc;
					dY2 = dB/2.0 - dDc;
					dZ2 = (dH/2.0-db);
				}
				else
				{
					dY1 = (dB/2.0-da) - ((dB/2.0-da)-(dB1/2.0-da1))*(dDc/dt);
					dZ1 = dH/2.0 - dDc;
					dY2 = dB/2.0 - dDc;
					dZ2 = (dH/2.0-db) - ((dH/2.0-db)-(dH1/2.0-db1))*(dDc/dt);
				}
				iGrup = 1;
				Create_LineRbarChk(iGrup,  dY1+dCenY, -dY1+dCenY,  dZ1+dCenZ,  dZ1+dCenZ, TRUE,
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 2;
				Create_LineRbarChk(iGrup, -dY1+dCenY, -dY2+dCenY,  dZ1+dCenZ,  dZ2+dCenZ, FALSE,
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 3;
				Create_LineRbarChk(iGrup, -dY2+dCenY, -dY2+dCenY,  dZ2+dCenZ, -dZ2+dCenZ, TRUE,
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 2;
				Create_LineRbarChk(iGrup, -dY2+dCenY, -dY1+dCenY, -dZ2+dCenZ, -dZ1+dCenZ, FALSE,
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 1;
				Create_LineRbarChk(iGrup, -dY1+dCenY,  dY1+dCenY, -dZ1+dCenZ, -dZ1+dCenZ, TRUE,
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 2;
				Create_LineRbarChk(iGrup,  dY1+dCenY,  dY2+dCenY, -dZ1+dCenZ, -dZ2+dCenZ, FALSE,
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 3;
				Create_LineRbarChk(iGrup,  dY2+dCenY,  dY2+dCenY, -dZ2+dCenZ,  dZ2+dCenZ, TRUE,
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 2;
				Create_LineRbarChk(iGrup,  dY2+dCenY,  dY1+dCenY,  dZ2+dCenZ,  dZ1+dCenZ, FALSE,
					                 RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
													 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
			}
		}		
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_SOCT)	// SOCT.
	{
		// Column, Brace.		
    double dH=0., dB=0., da=0., db=0.;
    dH = dSize[0];
		dB = dSize[1];
		da = dSize[2];
		db = dSize[3];
    
		double dCenY = dB/2.0;
		double dCenZ = dH/2.0;

		double dPi = 4.0*atan(1.0);
		double dTheta1 = (da==0.0) ? dPi/2.0 - cUMDRC_Zero : (dPi-atan(db/da))/2.0;	// Radian.
		double dTheta2 = (db==0.0) ? dPi/2.0 - cUMDRC_Zero : (dPi-atan(da/db))/2.0;	// Radian.
		if(bDesign)	// Design.
		{
			// Covering.
			double dDc = dUserDcDgn;
			// Top, Diag, Right, Diag, Bot, Diag, Left, Diag (StaY, EndY, StaZ, EndZ).
			double dY1 = (dB/2.0-da) - dDc/tan(dTheta1);
			double dZ1 = dH/2.0 - dDc;
			double dY2 = dB/2.0 - dDc;
			double dZ2 = (dH/2.0-db) - dDc/tan(dTheta2);
			Create_LineRbarDgn(1,  dY1+dCenY, -dY1+dCenY,  dZ1+dCenZ,  dZ1+dCenZ, UmdRcRbar);
			Create_LineRbarDgn(2, -dY1+dCenY, -dY2+dCenY,  dZ1+dCenZ,  dZ2+dCenZ, UmdRcRbar);
			Create_LineRbarDgn(3, -dY2+dCenY, -dY2+dCenY,  dZ2+dCenZ, -dZ2+dCenZ, UmdRcRbar);
			Create_LineRbarDgn(2, -dY2+dCenY, -dY1+dCenY, -dZ2+dCenZ, -dZ1+dCenZ, UmdRcRbar);
			Create_LineRbarDgn(1, -dY1+dCenY,  dY1+dCenY, -dZ1+dCenZ, -dZ1+dCenZ, UmdRcRbar);
			Create_LineRbarDgn(2,  dY1+dCenY,  dY2+dCenY, -dZ1+dCenZ, -dZ2+dCenZ, UmdRcRbar);
			Create_LineRbarDgn(3,  dY2+dCenY,  dY2+dCenY, -dZ2+dCenZ,  dZ2+dCenZ, UmdRcRbar);
			Create_LineRbarDgn(2,  dY2+dCenY,  dY1+dCenY,  dZ2+dCenZ,  dZ1+dCenZ, UmdRcRbar);
		}
		else	// Checking.
		{
			for(int i=0; i<CONST_UMDRC_iBAR_LAY; i++)
			{
				// Covering.
				double dDc = RbarData.dDc[i];
				if(dDc <= 0.0)	break;
				// Top, Diag, Right, Diag, Bot, Diag, Left, Diag (StaY, EndY, StaZ, EndZ).
				double dY1 = (dB/2.0-da) - dDc/tan(dTheta1);
				double dZ1 = dH/2.0 - dDc;
				double dY2 = dB/2.0 - dDc;
				double dZ2 = (dH/2.0-db) - dDc/tan(dTheta2);
				iGrup = 1;
				Create_LineRbarChk(iGrup,  dY1+dCenY, -dY1+dCenY,  dZ1+dCenZ,  dZ1+dCenZ, TRUE,
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 2;
				Create_LineRbarChk(iGrup, -dY1+dCenY, -dY2+dCenY,  dZ1+dCenZ,  dZ2+dCenZ, FALSE, 
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 3;
				Create_LineRbarChk(iGrup, -dY2+dCenY, -dY2+dCenY,  dZ2+dCenZ, -dZ2+dCenZ, TRUE,  
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 2;
				Create_LineRbarChk(iGrup, -dY2+dCenY, -dY1+dCenY, -dZ2+dCenZ, -dZ1+dCenZ, FALSE, 
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 1;
				Create_LineRbarChk(iGrup, -dY1+dCenY,  dY1+dCenY, -dZ1+dCenZ, -dZ1+dCenZ, TRUE,  
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 2;
				Create_LineRbarChk(iGrup,  dY1+dCenY,  dY2+dCenY, -dZ1+dCenZ, -dZ2+dCenZ, FALSE, 
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 3;
				Create_LineRbarChk(iGrup,  dY2+dCenY,  dY2+dCenY, -dZ2+dCenZ,  dZ2+dCenZ, TRUE,  
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 2;
				Create_LineRbarChk(iGrup,  dY2+dCenY,  dY1+dCenY,  dZ2+dCenZ,  dZ1+dCenZ, FALSE, 
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);						
			}
		}		
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_TRK)	// TRK
	{
		// Column, Brace.		
    double dH=0., dB=0., dt=0.;
    dH = dSize[0];
		dB = dSize[1];
		dt = dSize[2];
    

		double dCenY = dB/2.0;
		double dCenZ = dH/2.0;
		if(bDesign)	// Design.
		{
			// Covering.
			double dDc = dUserDcDgn;
			// Top, Left, Bot, Right (StaY, EndY, StaZ, EndZ).
			Create_LineRbarDgn(1,  ((dB-dH)/2.0)+dCenY, -((dB-dH)/2.0)+dCenY, dH-dDc, dH-dDc,        UmdRcRbar);
			Create_CirLRbarDgn(2,                dCenY,                dCenY, dH-dDc,  dCenZ, dB-dH, UmdRcRbar);
			Create_LineRbarDgn(1, -((dB-dH)/2.0)+dCenY,  ((dB-dH)/2.0)+dCenY,    dDc,    dDc,        UmdRcRbar);
			Create_CirRRbarDgn(2,                dCenY,                dCenY,    dDc,  dCenZ, dB-dH, UmdRcRbar);
		}
		else	// Checking.
		{
			for(int i=0; i<CONST_UMDRC_iBAR_LAY; i++)
			{
				// Covering.
				double dDc = RbarData.dDc[i];
				if(dDc <= 0.0)	break;
				// Top, Left, Bot, Right (StaY, EndY, StaZ, EndZ).
				iGrup = 1;
				Create_LineRbarChk(iGrup,  ((dB-dH)/2.0)+dCenY, -((dB-dH)/2.0)+dCenY, dH-dDc, dH-dDc,        TRUE, 
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);				
				iGrup = 2;
				Create_CirLRbarChk(iGrup,                dCenY,                dCenY, dH-dDc,  dCenZ, dB-dH,       
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);				
				iGrup = 1;
				Create_LineRbarChk(iGrup, -((dB-dH)/2.0)+dCenY,  ((dB-dH)/2.0)+dCenY,    dDc,    dDc,        TRUE, 
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);				
				iGrup = 2;
				Create_CirRRbarChk(iGrup,                dCenY,                dCenY,    dDc,  dCenZ, dB-dH,       
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);				
			}
		}		
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_STRK)	// STRK
	{
		// Column, Brace.		
    double dH=0., dB=0.;
    dH = dSize[0];
		dB = dSize[1];
    
		double dCenY = dB/2.0;
		double dCenZ = dH/2.0;
		if(bDesign)	// Design.
		{
			// Covering.
			double dDc = dUserDcDgn;
			// Top, Left, Bot, Right (StaY, EndY, StaZ, EndZ).
			Create_LineRbarDgn(1,  ((dB-dH)/2.0)+dCenY, -((dB-dH)/2.0)+dCenY, dH-dDc, dH-dDc,        UmdRcRbar);
			Create_CirLRbarDgn(2,                dCenY,                dCenY, dH-dDc,  dCenZ, dB-dH, UmdRcRbar);
			Create_LineRbarDgn(1, -((dB-dH)/2.0)+dCenY,  ((dB-dH)/2.0)+dCenY,    dDc,    dDc,        UmdRcRbar);
			Create_CirRRbarDgn(2,                dCenY,                dCenY,    dDc,  dCenZ, dB-dH, UmdRcRbar);
		}
		else	// Checking.
		{
			for(int i=0; i<CONST_UMDRC_iBAR_LAY; i++)
			{
				// Covering.
				double dDc = RbarData.dDc[i];
				if(dDc <= 0.0)	break;
				// Top, Left, Bot, Right (StaY, EndY, StaZ, EndZ).
				iGrup = 1;
				Create_LineRbarChk(iGrup,  ((dB-dH)/2.0)+dCenY, -((dB-dH)/2.0)+dCenY, dH-dDc, dH-dDc,        TRUE,  
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 2;
				Create_CirLRbarChk(iGrup,                dCenY,                dCenY, dH-dDc,  dCenZ, dB-dH,        
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 1;
				Create_LineRbarChk(iGrup, -((dB-dH)/2.0)+dCenY,  ((dB-dH)/2.0)+dCenY,    dDc,    dDc,        TRUE,  
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 2;
				Create_CirRRbarChk(iGrup,                dCenY,                dCenY,    dDc,  dCenZ, dB-dH,        
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
			}
		}		
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_HTRK)	// HTRK
	{
		// Column, Brace.
		double dH=0., dB=0.;
    dH = dSize[0];
		dB = dSize[1];
    
		double dCenY = (dB-dH/2.0)/2.0;
		double dCenZ = dH/2.0;
		if(bDesign)	// Design.
		{
			// Covering.
			double dDc = dUserDcDgn;
			// Top, Left, Bot, Right (StaY, EndY, StaZ, EndZ).
			Create_LineRbarDgn(1, (dB-dH/2.0),         dDc, dH-dDc, dH-dDc,            UmdRcRbar);
			Create_LineRbarDgn(3,         dDc,         dDc, dH-dDc,    dDc,            UmdRcRbar);
			Create_LineRbarDgn(1,         dDc, (dB-dH/2.0),    dDc,    dDc,            UmdRcRbar);
			Create_CirRRbarDgn(2,       dCenY,       dCenY,    dDc,  dCenZ, dB-dH/2.0, UmdRcRbar);
		}
		else // Checking.
		{
			for(int i=0; i<CONST_UMDRC_iBAR_LAY; i++)
			{
				// Covering.
				double dDc = RbarData.dDc[i];
				if(dDc <= 0.0)	break;
				// Top, Left, Bot, Right (StaY, EndY, StaZ, EndZ).
				iGrup = 1;
				Create_LineRbarChk(iGrup, (dB-dH/2.0),         dDc, dH-dDc, dH-dDc,            TRUE,  
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 3;
				Create_LineRbarChk(iGrup,         dDc,         dDc, dH-dDc,    dDc,            FALSE, 
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 1;
				Create_LineRbarChk(iGrup,         dDc, (dB-dH/2.0),    dDc,    dDc,            TRUE,  
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
				iGrup = 2;
				Create_CirRRbarChk(iGrup,       dCenY,       dCenY,    dDc,  dCenZ, dB-dH/2.0,        
					               RbarData.iBarNum[iGrup-1][i], RbarData.dBarDiaNa1[iGrup-1][i], RbarData.dBarDiaNa2[iGrup-1][i], 
												 RbarData.dBarAreaNa1[iGrup-1][i], RbarData.dBarAreaNa2[iGrup-1][i], UmdRcRbar);	
			}
		}	
	}
	else	ASSERT(0);

	int iRbarNum=UmdRcRbar.arRbarUnit.GetCount();
	
	if(iRbarNum > 0)	return TRUE;
	else							return FALSE;

	return FALSE;
}

BOOL CMdgnSectTool::MakeRbarDataGen(int iShape, const _UMD_RC_COL_MAINRBAR &RbarData, _UMD_RC_RBAR &UmdRcRbar)
{
	if(iShape!=DGN_SECT_SHAPE_INDEX_REG_GEN) { ASSERT(0); return FALSE; }

	UmdRcRbar.Init();

	int nSizeRbar = RbarData.arGenRbar.GetSize();
	for(int i=0; i<nSizeRbar; ++i)
	{
		_UMD_RC_GEN_RBAR_UNIT GenRbar = RbarData.arGenRbar[i];
		
		_UMD_RC_RBAR_UNIT RbarUnit;
		RbarUnit.dArea = GenRbar.dArea;
		RbarUnit.dDia  = GenRbar.dDia;
		RbarUnit.dyz[0]= GenRbar.dyz[0];
		RbarUnit.dyz[1]= GenRbar.dyz[1];		
		UmdRcRbar.arRbarUnit.SetAt(i, RbarUnit);
	}
	
	return TRUE;
}


double CMdgnSectTool::Get_CvrThkOnColmDgn(int iShape, double dSize[8], double dDcDgn)
{
	double dDc = dDcDgn;
	
  if(iShape==DGN_SECT_SHAPE_INDEX_REG_B)	// B.
	{
    double dTw, dTf1, dTf2;
		dTw = dSize[2];
		dTf1= dSize[3];
		dTf2= dSize[5];  
		dDc = (dDc==0.0 ? min(dTw,min(dTf1,dTf2))/2.0 : dDc);
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_P)	// P.
	{
    double dHc, dB1;
    dHc = dSize[0];
		dB1 = dSize[1];          
		dDc = (dDc==0.0 ? dB1/2.0 : dDc);
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_SR)	// SR.
	{		
    double dHc;
    dHc = dSize[0];    
		dDc = (dDc==0.0 ? dHc/8.0 : dDc);
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_SB)	// SB.
	{
		// Beam, Column, Brace.
    double dHc, dB1;
    dHc = dSize[0];
		dB1 = dSize[1];      
    dDc = (dDc==0.0 ? min(dHc,dB1)/8.0 : dDc);
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_OCT)	// OCT.
	{
    double dt;
    dt = dSize[4];
    dDc = (dDc==0.0 ? dt/2.0 : dDc);
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_SOCT)	// SOCT.
	{
    double dH, dB;
    dH = dSize[0];
		dB = dSize[1];      
    dDc = (dDc==0.0 ? min(dH,dB)/8.0 : dDc);
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_TRK)	// TRK
	{
    double dH, dB, dt;
    dH = dSize[0];
		dB = dSize[1];      
    dt = dSize[2];
    dDc = (dDc==0.0 ? dt/2.0 : dDc);
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_STRK)	// STRK
	{
    double dH, dB;
		dH = dSize[0];
		dB = dSize[1];      
    dDc = (dDc==0.0 ? min(dH,dB)/8.0 : dDc);
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_HTRK)	// HTRK
	{
    double dH, dB;
    dH = dSize[0];
		dB = dSize[1];      
    dDc = (dDc==0.0 ? min(dH,dB)/8.0 : dDc);
	}
	else if(iShape==DGN_SECT_SHAPE_INDEX_REG_GEN)
	{
		double dH, dB;
    dH = dSize[0];
		dB = dSize[1];      
    dDc = (dDc==0.0 ? min(dH,dB)/8.0 : dDc);
	}
  else ASSERT(0);
	
	return dDc;
}

BOOL CMdgnSectTool::Create_LineRbarDgn(int iGrup, double dStaY, double dEndY, double dStaZ, double dEndZ,
																										 _UMD_RC_RBAR& RbarData, int iBar_Div)
{
	if(iBar_Div <= 0) iBar_Div = CONST_UMDRC_iBAR_DIV;
	double dTotLen = sqrt(pow(dStaY-dEndY,2) + pow(dStaZ-dEndZ,2));
	if(fabs(dTotLen) < cUMDRC_Zero)	return FALSE;
	
	_UMD_RC_RBAR_UNIT RbarUnit;
	int iRbarCount = RbarData.arRbarUnit.GetCount();
	for(int i=0; i<iBar_Div; i++)
	{
		RbarUnit.Init();
		RbarUnit.iGrupNo = iGrup;
		RbarUnit.dLen		 = dTotLen / iBar_Div;
		RbarUnit.dyz[0]  = dStaY + (dEndY-dStaY)*((i+0.5)/iBar_Div);
		RbarUnit.dyz[1]  = dStaZ + (dEndZ-dStaZ)*((i+0.5)/iBar_Div);
		iRbarCount++;
		RbarData.arRbarUnit.SetAt(iRbarCount, RbarUnit);
	}
	return TRUE;
}

BOOL CMdgnSectTool::Create_LineRbarChk(int iGrup,	double dStaY, double dEndY, double dStaZ, double dEndZ, BOOL bIncEdge,
																										 int iDivNum, double dDiaNa1, double dDiaNa2, double dAreaNa1, double dAreaNa2, _UMD_RC_RBAR& RbarData)
{
	double dTotLen = sqrt(pow(dStaY-dEndY,2) + pow(dStaZ-dEndZ,2));
	if(fabs(dTotLen) < cUMDRC_Zero)	return FALSE;
	
	_UMD_RC_RBAR_UNIT RbarUnit;
	
	if(dDiaNa1==0.0 || dAreaNa1==0.0)	return FALSE;
	if(dDiaNa2==0.0)  dDiaNa2  = dDiaNa1;
	if(dAreaNa2==0.0)	dAreaNa2 = dAreaNa1;
	
	int iRbarCount = RbarData.arRbarUnit.GetCount();
	
	//int iTotDivNum = (IsBeamChk() ? iDivNum-1 : iDivNum);
	
  // changing the setting rebar form by SeungJun (`06.02.24)
  // if bIncEdge is TURE, with the rebar of the edge.
  // if bIncEdge is FALSE, without the rebar of the edge.
  int iTotDivNum = (bIncEdge ? iDivNum-1 : iDivNum+1);
  double RebarSeq = 0.0;
	
	if(iTotDivNum <= 0)	return FALSE;
	for(int i=0; i<iDivNum; i++)
	{
		RbarUnit.Init();
		RbarUnit.iGrupNo	= iGrup;
		RbarUnit.dDia	 = (i%2==0 ? dDiaNa1  : dDiaNa2);
		RbarUnit.dArea = (i%2==0 ? dAreaNa1 : dAreaNa2);
    RebarSeq = (bIncEdge ? (double)i : (double)i+1);
		RbarUnit.dyz[0]		= dStaY + (dEndY-dStaY)*(RebarSeq/iTotDivNum);
		RbarUnit.dyz[1]		= dStaZ + (dEndZ-dStaZ)*(RebarSeq/iTotDivNum);
		iRbarCount++;
		RbarData.arRbarUnit.SetAt(iRbarCount, RbarUnit);
	}
	return TRUE;
}

BOOL CMdgnSectTool::Create_CircRbarDgn(int iGrup, double dStaY, double dCenY, double dStaZ, double dCenZ,
																										 _UMD_RC_RBAR& RbarData)
{
	double dRadius = sqrt(pow(dStaY-dCenY,2) + pow(dStaZ-dCenZ,2));
	if(fabs(dRadius) < cUMDRC_Zero)	return FALSE;
	
	double dPi = 4.*atan(1.0);
	int iDivNum = 2*CONST_UMDRC_iBAR_DIV;
	
	_UMD_RC_RBAR_UNIT RbarUnit;
	int iRbarCount = RbarData.arRbarUnit.GetCount();
	for(int i=0; i<iDivNum; i++)
	{
		// Counterclockwise (from Bottom).
		RbarUnit.Init();
		RbarUnit.iGrupNo = iGrup;
		RbarUnit.dLen		 = 2*dPi*dRadius / iDivNum;
		RbarUnit.dyz[0]  = dCenY + (-1)*sin(2*dPi*((i+0.5)/iDivNum)) * dRadius;
		RbarUnit.dyz[1]  = dCenZ +      cos(2*dPi*((i+0.5)/iDivNum)) * dRadius;
		iRbarCount++;
		RbarData.arRbarUnit.SetAt(iRbarCount, RbarUnit);
	}
	return TRUE;
}

BOOL CMdgnSectTool::Create_CircRbarChk(int iGrup, double dStaY, double dCenY, double dStaZ, double dCenZ,
																										 int iDivNum, double dDiaNa1, double dDiaNa2, double dAreaNa1, double dAreaNa2, _UMD_RC_RBAR& RbarData)
{
	double dRadius = sqrt(pow(dStaY-dCenY,2) + pow(dStaZ-dCenZ,2));
	if(fabs(dRadius) < cUMDRC_Zero)	return FALSE;
	
	double dPi = 4.*atan(1.0);
	_UMD_RC_RBAR_UNIT RbarUnit;
	
	if(dDiaNa1==0.0 || dAreaNa1==0.0)	return FALSE;
	if(dDiaNa2==0.0)  dDiaNa2  = dDiaNa1;
	if(dAreaNa2==0.0)	dAreaNa2 = dAreaNa1;
	
	int iRbarCount = RbarData.arRbarUnit.GetCount();
	for(int i=0; i<iDivNum; i++)
	{
		// Counterclockwise (from Top).
		RbarUnit.Init();
		RbarUnit.iGrupNo	= iGrup;
		RbarUnit.dDia	 = (i%2==0 ? dDiaNa1  : dDiaNa2);
		RbarUnit.dArea = (i%2==0 ? dAreaNa1 : dAreaNa2);
		RbarUnit.dyz[0]		= dCenY + (-1)*sin(2*dPi*((double)i/iDivNum)) * dRadius;
		RbarUnit.dyz[1]		= dCenZ +      cos(2*dPi*((double)i/iDivNum)) * dRadius;
		iRbarCount++;
		RbarData.arRbarUnit.SetAt(iRbarCount, RbarUnit);
	}
	return TRUE;
}

BOOL CMdgnSectTool::Create_CirLRbarDgn(int iGrup,
																										 double dStaY, double dCenY, double dStaZ, double dCenZ, double dCenDist,
																										 _UMD_RC_RBAR& RbarData)
{
	double dRadius = sqrt(pow(dStaY-dCenY,2) + pow(dStaZ-dCenZ,2));
	if(fabs(dRadius) < cUMDRC_Zero)	return FALSE;
	
	double dPi = 4.*atan(1.0);
	int iDivNum = 1*CONST_UMDRC_iBAR_DIV;
	
	_UMD_RC_RBAR_UNIT RbarUnit;
	int iRbarCount = RbarData.arRbarUnit.GetCount();
	for(int i=0; i<iDivNum; i++)
	{
		// Counterclockwise (from Top).
		RbarUnit.Init();
		RbarUnit.iGrupNo = iGrup;
		RbarUnit.dLen		 = dPi*dRadius / iDivNum;
		RbarUnit.dyz[0]  = dCenY + (-1)*sin(dPi*((i+0.5)/iDivNum))*dRadius - (dCenDist/2.0);
		RbarUnit.dyz[1]  = dCenZ +      cos(dPi*((i+0.5)/iDivNum))*dRadius;
		iRbarCount++;
		RbarData.arRbarUnit.SetAt(iRbarCount, RbarUnit);
	}
	return TRUE;
}

BOOL CMdgnSectTool::Create_CirRRbarDgn(int iGrup,
																										 double dStaY, double dCenY, double dStaZ, double dCenZ, double dCenDist,
																										 _UMD_RC_RBAR& RbarData)
{
	double dRadius = sqrt(pow(dStaY-dCenY,2) + pow(dStaZ-dCenZ,2));
	if(fabs(dRadius) < cUMDRC_Zero)	return FALSE;
	
	double dPi = 4.*atan(1.0);
	int iDivNum = 1*CONST_UMDRC_iBAR_DIV;
	
	_UMD_RC_RBAR_UNIT RbarUnit;
	int iRbarCount = RbarData.arRbarUnit.GetCount();
	for(int i=0; i<iDivNum; i++)
	{
		// Counterclockwise (from Bot).
		RbarUnit.Init();
		RbarUnit.iGrupNo = iGrup;
		RbarUnit.dLen		 = dPi*dRadius / iDivNum;
		RbarUnit.dyz[0]  = dCenY +      sin(dPi*((i+0.5)/iDivNum))*dRadius + (dCenDist/2.0);
		RbarUnit.dyz[1]  = dCenZ + (-1)*cos(dPi*((i+0.5)/iDivNum))*dRadius;
		iRbarCount++;
		RbarData.arRbarUnit.SetAt(iRbarCount, RbarUnit);
	}
	return TRUE;
}

BOOL CMdgnSectTool::Create_CirLRbarChk(int iGrup,
																										 double dStaY, double dCenY, double dStaZ, double dCenZ, double dCenDist,
																										 int iDivNum, double dDiaNa1, double dDiaNa2, double dAreaNa1, double dAreaNa2, _UMD_RC_RBAR& RbarData)
{
	double dRadius = sqrt(pow(dStaY-dCenY,2) + pow(dStaZ-dCenZ,2));
	if(fabs(dRadius) < cUMDRC_Zero)	return FALSE;
	
	double dPi = 4.*atan(1.0);
	_UMD_RC_RBAR_UNIT RbarUnit;
	
	if(dDiaNa1==0.0 || dAreaNa1==0.0)	return FALSE;
	if(dDiaNa2==0.0)  dDiaNa2  = dDiaNa1;
	if(dAreaNa2==0.0)	dAreaNa2 = dAreaNa1;
	
	// changing the setting rebar form by SeungJun (`06.02.24)
  // left circle, always withour the rebar of the edge
  // so iDivNum+1
  int iRbarCount = RbarData.arRbarUnit.GetCount();
	for(int i=0; i<iDivNum; i++)
	{
		// Counterclockwise (from Top).
		RbarUnit.Init();
		RbarUnit.iGrupNo	= iGrup;
		RbarUnit.dDia	 = (i%2==0 ? dDiaNa1  : dDiaNa2);
		RbarUnit.dArea = (i%2==0 ? dAreaNa1 : dAreaNa2);
		RbarUnit.dyz[0]		= dCenY + (-1)*sin(dPi*(((double)i+1.)/(iDivNum+1))) * dRadius - (dCenDist/2.0);
		RbarUnit.dyz[1]		= dCenZ +      cos(dPi*(((double)i+1.)/(iDivNum+1))) * dRadius;
		iRbarCount++;
		RbarData.arRbarUnit.SetAt(iRbarCount, RbarUnit);
	}
	return TRUE;
}

BOOL CMdgnSectTool::Create_CirRRbarChk(int iGrup,
																										 double dStaY, double dCenY, double dStaZ, double dCenZ, double dCenDist,
																										 int iDivNum, double dDiaNa1, double dDiaNa2, double dAreaNa1, double dAreaNa2, _UMD_RC_RBAR& RbarData)
{
	double dRadius = sqrt(pow(dStaY-dCenY,2) + pow(dStaZ-dCenZ,2));
	if(fabs(dRadius) < cUMDRC_Zero)	return FALSE;
	
	double dPi = 4.*atan(1.0);
	_UMD_RC_RBAR_UNIT RbarUnit;
	
	if(dDiaNa1==0.0 || dAreaNa1==0.0)	return FALSE;
	if(dDiaNa2==0.0)  dDiaNa2  = dDiaNa1;
	if(dAreaNa2==0.0)	dAreaNa2 = dAreaNa1;
	
  // changing the setting rebar form by SeungJun (`06.02.24)
  // left circle, always withour the rebar of the edge
  // so iDivNum+1
	int iRbarCount = RbarData.arRbarUnit.GetCount();
	for(int i=0; i<iDivNum; i++)
	{
		// Counterclockwise (from Bot).
		RbarUnit.Init();
		RbarUnit.iGrupNo	= iGrup;
		RbarUnit.dDia	 = (i%2==0 ? dDiaNa1  : dDiaNa2);
		RbarUnit.dArea = (i%2==0 ? dAreaNa1 : dAreaNa2);
		RbarUnit.dyz[0]		= dCenY +      sin(dPi*(((double)i+1.)/(iDivNum+1))) * dRadius + (dCenDist/2.0);
		RbarUnit.dyz[1]		= dCenZ + (-1)*cos(dPi*(((double)i+1.)/(iDivNum+1))) * dRadius;
		iRbarCount++;
		RbarData.arRbarUnit.SetAt(iRbarCount, RbarUnit);
	}
	return TRUE;
}

BOOL CMdgnSectTool::Calc_SectGeneral(_UMD_RC_GSEC_POLYGON_LIST &aOutPolyData,  _UMD_RC_GSEC_POLYGON_LIST &aInPolyData, double &rdArea, double &rdYbar, double &rdZbar)
{
	rdArea = 0.0;
	rdYbar = 0.0;
	rdZbar = 0.0;
	
	// 도심점 및 경계값 산출
	double dXmax, dXmin, dYmax, dYmin;
	dXmax = dXmin = dYmax = dYmin = 0.0;
	
	BOOL bFirst = TRUE;
	int nOutSize = aOutPolyData.GetSize();
	int nInSize  = aInPolyData.GetSize();
	
	CArray<double, double>                    arOutArea;
	CArray<double, double>                    arInArea;
	CArray<_UMD_RC_GSEC_VERTEX, _UMD_RC_GSEC_VERTEX&> arOutCenter;
	CArray<_UMD_RC_GSEC_VERTEX, _UMD_RC_GSEC_VERTEX&> arInCenter;
	arOutArea.SetSize(nOutSize);  arOutCenter.SetSize(nOutSize);
	arInArea.SetSize(nOutSize);   arInCenter.SetSize(nOutSize);
	
	double dTotalArea = 0.0;
	double dCenterX  = 0.0; 
	double dCenterY  = 0.0;
	
double dPolyArea, dPolyXmax, dPolyXmin, dPolyYmax, dPolyYmin;
	_UMD_RC_GSEC_VERTEX vertexD;
	for(int i=0 ; i<nOutSize ; i++)
	{
		dPolyArea = vertexD.dX = vertexD.dY = 0.0;
		
		Get_AreaCenter(vertexD, dPolyArea, aOutPolyData[i]);		
		aOutPolyData[i].GetBoundary(dPolyXmax, dPolyXmin, dPolyYmax, dPolyYmin);
		
		if(bFirst) {dXmax = dPolyXmax;             dXmin = dPolyXmin;             dYmax = dPolyYmax;             dYmin = dPolyYmin;          bFirst = FALSE;}
		else       {dXmax = max(dXmax,dPolyXmax);  dXmin = min(dXmin,dPolyXmin);  dYmax = max(dYmax,dPolyYmax);  dYmin = min(dYmin,dPolyYmin);}
		
		dTotalArea += dPolyArea;
		dCenterX  += dPolyArea*vertexD.dX;		
		dCenterY  += dPolyArea*vertexD.dY;
	}
	for(int i=0 ; i<nInSize ; i++)
	{
		dPolyArea = vertexD.dX = vertexD.dY = 0.0;
		
		Get_AreaCenter(vertexD, dPolyArea, aInPolyData[i]);	
		aInPolyData[i].GetBoundary(dPolyXmax, dPolyXmin, dPolyYmax, dPolyYmin);
		
		if(bFirst) {dXmax = dPolyXmax;             dXmin = dPolyXmin;             dYmax = dPolyYmax;             dYmin = dPolyYmin;          bFirst = FALSE;}
		else       {dXmax = max(dXmax,dPolyXmax);  dXmin = min(dXmin,dPolyXmin);  dYmax = max(dYmax,dPolyYmax);  dYmin = min(dYmin,dPolyYmin);}
		
		dTotalArea -= dPolyArea;
		dCenterX  -= dPolyArea*vertexD.dX;		
		dCenterY  -= dPolyArea*vertexD.dY;
	}
	
	if(dTotalArea > 0.0) { dCenterX /= dTotalArea;  dCenterY /= dTotalArea; }
	else                 { dCenterX = dCenterY = 0.0; }
	
	rdYbar = dCenterX-dXmin;//좌측으로부터의 거리
	rdZbar = dCenterY-dYmin;//하단으로부터의 거리
	rdArea = dTotalArea;	
		
	return TRUE;
}

BOOL CMdgnSectTool::Get_AreaCenter(_UMD_RC_GSEC_VERTEX& PointCenV, double& dArea, _UMD_RC_GSEC_POLYGON& polygonD)
{
	dArea = 0.0;
	PointCenV.dX = PointCenV.dY = 0.0;
	
	int nPointSize = polygonD.aVertex.GetSize();
	
	if(nPointSize >= 3)
	{
		double  *arX;  arX = new double[nPointSize];
		double  *arY;  arY = new double[nPointSize];
		
		for(int i=0 ; i<nPointSize ; i++)
		{
			arX[i] = polygonD.aVertex[i].dX;  arY[i] = polygonD.aVertex[i].dY;
		}
		CMathFunc::mathPolyCentroid(nPointSize, arX, arY, PointCenV.dX, PointCenV.dY, dArea);
		delete arX;
		delete arY;
		return TRUE;
	}		
	return FALSE;
}