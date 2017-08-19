// MdgnSectTool.h: interface for the CMdgnSectTool class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MDGNSECTTOOL_H__89B02AA1_2428_41FC_A5EA_BC5F21003957__INCLUDED_)
#define AFX_MDGNSECTTOOL_H__89B02AA1_2428_41FC_A5EA_BC5F21003957__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StructDef.h"
#include "PolyMaker.h"

class CMdgnSectTool  
{
public:
	CMdgnSectTool();
	virtual ~CMdgnSectTool();

public:
	BOOL Get_RotateSectData(int nType, double center_y, double center_z, int nSliceNum, int nDivisionNumber_90deg, 
		_UMD_RC_CON_POLY_DATA& PolyData, _UMD_RC_RBAR& RbarData, _UMD_RC_COL_ROTATE_SECT& RotSectData);
	
	int Get_CheckedSymmetryType(_UMD_RC_CON_POLY_DATA& ConPolyData, _UMD_RC_RBAR& RbarUnitData, double dCenter_y, double dCenter_z);

	BOOL Get_RotateRebarData(double center_y, double center_z, double rotation, _UMD_RC_RBAR& RbarData, _UMD_RC_RBAR_LIST& RotateRbarData);

	BOOL Get_RotateConPolyData(double center_y, double center_z, double rotation, _UMD_RC_CON_POLY_DATA& PolyData, double& max_y, double& min_y, double& max_z, double& min_z, CPolyMaker* pPolyMaker);

	BOOL Get_RotateConUnitPartList(int nType, double center_y, double center_z, double rotation, int nSliceNum, _UMD_RC_CON_POLY_DATA& PolyData, _UMD_RC_CON_PART_LIST& arConUnit);

	// 콘크리트 및 철근의 UnitList를 보유할 갯수를 넘겨줌 //참고문서 : 신봉진/RC_Column_MPhiCurve자료.ppt (2Page) 
  int Get_HoldNumberToUnitList(int iSymmetryType, int nDivisionNumber_90deg);
	
	// iRot번째에 해당하는 중립축의 회전각도를 넘겨줌
  //   iRot : 단면의 몇번째 회전시켰는지 횟수(0 ~ 2*분할수-1)  
	double Get_RotationAngle(int iRot, int iSymmetryType, int nDivisionNumber_90deg);

	BOOL Get_PolyDataCentroid(double dZ1, double dZ2, double dYmax, double dYmin, CPolyMaker* pPolyMaker, double& center_y, double& center_z, double& dArea_Inner);

	//------------------------------------------------------------------------------------------------------------
	// Make Section
	BOOL MakeConDataReg(int iShape, double dSize[8], _UMD_RC_CONC& UmdRcConc, int iDivNum);

	BOOL MakeRbarDataReg(int iShape, double dSize[8], const _UMD_RC_COL_MAINRBAR &RbarData, _UMD_RC_RBAR& UmdRcRbar);
	BOOL MakeRbarDataGen(int iShape, const _UMD_RC_COL_MAINRBAR &RbarData, _UMD_RC_RBAR &UmdRcRbar);
	

	double Get_CvrThkOnColmDgn(int iShape, double dSize[8], double dDcDgn);

	BOOL CalcSectRegular_InOCT(double dSize[8], double& dH1, double& dB1, double& da1, double& db1, double& dDist0, double& dTheta1);

	BOOL Create_LineRbarDgn(int iGrup, double dStaY, double dEndY, double dStaZ, double dEndZ,                  _UMD_RC_RBAR& RbarData, int iBar_Div = -1);

	BOOL Create_LineRbarChk(int iGrup,	double dStaY, double dEndY, double dStaZ, double dEndZ, BOOL bIncEdge,
																										 int iDivNum, double dDiaNa1, double dDiaNa2, double dAreaNa1, double dAreaNa2, _UMD_RC_RBAR& RbarData);

	BOOL Create_CircRbarDgn(int iGrup, double dStaY, double dCenY, double dStaZ, double dCenZ, _UMD_RC_RBAR& RbarData);
	BOOL Create_CircRbarChk(int iGrup, double dStaY, double dCenY, double dStaZ, double dCenZ,
																										 int iDivNum, double dDiaNa1, double dDiaNa2, double dAreaNa1, double dAreaNa2, _UMD_RC_RBAR& RbarData);


	BOOL Create_CirLRbarDgn(int iGrup, double dStaY, double dCenY, double dStaZ, double dCenZ, double dCenDist, _UMD_RC_RBAR& RbarData);
	BOOL Create_CirRRbarDgn(int iGrup, double dStaY, double dCenY, double dStaZ, double dCenZ, double dCenDist, _UMD_RC_RBAR& RbarData);

	BOOL Create_CirLRbarChk(int iGrup,   double dStaY, double dCenY, double dStaZ, double dCenZ, double dCenDist, 
		int iDivNum, double dDiaNa1, double dDiaNa2, double dAreaNa1, double dAreaNa2, _UMD_RC_RBAR& RbarData);
	BOOL Create_CirRRbarChk(int iGrup,   double dStaY, double dCenY, double dStaZ, double dCenZ, double dCenDist, 
		int iDivNum, double dDiaNa1, double dDiaNa2, double dAreaNa1, double dAreaNa2, _UMD_RC_RBAR& RbarData);

	// temp.
	BOOL Calc_SectGeneral(_UMD_RC_GSEC_POLYGON_LIST &aOutPolyData, _UMD_RC_GSEC_POLYGON_LIST &aInPolyData, double &rdArea, double &rdYbar, double &rdZbar);
	BOOL Get_AreaCenter(_UMD_RC_GSEC_VERTEX& PointCenV, double& dArea, _UMD_RC_GSEC_POLYGON& polygonD);

};

#endif // !defined(AFX_MDGNSECTTOOL_H__89B02AA1_2428_41FC_A5EA_BC5F21003957__INCLUDED_)
