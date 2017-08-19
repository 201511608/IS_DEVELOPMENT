// MdgnPmTool.h: interface for the CMdgnPmTool class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MDGNPMTOOL_H__A8A50FB1_A2F6_46A6_B91E_B24EBF392B8D__INCLUDED_)
#define AFX_MDGNPMTOOL_H__A8A50FB1_A2F6_46A6_B91E_B24EBF392B8D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StructDef.h"

class CMdgnPmTool  
{
public:
	CMdgnPmTool();
	virtual ~CMdgnPmTool();

	
public:
	double m_dTolAngle;// Tolerance for Angle (정밀 PM-Curve 산출시 수렴오차 default=0.1˚)
	double m_dTolM0;   // 최대 모멘트 산출시의 인접 중립축위치 산출및 모멘트 오차비(모멘트 오차 = m_dTolM0*최대모멘트). (default=0.001)  
	double m_dTolP0;   // 축하중 오차비(축하중 오차 = m_dTolP0*최대축하중). (default=0.001)  

	void Set_Tolerance(double dTolAngle, double dTolM0, double dTolP0);
	
	// 설계결과 정보(T_PMCV_3D)를 CDgnCalcBase_3DPM_Tool용 정보(_UMD_3DPM_DATA)로 변환하여줌
	//   nType : Phi고려여부 (1:Phi미고려, 2:Phi고려)
	BOOL Get_ConvertToPMTool(int nType, T_PMCV_3D& InData, _UMD_3DPM_DATA& OutData);

	// 3D-PM전체 결과를 이용하여 외력과 동일방향(Ps,dMs_y,dMs_z 기준)에서의 PM-Curve상의 점을 찾아서 넘겨줌
	//   Umd3DPmData     : 3DPM 상관도 정보
	//   dPs, dMs_y, dMs_z : 외력방향의 시작점의 축하중, y축에 대한 모멘트, z축에 대한 모멘트(Umd3DPmData로 구성된 망의 내부에 있어야함)
  //   dPd, dMd_y, dMd_z : 외력의 축하중, y축에 대한 모멘트, z축에 대한 모멘트
  // < dPn, dMn_y, dMn_z, desi, dXn : 외력과 동일한 방향의 PM-Curve상의 축하중, y축에 대한 모멘트, z축에 대한 모멘트를 넘겨 받을 변수
  // < bCheck            : 외력에 견딜수 있는지 여부를 넘겨 받을 변수
	// < iLeftFitDirID, dModul_DirID : 교차점의 위치를 알려주는 Umd3DPmData의 ID와 영향률을 넘겨 받을 변수
	//   bChkInStartPMM : 시작점(Ps,dMs_y,dMs_z)이 내부에 있는지를 체크할지 여부	
	// < bResInStartPMM : 시작점(Ps,dMs_y,dMs_z)이 내부에 있는지를 체크결과	(bChkInStartPMM=TRUE일때만 검토함)
	// ※참고문서 : 신봉진/RC_Column_MPhiCurve자료.ppt (3Page) 
	// ※주의사항 : 꼭 확인할 것!!!
	//              -원점, 시작점, 작용외력이 MyMz평면의 동일 각도상에 있지 않을 경우 꼭 실제 3D-PM(Φ값이 적용된)를 사용하여야 한다.
	//              -시작점(Ps,dMs_y,dMs_z)은 3DPM 내부에 있어야만 바른 결과가 나옴, 
	//               시작점이 외부에 있을지도 모를 경우에는 반드시 bChkInStartPMM를 TRUE로 두어 Check를 수행하여야 함
  BOOL Get_PmcvUnitAt3DPmcvData(_UMD_3DPM_DATA& Umd3DPmData, double dPs, double dMs_y, double dMs_z, 
		double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, 
		BOOL& bCheck, int& iLeftFitDirID, double& dModul_DirID, BOOL bChkInStartPMM, BOOL& bResInStartPMM);
	// 이미 계산된 M-Φ전체 결과를 이용하여 외력과 동일방향(원점기준)에서의 PM-Curve상의 점을 찾아서 넘겨줌
	BOOL Get_PmcvUnitAt3DPmcvData(_UMD_3DPM_DATA& Umd3DPmData, double dPd, double dMd_y, double dMd_z, 
		                            double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, BOOL& bCheck, int& iLeftFitDirID, double& dModul_DirID);

	BOOL Get_PmcvDataToDirID(_UMD_3DPM_DATA& Umd3DPmData, int iAxisDir_Global, double dModul_DirID, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcv, _UMD_PMCV_UNIT_BALANCE& PmBalance);

	BOOL Get_PiercePointForContactPlane(_UMD_3DPM_DATA& Umd3DPmData, double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, int iFitDirID, int iFitUnitID, int& iLeftFitDirID, double& dModul_DirID);

	BOOL Get_FitPmcvUnitID(_UMD_3DPM_DATA& Umd3DPmData, double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, int& iFitDirID, int& iFitUnitID);

	BOOL Get_PiercePontForDirection(_UMD_3DPM_DATA& Umd3DPmData, double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, int iDirID_1, int iDirID_2, double& dModul_DirID);

	void Get_ResultPmcvDataID(_UMD_3DPM_DATA& Umd3DPmData, int iAxisDir_Global, int& iPmcvDataID, double& dModul_My, double& dModul_Mz);

	_UMD_PMCV_UNIT_POLY Get_CngSignPmcvUnit(_UMD_PMCV_UNIT_POLY& pmcvUnit_In, double dModul_My, double dModul_Mz);  

	double Get_RotationGlobalAngle(int iDivisionNumber_90deg, int iAxisDir_Global);
	double Get_RotationGlobalAngle(int iDivisionNumber_90deg, int iAxisDir_Global, double dModul_DirID);

	// 해당 외력방향에 대한 Line이 arPmcvUnit내를 꾀뚫는지 여부 넘겨줌
	//   dPs, dMs_y, dMs_z : 외력방향의 시작점의 축하중, y축에 대한 모멘트, z축에 대한 모멘트
  //   dPd, dMd_y, dMd_z : 외력의 축하중, y축에 대한 모멘트, z축에 대한 모멘트
  // < arPmcvUnit: 해당 외력방향에 있는지 확인할 3D좌표망을 구성하는 배열
  // < RETURN: 관툥여부
  BOOL Cal_PierceCheckToPolyLine(double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit);
	// 원점에서 부터 해당 외력방향에 대한 Line이 arPmcvUnit내를 꾀뚫는지 여부 넘겨줌
	BOOL Cal_PierceCheckToPolyLine(double dPd, double dMd_y, double dMd_z, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit);

	// 해당 외력방향에 대한 Line이 3각형의 arPmcvUnit이 나타내는 평면과의 교차점을 계산하여 넘겨줌(방향성을 가짐)
	//   dPs, dMs_y, dMs_z : 외력방향의 시작점의 축하중, y축에 대한 모멘트, z축에 대한 모멘트
  //   dPd, dMd_y, dMd_z : 외력의 축하중, y축에 대한 모멘트, z축에 대한 모멘트
  // < dPn, dMn_y, dMn_z : 외력과 동일한 방향의 평면상의 교차점에 대한 축하중, y축에 대한 모멘트, z축에 대한 모멘트를 넘겨 받을 변수
	// < desi, dXn         : Pn,Mn시에서 최대 인장철근의 변형율 (압축시 '-'), 중립축의 직각방향으로의 위치(하단에서부터의 길이) 
  // < arPmcvUnit: 해당 외력방향과 교차점을 찾을 3D좌표망을 구성하는 3점의 배열
  // < RETURN: 교차점을 찾을수 없거나 교차점이 외력방향에 없을 경우 FALSE를 넘김 
  BOOL Cal_PiercePointToPolyLine(double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit);
	// 원점에서 부터 해당 외력방향에 대한 Line이 3각형의 arPmcvUnit이 나타내는 평면과의 교차점을 계산하여 넘겨줌(방향성을 가짐)
	BOOL Cal_PiercePointToPolyLine(double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit);
	
	// 해당 외력방향에 대한 Line이 3각형의 arPmcvUnit이 나타내는 평면과의 교차점을 계산하여 넘겨줌(Line의 교차점임)
  //   dPs, dMs_y, dMs_z : 외력방향의 시작점의 축하중, y축에 대한 모멘트, z축에 대한 모멘트
	//   dPd, dMd_y, dMd_z : 외력의 축하중, y축에 대한 모멘트, z축에 대한 모멘트
  // < dPn, dMn_y, dMn_z : 외력과 동일한 방향의 평면상의 교차점에 대한 축하중, y축에 대한 모멘트, z축에 대한 모멘트를 넘겨 받을 변수
	// < desi, dXn         : Pn,Mn시에서 최대 인장철근의 변형율 (압축시 '-'), 중립축의 직각방향으로의 위치(하단에서부터의 길이) 
  // < arPmcvUnit: 해당 외력방향과 교차점을 찾을 3D좌표망을 구성하는 3점의 배열
	BOOL Cal_CrossPointToPolyLine(double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit);
	// 원점에서 해당 외력방향에 대한 Line이 3각형의 arPmcvUnit이 나타내는 평면과의 교차점을 계산하여 넘겨줌(Line의 교차점임)
  BOOL Cal_CrossPointToPolyLine(double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit);

	BOOL Get_UnitVectorToMzMyP(double dP, double dMy, double dMz, double vectorMzMyP[3], double vectorMP[2], double vectorMzMy[2]);

	void Get_AngleMinIDToPmCurve_NonScale(_UMD_3DPM_DATA& Umd3DPmData, double nMzMyP_Org[3], double nVectorMzMyP[3], int iAxisDir_Global, double dPmax, double dMmax, int& iUnitID, double& dAngleMin);

	BOOL Cal_MathPierceCheckToPolyLine(double p0[3], double p1[3], const int nData, double polyLine[][3]);

	BOOL Cal_MathPierceCheckToPolyLine(double p1[3], const int nData, double polyLine[][3]);


};

#endif // !defined(AFX_MDGNPMTOOL_H__A8A50FB1_A2F6_46A6_B91E_B24EBF392B8D__INCLUDED_)
