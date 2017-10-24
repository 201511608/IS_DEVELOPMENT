// CalcPmcv.h: interface for the CCalcPmcv class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CALCPMCV_H__B8D92DCA_7770_4AA9_B39F_989A6821BE48__INCLUDED_)
#define AFX_CALCPMCV_H__B8D92DCA_7770_4AA9_B39F_989A6821BE48__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StructDef.h"
#include "PolyMaker.h"

class CCalcPmcv  
{
public:
	CCalcPmcv();
	virtual ~CCalcPmcv();

public:
	BOOL ExecuteCheck(BOOL bAllAxis);

	//----------------------------------------------------------------------------
	// SET INPUT DATA.
	//----------------------------------------------------------------------------

	//----------------------------------------------------------------------------
	// GET OUTPUT DATA.
	//----------------------------------------------------------------------------
	// [1] 3D P-M Curve Data
	BOOL Get3DPMCurveData(T_PMCV_3D &rPmcv3D);	
	// [2] 2D P-M Curve Data for specific rotate angle.
	BOOL Get2DPMCurveData(double dAngle, T_PMCV_2D &rPmcv2D);	
	// [3] (Pn, Mn), (phiPn, phiMn) for specific force.
	BOOL GetNomStrength(double dXn, double dAngle, T_PM_SECT_RES &rRes);
	BOOL GetDgnStrengthByEccnRatio(double dPu, double dMuy, double dMuz, const T_PMCV_2D &pmcv2D, T_PM_CHK_RES &rRes);	
	BOOL GetDgnStrengthByPu(double dPu, double dMuy, double dMuz, const T_PMCV_2D &pmcv2D, T_PM_CHK_RES &rRes);	

protected:
	void SetDesignFlag(BOOL bDesign) { m_bDesign = bDesign; }
	BOOL IsGeneralType();
	double GetRebarArea();
	BOOL Get3DPMCurveDataReg(double dAsThk, T_PMCV_3D& PmData3D);
	BOOL Get3DPMCurveDataGen(double dAsThk, T_PMCV_3D& PmData3D);

	BOOL GetRealPMCurveData(double dAsThk, T_PMCV_2D& PmcvData, T_PM_SECT_RES& CbPmData, double& dPhiPnmax, int nCheckType=0);
	// 중립축의 회전각에 대한 PM상관도를 계산하여줌
  BOOL Get_PMCurve(double dInRota, T_PMCV_2D& PmcvData, T_PM_SECT_RES& CbPmData, double& dPhiPnmax, double dAsThk);
	// General단면의 중립축의 회전각에 대한 PM상관도를 계산하여줌
	// ※주의 : 1회 계산시에만 사용 (반복계산시에는 속도가 저하됨)
	BOOL Get_PMCurve_Gen(double dInRota, T_PMCV_2D& PmcvData, T_PM_SECT_RES& CbPmData, double& dPhiPnmax, double dAsThk, BOOL bRbarDetail=FALSE);
	// 중립축의 회전각에 대한 PM상관도를 계산하여줌
  BOOL Get_PMCurve_Gen(int iAxisDir, T_PMCV_2D& PmcvData, T_PM_SECT_RES& CbPmData, T_PM_SECT_RES& P0PmData, double& dPhiPnmax, double dAsThk, BOOL bRbarDetail=FALSE);
	BOOL Check_PnMn(double dAsCur, T_PMCV_2D& PmcvData, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, double& dRatio, T_PM_SECT_RES* pResData=NULL);
	BOOL Check_PnMn_Gen(_UMD_3DPM_DATA& PmData3D, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, double& dRatio, T_PMCV_2D& PmcvData, T_PM_SECT_RES& CbPmData, T_PM_SECT_RES* pResData=NULL);

	// DB단면용 ETC 함수 ============================================================================================
	double Calc_Rota(double dInRota=-1.0);// dInRota(Radian)
	double Calc_Dmax(double dInRota=-1.0);// dInRota(Radian)

	// PM상관도로부터 작용력의 안전율 산출함
	double Calc_Ratio(double dPu, double dMuy, double dMuz, double dPhiPnmax, T_PMCV_2D& PmData, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, BOOL bPureMonemt=FALSE);
	double Calc_Ratio(double dPu, double dMuy, double dMuz, double dPhiPnmax, T_PMCV_2D& PmData, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, T_PM_SECT_RES& ResData, BOOL bPureMonemt=FALSE);
	double Calc_Ratio_ToPu(double dPu, double dMuy, double dMuz, double dPhiPnmax, T_PMCV_2D& PmData, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, T_PM_SECT_RES& ResData);

	// 해당 회전각에 대해 회전된 단면좌표의 유효깊이를 넘겨줌(단면좌표는 회전중심점을 중심으로 중립축 회전각의 역방향으로 회전한 값을 가짐) 
  //   dDmax     : 단면의 중립축 직각방향으로의 길이 
  //   dZmim     : 해당 회전각에 대해 회전된 단면좌표의 최소 Z값(단면좌표는 회전중심점을 중심으로 중립축 회전각의 역방향으로 회전한 값을 가짐) 
  //   arRbarUnit: 철근의 위치정보(중복을 피하기위해 원본을 넘겨줌:변동없음)
	//   bFixDgnRebarOnly : Dgn일 경우 고정단면적 철근만을 고려한 거리를 넘겨줄지 여부(Dgn이고 가상철근두께 t가 0.0일때에만 적용함. 고정단면적철근 적용시에만 사용함)
	double Calc_Deff_Gen(double dDmax, double dZmin, _UMD_RC_RBAR_LIST& arRbarUnit, BOOL bFixDgnRebarOnly = FALSE);
  // 유효깊이 계산해서 넘겨줌. 4개( 0=Y-Dir-Positive, 1=Y-Dir-Negative, 2=Z-Dir-Positive, 3=Z-Dir-Negative)
  double Calc_Deff_Gen(int iDir, BOOL bFixDgnRebarOnly = FALSE);    

	// General단면의 중립축의 위치,각도가 정해졌을때의 단면의 강도 계산(전체 계산시 : 일정 각도간격에 대한 반복계산시)
	// (비선형 응력분포를 중립축에 평행하게 Slice친 단면 정보를 이용하여 계산)
	//   dXn, dInRota : 중립축의 위치비(0.0~1.0, 0.0순수압축), 중립축의 회전각(계산된 강도를 실제축의 강도로 변환하기 위해서만 사용함)
	//   dAsThk       : Design상태(필요철근량사정)일때 철근띠의 두께
	// < ResData      : 단면의 강도정보
	//   arConcUnit   : 중립축의 회전각의 반대방향으로 회전된 Slice친 단면정보
	//   arRbarUnit   : 중립축의 회전각의 반대방향으로 회전된 철근정보
	// ※ 단면정보를 중립축의 회전각의 반대방향으로 회전시킴으로서 중립축은 수평으로 만들어 부재강도를 계산하며, 
	//    계산된 강도는 회전각(dInRota)에 대해 다시 변환하여 실제축의 강도값으로 변환하여 돌려주는 방식임
	BOOL   Calc_PnMn_Gen_SliceMethod(double dXn, double dInRota, double dAsThk, T_PM_SECT_RES& ResData, 
		_UMD_RC_CON_PART_LIST& arConcUnit, _UMD_RC_RBAR_LIST& arRbarUnit, BOOL bRbarDetail=FALSE);
	// General단면의 중립축의 위치,각도가 정해졌을때의 단면의 강도 계산
	// (직사각형등가응력분포를 단면좌표정보를 가진 CPolyMaker를 이용하여 계산)
	// ※ 단면정보를 중립축의 회전각의 반대방향으로 회전시킴으로서 중립축은 수평으로 만들어 부재강도를 계산하며, 
	//    계산된 강도는 회전각(dInRota)에 대해 다시 변환하여 실제축의 강도값으로 변환하여 돌려주는 방식임
	BOOL   Calc_PnMn_Gen_RectStress(double dXn, double dInRota, double dAsThk, T_PM_SECT_RES& ResData, 
		                           double dYmax, double dYmin, double dZmax, double dZmin, CPolyMaker* pPolyMaker, _UMD_RC_RBAR_LIST& arRbarUnit, BOOL bRbarDetail=FALSE);

	// PM상관도로부터 최대 압축 설계강도(αΦPnmax)를 산출함
	double Get_PhiPnmax(T_PMCV_2D& PmData);

	// PmcvData를 보유할 갯수를 넘겨줌 //참고문서 : 신봉진/RC_Column_MPhiCurve자료.ppt (2Page) 
  int Get_HoldNumberToPmcvData();

	// m_RbarData에 저장된 설계용 가상철근폭의 총값을 넘겨줌
	double Get_RbarLenDgn();
	// m_RbarData에 저장된 설계용 가상철근중 고정단면적의 총값을 넘겨줌
	double Get_RbarFixDgnArea();

	// DB단면의 중립축의 위치,각도가 정해졌을때의 단면의 강도 계산
	//   dXn, dInRota : 중립축의 위치비(0.0~1.0, '-'순수압축), 중립축의 회전각(Radian)
	//   dAsThk       : Design상태(필요철근량사정)일때 철근띠의 두께
	// < ResData      : 단면의 강도정보
	BOOL   Calc_PnMn(double dXn, double dInRota, double dAsThk, T_PM_SECT_RES& ResData, BOOL bRbarDetail=FALSE);
	BOOL   Calc_PnMn_Gen(double dXn, double dInRota, double dAsThk, T_PM_SECT_RES& ResData, BOOL bRbarDetail=FALSE);

		// 단면이 순수 압축을 받는지 여부
  //  dXn : 중립축의 직각방향으로의 위치(0~1의 값을 가짐, 0:최하단, 1:최상단)
	BOOL Is_PureComp(double dXn);
  // 단면이 순수 인장을 받는지 여부
  //  dXn : 중립축의 직각방향으로의 위치(0~1의 값을 가짐, 0:최하단, 1:최상단)
	BOOL Is_PureTens(double dXn);
	// E.T.C ///////////////////////////////////
	// 해당 변형률과 중립축 회전각 상태에서의 각 철근별 상세 정보를 계산하여 넘겨줌
	// ※ 주의 : 실제 강도계산과는 독립적으로 이루어 지므로 재정의된 함수가 내부에 들어갈 경우 값이 일치하지 않을 수 있으므로 주의 요망
	BOOL Get_Rebar_DatilData(double dXn, double dRotation, CArray<T_PM_RBAR_RES, T_PM_RBAR_RES&>& arRbarRes, double dAsThk);

		// 중립축의 위치비(dXn)와 회전각을 받아 최 외각 인장 철근까지의 거리를 넘겨줌 
	//   bFixDgnRebarOnly : Dgn일 경우 고정단면적 철근만을 고려한 거리를 넘겨줄지 여부(Dgn이고 가상철근두께 t가 0.0일때에만 적용함. 고정단면적철근 적용시에만 사용함)
	double Calc_MaxDistTensColmRbar(double dXn, double dRota, double dDmax, BOOL bFixDgnRebarOnly = FALSE);
	// dDy, dDz와 회전각을 직접 받아 최 외각 인장 철근까지의 거리를 넘겨줌
	//   bFixDgnRebarOnly : Dgn일 경우 고정단면적 철근만을 고려한 거리를 넘겨줄지 여부(Dgn이고 가상철근두께 t가 0.0일때에만 적용함. 고정단면적철근 적용시에만 사용함)
	double Calc_MaxDistTensColmRbar2(double dDy, double dDz, double dRota, BOOL bFixDgnRebarOnly = FALSE);
	// ConcUnit의 압축여부와 압축구간의 면적, 압축구간의 도심점이 ConcUnit의 도심점과의 차이를 계산함
	// ※관련 문서(DgnCalcBase.dll-Manual.doc참조)
	BOOL Is_ConcComp(int iIndexC, double dRota, _UMD_RC_CONC_UNIT& ConcUnit, double dDy2, double dDz2, double& dCalcAc, double& dyInc, double& dzInc);

	double Calc_esi(double dDcomp, double dDist);// 중립축에서 떨어진거리(dDist)에서의 단면의 변형률
	double Calc_Fsi(BOOL bComp, double desi);    // 중립축에서 떨어진거리(dDist)에서의 철근의 응력
	double Get_ecc();
	double Get_ecu() { return 0.003; } // 콘크리트의 극한 변형율
	double Get_esu() { return 0.0; } // 철근의 인장지배 변형률 한계
	// 철근의 항복시의 변형률
  double Get_ery();

	BOOL Calc_DistanceFromSectionToAngleLine(double dx, double dy, double dAngle, double& dDistMax, double& dDistMin);
	BOOL Cal_MaxDistRbar_Gen(double dXn, double dDmax, double dZmin, _UMD_RC_RBAR_LIST& arRbarUnit, double& dMaxTens, double& dMaxComp, BOOL bFixDgnRebarOnly);
	BOOL Cal_MaxDistRbar_Gen(double dXn, int iAxisDir, double dDmax, double& dMaxTens, double& dMaxComp, BOOL bFixDgnRebarOnly);
	double Get_Xn(int iCurDiv, int iDivNum, double dXb);
	BOOL Cng_AxisDir_Gen(int iAxisDir, _UMD_RC_CON_PART_LIST& arConUnit, _UMD_RC_RBAR_LIST& arRbarUnit);
	double Get_Db();
	double Get_RotationAngle(int iAxisDir);

protected:
	double Get_Sigcc();
	double Calc_SigcU(double dDcomp, double dDist, double dMaxDistCompRbar, double& dec);
	double Calc_SigrU(double dDcomp, double dDist, double dMaxDistCompRbar, double& der);
	double GetMaxRebarRatio(double dMaxRhoUser) { return 0.08; } //!!! code dependent.
	double GetMinRebarRatio()                   { return 0.01; }  //!!! code dependent.
	double Calc_Phi(double dPn, double dPb, double desi); // Pn에서 강도 감소계수φ
	// φPn과 φPb를 이용하여 φPn에서 강도 감소계수φ산출 (현재(07.11.05) 미리 계산된 3DPm-curve정보로 CDgnCalcBase_3DPM_Tool에서 값을 계산하였을 경우 φ를 이미 계산한 좌표망으로 계산하였을 경우 φ값을 재산출할때 사용)
	double Calc_Phi2(double dPhiPn, double dPhiPb, double desi) { return 0.0; }
	double Get_MaxPnReductionFactor() { return 0.8; } //!!! code dependent.
	// 균형파괴상태일때의 중립축의 위치
	double Get_BalancedXb(double dDmax, double dDeff);
	double Get_ReBarStrength(double dBarDia=0.0) { return m_dFyr; } //주철근의 항복응력
	double Calc_Alpha() { return 0.85; } 
	double Calc_Beta();/* { return 0.0; } */

public:
	BOOL SetCommInputData();
	BOOL SetSectInputData(int nSectShap, double dSize[8], const _UMD_RC_COL_MAINRBAR &MainRbar);
	BOOL SetSectInputDataGen(int nSectShap, double dSize[8], _UMD_RC_CON_POLY_DATA &ConPoly, _UMD_RC_COL_MAINRBAR &MainRbar);
	BOOL SetMatlInputData(double dFck, double dEc, double dFyr, double dEs);
	BOOL SetLoadInputData(double dPu, double dMuy, double dMuz, double dMu);
	
protected:
	//--------------------------------------------------------------------------------------------------
	// Input variables
	//--------------------------------------------------------------------------------------------------
	CArray<_UMD_RC_FORCE, _UMD_RC_FORCE&> m_arLoad;      //부재에 작용하는 설계부재력
	int    m_iPmDivNum;                    // Pm상관도 생성시 좌표의 수(생성시에만 적용)
	_UMD_RC_CONC  m_ConcData;    // 콘크리트의 분할된 단면정보	
	_UMD_RC_RBAR  m_RbarData;    // 배근된 주철근정보(Conc좌측하단 원점)
	_UMD_RC_CON_POLY_DATA   m_ConPolyData; // 콘크리트 좌표정보(Conc좌측하단 원점 : 회전되지 않은)	
	BOOL m_bDesign;
	double m_dArea;
	double m_dYbar;   // 도심위치의 Y방향 (좌측하단기준)
	double m_dZbar;   // 도심위치의 Z방향 (좌측하단기준)
	int m_iSectShp;   // Shape ID
	double m_dSize[8];// 크기정보(0:H, 1:B, ...)
	double m_dMaxRhoUser;
	BOOL   m_bEqSpecial; //지진특별기준 적용여부 
	int    m_iSeismicTypeLcom;
	double m_dPu, m_dMu, m_dVu;
	double m_dMuy, m_dMuz, m_dVuy, m_dVuz;
	int  m_iSymmetryType; // Section의 대칭 형태(1:비대칭  2:상하,좌우 대칭  3:z축에 대칭(좌우대칭)  4:y축에 대칭(상하대칭)  5:원점에 대칭) ※관련 문서(신봉진\RC_Column_MPhiCurve자료.ppt 2Page)	

	BOOL   m_bReduceRbar;       // 합성단면 계산시 단면공제 여부(FALSE=공제안함, TRUE=공제함). // 기둥의 순수 압축에서는 무조건 공제(07.10.30)
	int    m_iTypeStressStrain; // 1=Parabolic-Descending, 2=Parabolic-Plateau, 3=Whitney(직사각형 응력블럭), 4=Whitney(직사각형 응력블럭:SliceMethod) //현재(07.10.30) 기둥만 적용함
	int    m_iSliceNum;                    // Slice 개수 
	int    m_iDivisionNumber_90deg;        // Pm상관도 생성시 90의 분할수(생성시에만 적용)
	double m_dTolAngle;                    // Tolerance for Angle (정밀 PM-Curve 산출시 수렴오차 default=0.1˚)
	double m_dTolM0;                       // 최대 모멘트 산출시의 인접 중립축위치 산출및 모멘트 오차비(모멘트 오차 = m_dTolM0*최대모멘트). (default=0.001)  
	double m_dTolP0;                       // 축하중 오차비(축하중 오차 = m_dTolP0*최대축하중). (default=0.001)  
	int    m_iIterNum;
	BOOL   m_bFindMmax;                    // 균형파괴시를 Mmax지점(상세계산)을 사용할지 Get_BalancedXb로 계산된 중립축위치를 사용할지 여부(DB단면은 PM좌표중 최대값(간략계산)을 기준으로 Phi값 산출, Gen단면은 m_bFindMmax에 의해 결정된 균형강도값으로 Phi값 산출)
	
	int m_iDgnCode;   // Code ID
	
	double m_dFc;     // 콘크리트의 압축강도(fck)
	double m_dFyr;    // Main Rebar의 항복응력
	double m_dEsr;
	double m_dEc;
	int m_iHoopType;
	double m_dPhi[5];

	
	//--------------------------------------------------------------------------------------------------
	// Output variables
	//--------------------------------------------------------------------------------------------------
	CArray<T_PM_CHK_RES, T_PM_CHK_RES&> m_arPmmRes;
	T_PMCV_3D m_3DPmData;
	_UMD_RC_COL_ROTATE_SECT m_RotSectData; // 회전된 단면(콘크리트,철근)정보를 보관
		
	//GSD_CALC_FIMP_PROP m_NonLinearProp;  for GSD.

	CString GetFileName();

};

#endif // !defined(AFX_CALCPMCV_H__B8D92DCA_7770_4AA9_B39F_989A6821BE48__INCLUDED_)
