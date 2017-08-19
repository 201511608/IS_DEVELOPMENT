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
	// �߸����� ȸ������ ���� PM������� ����Ͽ���
  BOOL Get_PMCurve(double dInRota, T_PMCV_2D& PmcvData, T_PM_SECT_RES& CbPmData, double& dPhiPnmax, double dAsThk);
	// General�ܸ��� �߸����� ȸ������ ���� PM������� ����Ͽ���
	// ������ : 1ȸ ���ÿ��� ��� (�ݺ����ÿ��� �ӵ��� ���ϵ�)
	BOOL Get_PMCurve_Gen(double dInRota, T_PMCV_2D& PmcvData, T_PM_SECT_RES& CbPmData, double& dPhiPnmax, double dAsThk, BOOL bRbarDetail=FALSE);
	// �߸����� ȸ������ ���� PM������� ����Ͽ���
  BOOL Get_PMCurve_Gen(int iAxisDir, T_PMCV_2D& PmcvData, T_PM_SECT_RES& CbPmData, T_PM_SECT_RES& P0PmData, double& dPhiPnmax, double dAsThk, BOOL bRbarDetail=FALSE);
	BOOL Check_PnMn(double dAsCur, T_PMCV_2D& PmcvData, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, double& dRatio, T_PM_SECT_RES* pResData=NULL);
	BOOL Check_PnMn_Gen(_UMD_3DPM_DATA& PmData3D, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, double& dRatio, T_PMCV_2D& PmcvData, T_PM_SECT_RES& CbPmData, T_PM_SECT_RES* pResData=NULL);

	// DB�ܸ�� ETC �Լ� ============================================================================================
	double Calc_Rota(double dInRota=-1.0);// dInRota(Radian)
	double Calc_Dmax(double dInRota=-1.0);// dInRota(Radian)

	// PM������κ��� �ۿ���� ������ ������
	double Calc_Ratio(double dPu, double dMuy, double dMuz, double dPhiPnmax, T_PMCV_2D& PmData, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, BOOL bPureMonemt=FALSE);
	double Calc_Ratio(double dPu, double dMuy, double dMuz, double dPhiPnmax, T_PMCV_2D& PmData, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, T_PM_SECT_RES& ResData, BOOL bPureMonemt=FALSE);
	double Calc_Ratio_ToPu(double dPu, double dMuy, double dMuz, double dPhiPnmax, T_PMCV_2D& PmData, double& dPhi, double& dPhiPn, double& dPhiMny, double& dPhiMnz, T_PM_SECT_RES& ResData);

	// �ش� ȸ������ ���� ȸ���� �ܸ���ǥ�� ��ȿ���̸� �Ѱ���(�ܸ���ǥ�� ȸ���߽����� �߽����� �߸��� ȸ������ ���������� ȸ���� ���� ����) 
  //   dDmax     : �ܸ��� �߸��� �������������� ���� 
  //   dZmim     : �ش� ȸ������ ���� ȸ���� �ܸ���ǥ�� �ּ� Z��(�ܸ���ǥ�� ȸ���߽����� �߽����� �߸��� ȸ������ ���������� ȸ���� ���� ����) 
  //   arRbarUnit: ö���� ��ġ����(�ߺ��� ���ϱ����� ������ �Ѱ���:��������)
	//   bFixDgnRebarOnly : Dgn�� ��� �����ܸ��� ö�ٸ��� ����� �Ÿ��� �Ѱ����� ����(Dgn�̰� ����ö�ٵβ� t�� 0.0�϶����� ������. �����ܸ���ö�� ����ÿ��� �����)
	double Calc_Deff_Gen(double dDmax, double dZmin, _UMD_RC_RBAR_LIST& arRbarUnit, BOOL bFixDgnRebarOnly = FALSE);
  // ��ȿ���� ����ؼ� �Ѱ���. 4��( 0=Y-Dir-Positive, 1=Y-Dir-Negative, 2=Z-Dir-Positive, 3=Z-Dir-Negative)
  double Calc_Deff_Gen(int iDir, BOOL bFixDgnRebarOnly = FALSE);    

	// General�ܸ��� �߸����� ��ġ,������ ������������ �ܸ��� ���� ���(��ü ���� : ���� �������ݿ� ���� �ݺ�����)
	// (���� ���º����� �߸��࿡ �����ϰ� Sliceģ �ܸ� ������ �̿��Ͽ� ���)
	//   dXn, dInRota : �߸����� ��ġ��(0.0~1.0, 0.0��������), �߸����� ȸ����(���� ������ �������� ������ ��ȯ�ϱ� ���ؼ��� �����)
	//   dAsThk       : Design����(�ʿ�ö�ٷ�����)�϶� ö�ٶ��� �β�
	// < ResData      : �ܸ��� ��������
	//   arConcUnit   : �߸����� ȸ������ �ݴ�������� ȸ���� Sliceģ �ܸ�����
	//   arRbarUnit   : �߸����� ȸ������ �ݴ�������� ȸ���� ö������
	// �� �ܸ������� �߸����� ȸ������ �ݴ�������� ȸ����Ŵ���μ� �߸����� �������� ����� ���簭���� ����ϸ�, 
	//    ���� ������ ȸ����(dInRota)�� ���� �ٽ� ��ȯ�Ͽ� �������� ���������� ��ȯ�Ͽ� �����ִ� �����
	BOOL   Calc_PnMn_Gen_SliceMethod(double dXn, double dInRota, double dAsThk, T_PM_SECT_RES& ResData, 
		_UMD_RC_CON_PART_LIST& arConcUnit, _UMD_RC_RBAR_LIST& arRbarUnit, BOOL bRbarDetail=FALSE);
	// General�ܸ��� �߸����� ��ġ,������ ������������ �ܸ��� ���� ���
	// (���簢������º����� �ܸ���ǥ������ ���� CPolyMaker�� �̿��Ͽ� ���)
	// �� �ܸ������� �߸����� ȸ������ �ݴ�������� ȸ����Ŵ���μ� �߸����� �������� ����� ���簭���� ����ϸ�, 
	//    ���� ������ ȸ����(dInRota)�� ���� �ٽ� ��ȯ�Ͽ� �������� ���������� ��ȯ�Ͽ� �����ִ� �����
	BOOL   Calc_PnMn_Gen_RectStress(double dXn, double dInRota, double dAsThk, T_PM_SECT_RES& ResData, 
		                           double dYmax, double dYmin, double dZmax, double dZmin, CPolyMaker* pPolyMaker, _UMD_RC_RBAR_LIST& arRbarUnit, BOOL bRbarDetail=FALSE);

	// PM������κ��� �ִ� ���� ���谭��(���Pnmax)�� ������
	double Get_PhiPnmax(T_PMCV_2D& PmData);

	// PmcvData�� ������ ������ �Ѱ��� //������ : �ź���/RC_Column_MPhiCurve�ڷ�.ppt (2Page) 
  int Get_HoldNumberToPmcvData();

	// m_RbarData�� ����� ����� ����ö������ �Ѱ��� �Ѱ���
	double Get_RbarLenDgn();
	// m_RbarData�� ����� ����� ����ö���� �����ܸ����� �Ѱ��� �Ѱ���
	double Get_RbarFixDgnArea();

	// DB�ܸ��� �߸����� ��ġ,������ ������������ �ܸ��� ���� ���
	//   dXn, dInRota : �߸����� ��ġ��(0.0~1.0, '-'��������), �߸����� ȸ����(Radian)
	//   dAsThk       : Design����(�ʿ�ö�ٷ�����)�϶� ö�ٶ��� �β�
	// < ResData      : �ܸ��� ��������
	BOOL   Calc_PnMn(double dXn, double dInRota, double dAsThk, T_PM_SECT_RES& ResData, BOOL bRbarDetail=FALSE);
	BOOL   Calc_PnMn_Gen(double dXn, double dInRota, double dAsThk, T_PM_SECT_RES& ResData, BOOL bRbarDetail=FALSE);

		// �ܸ��� ���� ������ �޴��� ����
  //  dXn : �߸����� �������������� ��ġ(0~1�� ���� ����, 0:���ϴ�, 1:�ֻ��)
	BOOL Is_PureComp(double dXn);
  // �ܸ��� ���� ������ �޴��� ����
  //  dXn : �߸����� �������������� ��ġ(0~1�� ���� ����, 0:���ϴ�, 1:�ֻ��)
	BOOL Is_PureTens(double dXn);
	// E.T.C ///////////////////////////////////
	// �ش� �������� �߸��� ȸ���� ���¿����� �� ö�ٺ� �� ������ ����Ͽ� �Ѱ���
	// �� ���� : ���� ���������� ���������� �̷�� ���Ƿ� �����ǵ� �Լ��� ���ο� �� ��� ���� ��ġ���� ���� �� �����Ƿ� ���� ���
	BOOL Get_Rebar_DatilData(double dXn, double dRotation, CArray<T_PM_RBAR_RES, T_PM_RBAR_RES&>& arRbarRes, double dAsThk);

		// �߸����� ��ġ��(dXn)�� ȸ������ �޾� �� �ܰ� ���� ö�ٱ����� �Ÿ��� �Ѱ��� 
	//   bFixDgnRebarOnly : Dgn�� ��� �����ܸ��� ö�ٸ��� ����� �Ÿ��� �Ѱ����� ����(Dgn�̰� ����ö�ٵβ� t�� 0.0�϶����� ������. �����ܸ���ö�� ����ÿ��� �����)
	double Calc_MaxDistTensColmRbar(double dXn, double dRota, double dDmax, BOOL bFixDgnRebarOnly = FALSE);
	// dDy, dDz�� ȸ������ ���� �޾� �� �ܰ� ���� ö�ٱ����� �Ÿ��� �Ѱ���
	//   bFixDgnRebarOnly : Dgn�� ��� �����ܸ��� ö�ٸ��� ����� �Ÿ��� �Ѱ����� ����(Dgn�̰� ����ö�ٵβ� t�� 0.0�϶����� ������. �����ܸ���ö�� ����ÿ��� �����)
	double Calc_MaxDistTensColmRbar2(double dDy, double dDz, double dRota, BOOL bFixDgnRebarOnly = FALSE);
	// ConcUnit�� ���࿩�ο� ���౸���� ����, ���౸���� �������� ConcUnit�� ���������� ���̸� �����
	// �ذ��� ����(DgnCalcBase.dll-Manual.doc����)
	BOOL Is_ConcComp(int iIndexC, double dRota, _UMD_RC_CONC_UNIT& ConcUnit, double dDy2, double dDz2, double& dCalcAc, double& dyInc, double& dzInc);

	double Calc_esi(double dDcomp, double dDist);// �߸��࿡�� �������Ÿ�(dDist)������ �ܸ��� ������
	double Calc_Fsi(BOOL bComp, double desi);    // �߸��࿡�� �������Ÿ�(dDist)������ ö���� ����
	double Get_ecc();
	double Get_ecu() { return 0.003; } // ��ũ��Ʈ�� ���� ������
	double Get_esu() { return 0.0; } // ö���� �������� ������ �Ѱ�
	// ö���� �׺����� ������
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
	double Calc_Phi(double dPn, double dPb, double desi); // Pn���� ���� ���Ұ����
	// ��Pn�� ��Pb�� �̿��Ͽ� ��Pn���� ���� ���Ұ�������� (����(07.11.05) �̸� ���� 3DPm-curve������ CDgnCalcBase_3DPM_Tool���� ���� ����Ͽ��� ��� ���� �̹� ����� ��ǥ������ ����Ͽ��� ��� ������ ������Ҷ� ���)
	double Calc_Phi2(double dPhiPn, double dPhiPb, double desi) { return 0.0; }
	double Get_MaxPnReductionFactor() { return 0.8; } //!!! code dependent.
	// �����ı������϶��� �߸����� ��ġ
	double Get_BalancedXb(double dDmax, double dDeff);
	double Get_ReBarStrength(double dBarDia=0.0) { return m_dFyr; } //��ö���� �׺�����
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
	CArray<_UMD_RC_FORCE, _UMD_RC_FORCE&> m_arLoad;      //���翡 �ۿ��ϴ� ��������
	int    m_iPmDivNum;                    // Pm����� ������ ��ǥ�� ��(�����ÿ��� ����)
	_UMD_RC_CONC  m_ConcData;    // ��ũ��Ʈ�� ���ҵ� �ܸ�����	
	_UMD_RC_RBAR  m_RbarData;    // ��ٵ� ��ö������(Conc�����ϴ� ����)
	_UMD_RC_CON_POLY_DATA   m_ConPolyData; // ��ũ��Ʈ ��ǥ����(Conc�����ϴ� ���� : ȸ������ ����)	
	BOOL m_bDesign;
	double m_dArea;
	double m_dYbar;   // ������ġ�� Y���� (�����ϴܱ���)
	double m_dZbar;   // ������ġ�� Z���� (�����ϴܱ���)
	int m_iSectShp;   // Shape ID
	double m_dSize[8];// ũ������(0:H, 1:B, ...)
	double m_dMaxRhoUser;
	BOOL   m_bEqSpecial; //����Ư������ ���뿩�� 
	int    m_iSeismicTypeLcom;
	double m_dPu, m_dMu, m_dVu;
	double m_dMuy, m_dMuz, m_dVuy, m_dVuz;
	int  m_iSymmetryType; // Section�� ��Ī ����(1:���Ī  2:����,�¿� ��Ī  3:z�࿡ ��Ī(�¿��Ī)  4:y�࿡ ��Ī(���ϴ�Ī)  5:������ ��Ī) �ذ��� ����(�ź���\RC_Column_MPhiCurve�ڷ�.ppt 2Page)	

	BOOL   m_bReduceRbar;       // �ռ��ܸ� ���� �ܸ���� ����(FALSE=��������, TRUE=������). // ����� ���� ���࿡���� ������ ����(07.10.30)
	int    m_iTypeStressStrain; // 1=Parabolic-Descending, 2=Parabolic-Plateau, 3=Whitney(���簢�� ���º�), 4=Whitney(���簢�� ���º�:SliceMethod) //����(07.10.30) ��ո� ������
	int    m_iSliceNum;                    // Slice ���� 
	int    m_iDivisionNumber_90deg;        // Pm����� ������ 90�� ���Ҽ�(�����ÿ��� ����)
	double m_dTolAngle;                    // Tolerance for Angle (���� PM-Curve ����� ���ſ��� default=0.1��)
	double m_dTolM0;                       // �ִ� ���Ʈ ������� ���� �߸�����ġ ����� ���Ʈ ������(���Ʈ ���� = m_dTolM0*�ִ���Ʈ). (default=0.001)  
	double m_dTolP0;                       // ������ ������(������ ���� = m_dTolP0*�ִ�������). (default=0.001)  
	int    m_iIterNum;
	BOOL   m_bFindMmax;                    // �����ı��ø� Mmax����(�󼼰��)�� ������� Get_BalancedXb�� ���� �߸�����ġ�� ������� ����(DB�ܸ��� PM��ǥ�� �ִ밪(�������)�� �������� Phi�� ����, Gen�ܸ��� m_bFindMmax�� ���� ������ �������������� Phi�� ����)
	
	int m_iDgnCode;   // Code ID
	
	double m_dFc;     // ��ũ��Ʈ�� ���భ��(fck)
	double m_dFyr;    // Main Rebar�� �׺�����
	double m_dEsr;
	double m_dEc;
	int m_iHoopType;
	double m_dPhi[5];

	
	//--------------------------------------------------------------------------------------------------
	// Output variables
	//--------------------------------------------------------------------------------------------------
	CArray<T_PM_CHK_RES, T_PM_CHK_RES&> m_arPmmRes;
	T_PMCV_3D m_3DPmData;
	_UMD_RC_COL_ROTATE_SECT m_RotSectData; // ȸ���� �ܸ�(��ũ��Ʈ,ö��)������ ����
		
	//GSD_CALC_FIMP_PROP m_NonLinearProp;  for GSD.

	CString GetFileName();

};

#endif // !defined(AFX_CALCPMCV_H__B8D92DCA_7770_4AA9_B39F_989A6821BE48__INCLUDED_)
