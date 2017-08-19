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
	double m_dTolAngle;// Tolerance for Angle (���� PM-Curve ����� ���ſ��� default=0.1��)
	double m_dTolM0;   // �ִ� ���Ʈ ������� ���� �߸�����ġ ����� ���Ʈ ������(���Ʈ ���� = m_dTolM0*�ִ���Ʈ). (default=0.001)  
	double m_dTolP0;   // ������ ������(������ ���� = m_dTolP0*�ִ�������). (default=0.001)  

	void Set_Tolerance(double dTolAngle, double dTolM0, double dTolP0);
	
	// ������ ����(T_PMCV_3D)�� CDgnCalcBase_3DPM_Tool�� ����(_UMD_3DPM_DATA)�� ��ȯ�Ͽ���
	//   nType : Phi������� (1:Phi�̰��, 2:Phi���)
	BOOL Get_ConvertToPMTool(int nType, T_PMCV_3D& InData, _UMD_3DPM_DATA& OutData);

	// 3D-PM��ü ����� �̿��Ͽ� �ܷ°� ���Ϲ���(Ps,dMs_y,dMs_z ����)������ PM-Curve���� ���� ã�Ƽ� �Ѱ���
	//   Umd3DPmData     : 3DPM ����� ����
	//   dPs, dMs_y, dMs_z : �ܷ¹����� �������� ������, y�࿡ ���� ���Ʈ, z�࿡ ���� ���Ʈ(Umd3DPmData�� ������ ���� ���ο� �־����)
  //   dPd, dMd_y, dMd_z : �ܷ��� ������, y�࿡ ���� ���Ʈ, z�࿡ ���� ���Ʈ
  // < dPn, dMn_y, dMn_z, desi, dXn : �ܷ°� ������ ������ PM-Curve���� ������, y�࿡ ���� ���Ʈ, z�࿡ ���� ���Ʈ�� �Ѱ� ���� ����
  // < bCheck            : �ܷ¿� �ߵ��� �ִ��� ���θ� �Ѱ� ���� ����
	// < iLeftFitDirID, dModul_DirID : �������� ��ġ�� �˷��ִ� Umd3DPmData�� ID�� ������� �Ѱ� ���� ����
	//   bChkInStartPMM : ������(Ps,dMs_y,dMs_z)�� ���ο� �ִ����� üũ���� ����	
	// < bResInStartPMM : ������(Ps,dMs_y,dMs_z)�� ���ο� �ִ����� üũ���	(bChkInStartPMM=TRUE�϶��� ������)
	// �������� : �ź���/RC_Column_MPhiCurve�ڷ�.ppt (3Page) 
	// �����ǻ��� : �� Ȯ���� ��!!!
	//              -����, ������, �ۿ�ܷ��� MyMz����� ���� ������ ���� ���� ��� �� ���� 3D-PM(�հ��� �����)�� ����Ͽ��� �Ѵ�.
	//              -������(Ps,dMs_y,dMs_z)�� 3DPM ���ο� �־�߸� �ٸ� ����� ����, 
	//               �������� �ܺο� �������� �� ��쿡�� �ݵ�� bChkInStartPMM�� TRUE�� �ξ� Check�� �����Ͽ��� ��
  BOOL Get_PmcvUnitAt3DPmcvData(_UMD_3DPM_DATA& Umd3DPmData, double dPs, double dMs_y, double dMs_z, 
		double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, 
		BOOL& bCheck, int& iLeftFitDirID, double& dModul_DirID, BOOL bChkInStartPMM, BOOL& bResInStartPMM);
	// �̹� ���� M-����ü ����� �̿��Ͽ� �ܷ°� ���Ϲ���(��������)������ PM-Curve���� ���� ã�Ƽ� �Ѱ���
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

	// �ش� �ܷ¹��⿡ ���� Line�� arPmcvUnit���� �Ҷմ��� ���� �Ѱ���
	//   dPs, dMs_y, dMs_z : �ܷ¹����� �������� ������, y�࿡ ���� ���Ʈ, z�࿡ ���� ���Ʈ
  //   dPd, dMd_y, dMd_z : �ܷ��� ������, y�࿡ ���� ���Ʈ, z�࿡ ���� ���Ʈ
  // < arPmcvUnit: �ش� �ܷ¹��⿡ �ִ��� Ȯ���� 3D��ǥ���� �����ϴ� �迭
  // < RETURN: ��������
  BOOL Cal_PierceCheckToPolyLine(double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit);
	// �������� ���� �ش� �ܷ¹��⿡ ���� Line�� arPmcvUnit���� �Ҷմ��� ���� �Ѱ���
	BOOL Cal_PierceCheckToPolyLine(double dPd, double dMd_y, double dMd_z, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit);

	// �ش� �ܷ¹��⿡ ���� Line�� 3������ arPmcvUnit�� ��Ÿ���� ������ �������� ����Ͽ� �Ѱ���(���⼺�� ����)
	//   dPs, dMs_y, dMs_z : �ܷ¹����� �������� ������, y�࿡ ���� ���Ʈ, z�࿡ ���� ���Ʈ
  //   dPd, dMd_y, dMd_z : �ܷ��� ������, y�࿡ ���� ���Ʈ, z�࿡ ���� ���Ʈ
  // < dPn, dMn_y, dMn_z : �ܷ°� ������ ������ ������ �������� ���� ������, y�࿡ ���� ���Ʈ, z�࿡ ���� ���Ʈ�� �Ѱ� ���� ����
	// < desi, dXn         : Pn,Mn�ÿ��� �ִ� ����ö���� ������ (����� '-'), �߸����� �������������� ��ġ(�ϴܿ��������� ����) 
  // < arPmcvUnit: �ش� �ܷ¹���� �������� ã�� 3D��ǥ���� �����ϴ� 3���� �迭
  // < RETURN: �������� ã���� ���ų� �������� �ܷ¹��⿡ ���� ��� FALSE�� �ѱ� 
  BOOL Cal_PiercePointToPolyLine(double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit);
	// �������� ���� �ش� �ܷ¹��⿡ ���� Line�� 3������ arPmcvUnit�� ��Ÿ���� ������ �������� ����Ͽ� �Ѱ���(���⼺�� ����)
	BOOL Cal_PiercePointToPolyLine(double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit);
	
	// �ش� �ܷ¹��⿡ ���� Line�� 3������ arPmcvUnit�� ��Ÿ���� ������ �������� ����Ͽ� �Ѱ���(Line�� ��������)
  //   dPs, dMs_y, dMs_z : �ܷ¹����� �������� ������, y�࿡ ���� ���Ʈ, z�࿡ ���� ���Ʈ
	//   dPd, dMd_y, dMd_z : �ܷ��� ������, y�࿡ ���� ���Ʈ, z�࿡ ���� ���Ʈ
  // < dPn, dMn_y, dMn_z : �ܷ°� ������ ������ ������ �������� ���� ������, y�࿡ ���� ���Ʈ, z�࿡ ���� ���Ʈ�� �Ѱ� ���� ����
	// < desi, dXn         : Pn,Mn�ÿ��� �ִ� ����ö���� ������ (����� '-'), �߸����� �������������� ��ġ(�ϴܿ��������� ����) 
  // < arPmcvUnit: �ش� �ܷ¹���� �������� ã�� 3D��ǥ���� �����ϴ� 3���� �迭
	BOOL Cal_CrossPointToPolyLine(double dPs, double dMs_y, double dMs_z, double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit);
	// �������� �ش� �ܷ¹��⿡ ���� Line�� 3������ arPmcvUnit�� ��Ÿ���� ������ �������� ����Ͽ� �Ѱ���(Line�� ��������)
  BOOL Cal_CrossPointToPolyLine(double dPd, double dMd_y, double dMd_z, double& dPn, double& dMn_y, double& dMn_z, double& desi, double& dXn, CArray<_UMD_PMCV_UNIT_POLY, _UMD_PMCV_UNIT_POLY&>& arPmcvUnit);

	BOOL Get_UnitVectorToMzMyP(double dP, double dMy, double dMz, double vectorMzMyP[3], double vectorMP[2], double vectorMzMy[2]);

	void Get_AngleMinIDToPmCurve_NonScale(_UMD_3DPM_DATA& Umd3DPmData, double nMzMyP_Org[3], double nVectorMzMyP[3], int iAxisDir_Global, double dPmax, double dMmax, int& iUnitID, double& dAngleMin);

	BOOL Cal_MathPierceCheckToPolyLine(double p0[3], double p1[3], const int nData, double polyLine[][3]);

	BOOL Cal_MathPierceCheckToPolyLine(double p1[3], const int nData, double polyLine[][3]);


};

#endif // !defined(AFX_MDGNPMTOOL_H__A8A50FB1_A2F6_46A6_B91E_B24EBF392B8D__INCLUDED_)
