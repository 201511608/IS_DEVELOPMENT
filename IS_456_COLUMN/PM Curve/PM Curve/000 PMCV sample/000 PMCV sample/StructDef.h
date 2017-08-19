#if !defined(__DGNCALC_KR_RCSC_COLUMNSTRUCT_H__)
#define __DGNCALC_KR_RCSC_COLUMNSTRUCT_H__

#include <afxtempl.h>

const double cUMDRC_Zero	         = 1.0E-07;
const double cUMDRC_Upon	         = 1.0E+28;

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
#define DGN_SECT_SHAPE_INDEX_REG_GEN   21  // General Section (Value������ ���)
#define DGN_SECT_SHAPE_INDEX_REG_CR    22  // Cross (+�� �ܸ�) // Only Building
#define DGN_SECT_SHAPE_INDEX_REG_UDT   23  // Upside-Down T
#define DGN_SECT_SHAPE_INDEX_REG_CB    24

const double CONST_UMDRC_dRATLIM2  = 0.025;	// Limit for Neutral Axis (Column).
const double CONST_UMDRC_dRATLIM3  = 0.050;	// Limit for Rotated Angle[Deg] (Muy/Muz).
const int    CONST_UMDRC_iRECT		 =  4;
const int    CONST_UMDRC_iBAR_LAY  =  5;
const int    CONST_UMDRC_iBAR_DIV  = 20;
 
enum EN_MDGN_CON_STRESS_STRAIN_TYPE // ��ũ��Ʈ�� ����-������ ���� TYPE.
{
	EN_MDGN_CON_STRESS_STRAIN_PARABOLIC_DESCENDING = 1,
	EN_MDGN_CON_STRESS_STRAIN_PARABOLIC_PLATEAU,
	EN_MDGN_CON_STRESS_STRAIN_WHITNEY,        // Whitney(���簢�� ���º�)
	EN_MDGN_CON_STRESS_STRAIN_WHITNEY_SLICE,  // Whitney(���簢�� ���º�:SliceMethod)
	EN_MDGN_CON_STRESS_STRAIN_NONLINEAR,
	EN_MDGN_CON_STRESS_STRAIN_USER,
};

struct T_PM_UNIT
{
	double dPhi;
	double dXn;
	double dPn;
	double dMny;
  double dMnz;
  double desi;
	void Init() 
	{
		dPhi = 0.0;
		dXn  = 0.0;
		dPn  = 0.0;
		dMny = 0.0;
		dMnz = 0.0;
		desi = 0.0;
	}
	T_PM_UNIT()	{}
	T_PM_UNIT& operator= (const T_PM_UNIT &src)
	{
		dPhi = src.dPhi;
		dXn  = src.dXn;
		dPn  = src.dPn;
		dMny = src.dMny;
		dMnz = src.dMnz;
		desi = src.desi;
		return *this; 
	}
};

struct T_PMCV_2D// _UMD_RC_PM
{
	double dRotate;// �߸��� ȸ����(Radian)
	double dDmax;
	double dDeff;
	T_PM_UNIT Pmemin;
	CArray<T_PM_UNIT,T_PM_UNIT&>	arPmUnit;
	
	void Init() 
	{
		dRotate = 0.0;
		dDmax   = 0.0;
		dDeff   = 0.0;
		Pmemin.Init();
		arPmUnit.RemoveAll();
	}
	T_PMCV_2D()	{}
	T_PMCV_2D(const T_PMCV_2D &src) { *this = src; }
	T_PMCV_2D& operator = (const T_PMCV_2D &src)
	{
		dRotate = src.dRotate;
		dDmax   = src.dDmax;
		dDeff   = src.dDeff;
		Pmemin  = src.Pmemin;
		arPmUnit.RemoveAll();
		arPmUnit.Copy(src.arPmUnit);
		return *this; 
	}
};

struct T_PM_RBAR_RES // _UMD_RC_COL_RBAR_RES
{
	double dyz[2];	// 0:y, 1:z.
	double dAs;
	double dStress;
	double dForce;
  double dds;
  double desi;
	void Init() 
	{
		dyz[0] = dyz[1] = 0.0;
		dAs = 0.0;
		dStress = 0.0;
		dForce = 0.0;
		dds = 0.0;
		desi = 0.0;
	}	
	T_PM_RBAR_RES()	{}
	T_PM_RBAR_RES& operator= (const T_PM_RBAR_RES &src) 
	{
		dyz[0] = src.dyz[0];
		dyz[1] = src.dyz[1];
		dAs = src.dAs;
		dStress = src.dStress;
		dForce  = src.dForce;
		dds  = src.dds;
		desi = src.desi;
		return *this; 
	}
};

struct T_PM_SECT_RES// T_PM_SECT_RES
{
	int    nType;    //Data���� (0:��������  1:��  2:����  3:��+ö�ٻ�)
	// ��������
	double dNARotate;//�߸��� ȸ����(Radian)
  double dCb;
	double dXn;
  double dPn;
	double dMny;
	double dMnz;
	double desiMax;
	// ������
  double dPncc;
  double dPnsc;
  double dPnst;  
  double dMnccy;
  double dMnscy;
  double dMnsty;  
  double dMnccz;
  double dMnscz;
  double dMnstz;	
	double dCForce;
	double dTForce;	
	// ö�ٻ�����
	CArray<T_PM_RBAR_RES, T_PM_RBAR_RES&> arRbarRes;
  void Init()
	{
		nType = 0;
		dNARotate = 0.0;
		dCb = 0.0;
		dXn = 0.0;
		dPn = 0.0;
		dMny= 0.0;
		dMnz= 0.0;
		desiMax = 0.0;
		dPncc = 0.0;
		dPnsc = 0.0;
		dPnst = 0.0;
		dMnccy = 0.0;
		dMnscy = 0.0;
		dMnsty = 0.0;
		dMnccz = 0.0;
		dMnscz = 0.0;
		dMnstz = 0.0;
		dCForce = 0.0;
		dTForce = 0.0;
		arRbarRes.RemoveAll();
	}
  T_PM_SECT_RES()	{}	
	T_PM_SECT_RES(const T_PM_SECT_RES &src) { *this = src; }
	T_PM_SECT_RES& operator = (const T_PM_SECT_RES &src) 
	{
		nType = src.nType;
		dCb = src.dCb;
		dXn = src.dXn;
		dPn = src.dPn;
		dMny= src.dMny;
		dMnz= src.dMnz;
		desiMax = src.desiMax;
		dPncc = src.dPncc;
		dPnsc = src.dPnsc;
		dPnst = src.dPnst;
		dMnccy = src.dMnccy;
		dMnscy = src.dMnscy;
		dMnsty = src.dMnsty;
		dMnccz = src.dMnccz;
		dMnscz = src.dMnscz;
		dMnstz = src.dMnstz;
		dCForce = src.dCForce;
		dTForce = src.dTForce;
		arRbarRes.RemoveAll();
		return *this; 
	}
};



struct T_PM_CHK_RES//_UMD_RC_COL_PMM_RES
{	
	BOOL   bUse     ;//���� ���� ����
	double dPhi     ;
	double dPhiPn   ;
	double dPhiMn   ;
	double dPhiMny  ;
	double dPhiMnz  ;
	double dPhiPnmax;
	double dRatioP  ;
	double dRatioM  ;
	double dAs      ;
	int    iAsRes   ;// 0:ö�ٺ���� ����, 1:�ּ�ö�ٺ� �̸�, 2:�ִ�ö�ٺ� �ʰ�
	double dRho     ;
	double dRhoMin  ;
	double dRhoMax  ;
	BOOL   bPmcvData;
	BOOL   bCbData  ;
	T_PMCV_2D           PmcvData;	
	T_PM_SECT_RES CbPmData;
  T_PM_SECT_RES PuPnData;

	void Init() 
	{
		bUse = FALSE;
		dPhi = 0.0;
		dPhiPn = 0.0;
		dPhiMn = 0.0;
		dPhiMny = 0.0;
		dPhiMnz = 0.0;
		dPhiPnmax = 0.0;
		dRatioP = 0.0;
		dRatioM = 0.0;
		dAs = 0.0;
		iAsRes = 0;
		dRho = 0.0;
		dRhoMin = 0.0;
		dRhoMax = 0.0;
		bPmcvData = FALSE;
		bCbData   = FALSE;
		PmcvData.Init();
		CbPmData.Init();
		PuPnData.Init();
	}
	T_PM_CHK_RES() {}
	T_PM_CHK_RES& operator= (const T_PM_CHK_RES &src)
	{
		bUse = src.bUse;
		dPhi = src.dPhi;
		dPhiPn = src.dPhiPn;
		dPhiMn = src.dPhiMn;
		dPhiMny = src.dPhiMny;
		dPhiMnz = src.dPhiMnz;
		dPhiPnmax = src.dPhiPnmax;
		dRatioP   = src.dRatioP;
		dRatioM   = src.dRatioM;
		dAs = src.dAs;
		iAsRes = src.iAsRes;
		dRho = src.dRho;
		dRhoMax = src.dRhoMax;
		dRhoMin = src.dRhoMin;
		bPmcvData = src.bPmcvData;
		bCbData   = src.bCbData;
		PmcvData  = src.PmcvData;
		CbPmData  = src.CbPmData;
		PuPnData  = src.PuPnData;
		return *this;
	}
};


struct T_PMCV_3D//T_PMCV_3D
{
	int iDivisionNumber_90deg;// 3D-PM����� 90���� ������ ���� (defult=6 : 15���������� ����)  
  int iSymmetryType; // Section�� ��Ī ����(1:���Ī  2:����,�¿� ��Ī  3:z�࿡ ��Ī(�¿��Ī)  4:y�࿡ ��Ī(���ϴ�Ī)  5:������ ��Ī) �ذ��� ����(�ź���\RC_Column_MPhiCurve�ڷ�.ppt 2Page)
	double dPhiPnmax;
	BOOL b3DData;
	BOOL bCbData;
	BOOL bP0Data;
	CArray<T_PMCV_2D,T_PMCV_2D&>	arPmcvData;
	CArray<T_PM_SECT_RES,T_PM_SECT_RES&>	arCbPmData;
	CArray<T_PM_SECT_RES,T_PM_SECT_RES&>	arP0PmData;
	void Init()
	{
		iDivisionNumber_90deg = 0;
		iSymmetryType = 0;
		dPhiPnmax = 0.0;
		b3DData = FALSE;
		bCbData = FALSE;
		bP0Data = FALSE;
		arPmcvData.RemoveAll();
		arCbPmData.RemoveAll();
		arP0PmData.RemoveAll();
	}
	T_PMCV_3D()	{}
	T_PMCV_3D& operator= (const T_PMCV_3D& src) 
	{
		iDivisionNumber_90deg = src.iDivisionNumber_90deg;
		iSymmetryType = src.iSymmetryType;
		dPhiPnmax  = src.dPhiPnmax;
		b3DData = src.b3DData;
		bCbData = src.bCbData;
		bP0Data = src.bP0Data;
		arPmcvData.RemoveAll(); arPmcvData.Copy(src.arPmcvData);
		arCbPmData.RemoveAll(); arCbPmData.Copy(src.arCbPmData);
		arP0PmData.RemoveAll(); arP0PmData.Copy(src.arP0PmData);
		return *this; 
	}
};

// PM-curve���� �ϳ��� ��ǥ������ ������ ����ü
struct _UMD_PMCV_UNIT_POLY
{
	double dPn;  // ������  
  double dMny; // y�࿡ ���� ���Ʈ(��ü��ǥ��)  
  double dMnz; // z�࿡ ���� ���Ʈ(��ü��ǥ��)  
  double dXn;  // �߸����� �������������� ��ġ(�ϴܿ��������� ����)  
  double desi; // �ִ� ����ö���� ������ (����� '-')
  // �ʱ�ȭ
  void Init() {}
	_UMD_PMCV_UNIT_POLY() { Init(); }
	_UMD_PMCV_UNIT_POLY& operator= (const _UMD_PMCV_UNIT_POLY &src)
	{
		dPn = src.dPn;
		dMny = src.dMny;
		dMnz = src.dMnz;
		dXn  = src.dXn;
		desi = src.desi;
		return *this;
	}
};


//ȸ���� �Ҽ��߽��࿡ ���� PM-curve ������ ������ ����ü
struct _UMD_PMCV_DATA_POLY
{ //����(07.11.01) �������� 0�϶��� ���� ������� �ʰ� �ֽ��ϴ�. 
  double dRotate;    // �Ҽ��߽����� ȸ���� �Ҽ��߽����� y������̰� ����ΰ� Z��'+'�����϶� 0���� ���� ������ �ݽð���� (Deg)
	double dXb;        // ���� �ı������϶��� �Ҽ��߽����� �������������� ��ġ(0~1�� ���� ����, 0:���ϴ�, 1:�ֻ��) �Ǵ� ����� ����(Get_PmcvData�Ϸ���)
	double dPb;        // ���� �ı������϶��� ������
  double dMby, dMbz; // ���� �ı������϶��� ���Ʈ(��ü ��ǥ��)
	double desib;      // ���� �ı������϶��� �ִ� ����ö���� ������ (����� '-')
  double dX0;        // �������� 0�϶��� �Ҽ��߽����� �������������� ��ġ �Ǵ� ����� ����(Get_PmcvData�Ϸ���) 
  //double de0;        // �������� 0�϶��� ����� ������
	double dM0y, dM0z; // �������� 0�϶��� ���Ʈ(��ü ��ǥ��)
	double dDmax;      // �ܸ����(�߸����� �������⿡ ���� ��ũ��Ʈ�� �ֻ��,���ϴܱ����� �Ÿ�)
	double dDeff;      // ��ȿ����(�߸����� �������⿡ ���� ���࿬�ܿ��� �ֿܰ� ����ö�ٱ����� �Ÿ�)
  // ��ü PM-curve ��ǥ�� ���̺�(1:��������  final:��������)
  CMap<int,int,_UMD_PMCV_UNIT_POLY,_UMD_PMCV_UNIT_POLY> arPmcvUnit;

	void Init() 
	{
		dRotate = 0.0;
		dXb = 0.0;
		dPb = 0.0;
		dMby = 0.0;
		dMbz = 0.0;
		desib = 0.0;
		dX0 = 0.0;
		dM0y = 0.0;
		dM0z = 0.0;
		dDmax = 0.0;
		dDeff = 0.0;
		arPmcvUnit.RemoveAll();
	}
	_UMD_PMCV_DATA_POLY() { Init(); }
	_UMD_PMCV_DATA_POLY(const _UMD_PMCV_DATA_POLY &src) { *this = src; }
	_UMD_PMCV_DATA_POLY& operator= (const _UMD_PMCV_DATA_POLY &src) 
	{
		dRotate = src.dRotate;
		dXb = src.dXb;
		dPb = src.dPb;
		dMby = src.dMby;
		dMbz = src.dMbz;
		desib = src.desib;
		dX0 = src.dX0;
		dM0y = src.dM0y;
		dM0z = src.dM0z;
		dDmax = src.dDmax;
		dDeff = src.dDeff;
		arPmcvUnit.RemoveAll();
		
		int rKey = 0;
		_UMD_PMCV_UNIT_POLY rData;

		arPmcvUnit.RemoveAll();
		POSITION pos = src.arPmcvUnit.GetStartPosition();
		while(pos)
		{			
			src.arPmcvUnit.GetNextAssoc(pos, rKey, rData);
			arPmcvUnit.SetAt(rKey, rData);
		}
		return *this; 
	}
};

struct _UMD_3DPM_DATA
{	
	int nType;                 // Phi������� (1:Phi���, 2:Phi�̰��)
	int iDivisionNumber_90deg; // 3D-PM����� 90���� ������ ���� (defult=6 : 15���������� ����)  
  int iSymmetryType;         // Section�� ��Ī ����(1:���Ī  2:����,�¿� ��Ī  3:z�࿡ ��Ī(�¿��Ī)  4:y�࿡ ��Ī(���ϴ�Ī)  5:������ ��Ī) �ذ��� ����(�ź���\RC_Column_MPhiCurve�ڷ�.ppt 2Page)
	CArray<_UMD_PMCV_DATA_POLY,_UMD_PMCV_DATA_POLY&>	PmcvData;	
	// �ʱ�ȭ
  void Init()
	{
		nType = 0;
		iDivisionNumber_90deg = 0;
		iSymmetryType = 0;
		PmcvData.RemoveAll();
	}  
  _UMD_3DPM_DATA() {}  
  _UMD_3DPM_DATA& operator= (const _UMD_3DPM_DATA &src) 
	{
		nType = src.nType;
		iDivisionNumber_90deg = src.iDivisionNumber_90deg;
		iSymmetryType = src.iSymmetryType;
		PmcvData.RemoveAll();
		PmcvData.Copy(src.PmcvData);
		return *this; 
	}
};

struct _UMD_RC_FORCE
{
	BOOL   bUse;
  BOOL   bEqSpecial;
	double dP, dM, dV;
	double dMy, dMz, dVy, dVz;	
	
	void Init() 
	{
		bUse = FALSE;
		bEqSpecial = FALSE;
		dP = dM = dV = 0.0;
		dMy= dMz= dVy= dVz = 0.0;
	}
	_UMD_RC_FORCE() { Init(); }
	_UMD_RC_FORCE& operator= (const _UMD_RC_FORCE &src) 
	{ 
		bUse = src.bUse;
		bEqSpecial = src.bEqSpecial;
		dP = src.dP;
		dM = src.dM;
		dV = src.dV;
		dMy= src.dMy;
		dMz= src.dMz;
		dVy= src.dVy;
		dVz= src.dVz;
		return *this; 
	}	
};


struct _UMD_RC_RBAR_UNIT
{
	double dyz[2];	// 0:y, 1:z.
	// Only Beam
	int iGrupNo;
	// Only Design.
	int    iDgnType; // Design�� ö�ٰ�����(0:Line���(t�� �����Ű�鼭 ö�پ� ���ϴ� ö��) 1:�����ܸ���(t�� ����� ������� ����)) // 1���׸��� ����(090811)Building������ �����
	double dLen;     // Design�� Line����� ö�ٿ��� ����Ǵ� ö���� ��ǥ��(��ü ö���� ������ ���� �ϸ� �����ؼ� �ش� ö���� �д��ϴ� ���� ����)
	double dDgnArea; // Design�� �����ܸ��� ö���� �ܸ� // iDgnType==1�϶����� ����
	// Only Checking.
	double dArea;
	double dDia;
	
	void Init() 
	{
		dyz[0] = dyz[1] = 0.0;
		iGrupNo = 0;
		iDgnType = 0;
		dLen = 0.0;
		dDgnArea = 0.0;
		dArea = 0.0;
		dDia  = 0.0;
	}
	_UMD_RC_RBAR_UNIT()	{}
	_UMD_RC_RBAR_UNIT& operator= (const _UMD_RC_RBAR_UNIT& src) 
	{
		dyz[0]   = src.dyz[0];
		dyz[1]   = src.dyz[1];
		iGrupNo  = src.iGrupNo;
		iDgnType = src.iDgnType;
		dLen     = src.dLen;
		dDgnArea = src.dDgnArea;
		dArea    = src.dArea;
		dDia     = src.dDia;
		return *this; 
	}	
};

struct _UMD_RC_RBAR
{
	CMap<int,int&,_UMD_RC_RBAR_UNIT,_UMD_RC_RBAR_UNIT&>	arRbarUnit;
	void Init()	
	{
		arRbarUnit.RemoveAll();
	}
	_UMD_RC_RBAR()	{}
	_UMD_RC_RBAR& operator = (const _UMD_RC_RBAR &src) 
	{
		arRbarUnit.RemoveAll();
		
		_UMD_RC_RBAR_UNIT Unit;
		int iIndex=0;
		POSITION Pos = src.arRbarUnit.GetStartPosition();
		while(Pos)
		{
			Unit.Init();
			src.arRbarUnit.GetNextAssoc(Pos, iIndex, Unit);
			arRbarUnit.SetAt(iIndex, Unit);
		}
		return *this; 
	}
};

// ���� ���Ʈ�� ������ ������ ����ü
struct _UMD_PMCV_UNIT_BALANCE
{
	double dPb;     // ���� ���Ʈ���� ������
	double dMby;    // ���� ���Ʈ���� My
	double dMbz;    // ���� ���Ʈ���� Mz
	double dDmax;   // ���� ���Ʈ���� �ܸ����
  double dDeffb;  // ���� ���Ʈ���� ��ȿ����
	double dCb;     // ���౸���� ����
	double dAngleb; // �߸���ȸ����(Radian)
	double dXb;     // �߸����� ��ġ
	double desib;   // ���� ���Ʈ���� �ִ� ����ö���� ������ (����� '-')
  // �ʱ�ȭ
  void Init()
  {
		dPb=dMby=dMbz=0.0;
    dDmax=dDeffb=dCb=dAngleb=dXb=desib=0.0;
  }
	_UMD_PMCV_UNIT_BALANCE() {}
	_UMD_PMCV_UNIT_BALANCE& operator= (const _UMD_PMCV_UNIT_BALANCE &src) { return *this; }
};


// �̼Ҹ����� ������ ������ ����ü(�����̽�ģ �̼ҿ��)
struct _UMD_RC_CON_PART_UNIT
{
  // �̼ҿ���� �߽���ǥ([0]:y  [1]:z)
	double dyz[2];	// Center Position.
  // �̼ҿ���� ����
  double dArea;
  // �ʱ�ȭ
	void Init()
	{
		dyz[0] = 0.0;
		dyz[1] = 0.0;
		dArea = 0.0;
	}
	_UMD_RC_CON_PART_UNIT()	{}
	_UMD_RC_CON_PART_UNIT(const _UMD_RC_CON_PART_UNIT &src)	{*this = src;}
	_UMD_RC_CON_PART_UNIT& operator= (const _UMD_RC_CON_PART_UNIT &src) 
	{
		dyz[0] = src.dyz[0];
		dyz[1] = src.dyz[1];
		dArea  = src.dArea;
		return *this; 
	}	
};

struct _UMD_RC_CON_PART_LIST   
{
	CArray<_UMD_RC_CON_PART_UNIT,_UMD_RC_CON_PART_UNIT&> List;
	double dYmax;
	double dZmax;
	double dYmin;
	double dZmin;
	void Init() {}
	_UMD_RC_CON_PART_LIST()	{}
	_UMD_RC_CON_PART_LIST(const _UMD_RC_CON_PART_LIST& rData)	{*this = rData;}
	_UMD_RC_CON_PART_LIST& operator= (const _UMD_RC_CON_PART_LIST &src)
	{
		List.RemoveAll();
		List.Copy(src.List);
		dYmax = src.dYmax;
		dZmax = src.dZmax;
		dYmin = src.dYmin;
		dZmin = src.dZmin;
		return *this; 
	}
};

struct _UMD_RC_RBAR_LIST
{
	CArray<_UMD_RC_RBAR_UNIT,_UMD_RC_RBAR_UNIT&> List;
	void Init()
	{
		List.RemoveAll();
	}
	_UMD_RC_RBAR_LIST()	{}
	_UMD_RC_RBAR_LIST(const _UMD_RC_RBAR_LIST &src)	{*this = src;}
	_UMD_RC_RBAR_LIST& operator= (const _UMD_RC_RBAR_LIST &src)
	{
		List.RemoveAll();
		List.Copy(src.List);
		return *this; 
	}	
};

struct _UMD_RC_COL_ROTATE_SECT
{
	int  nDataType;  // �����̽�ģ ��ũ��Ʈ�� �ܸ����������� �ִ��� ����(0:�ܸ����� ��, 1:�����̽����� ������, 2:�����̽����� ����)
	int  iSymmetryType; // Section�� ��Ī ����(1:���Ī  2:����,�¿� ��Ī  3:z�࿡ ��Ī(�¿��Ī)  4:y�࿡ ��Ī(���ϴ�Ī)  5:������ ��Ī) �ذ��� ����(�ź���\RC_Column_MPhiCurve�ڷ�.ppt 2Page)	
	int  iDivisionNumber_90deg;// M-�ջ���� 90���� ������ ���� (defult=6 : 15���������� ����)  
	CArray<_UMD_RC_CON_PART_LIST,_UMD_RC_CON_PART_LIST&>  arConUnitList; // �߸����� ���⿡ ���� ������ �����̽�ģ ��ũ��Ʈ�̼ҿ�������� �迭(0��~90-360�������� ����:��Ī���¿� ���� �ٸ�)  
	CArray<_UMD_RC_RBAR_LIST,_UMD_RC_RBAR_LIST&>	arRbarUnitList;// �߸����� ���⿡ ���� ������ ö�ٿ�������� �迭(0��~90-360�������� ����:��Ī���¿� ���� �ٸ�)  
	void Init()
	{
		nDataType = 0;
		iSymmetryType = 0;
		iDivisionNumber_90deg = 0;
		arConUnitList.RemoveAll();
		arRbarUnitList.RemoveAll();
	}
	_UMD_RC_COL_ROTATE_SECT()	{}
	_UMD_RC_COL_ROTATE_SECT(const _UMD_RC_COL_ROTATE_SECT &src)	{*this = src;}
	_UMD_RC_COL_ROTATE_SECT& operator = (const _UMD_RC_COL_ROTATE_SECT &src) 
	{ 
		nDataType = src.nDataType;
		iSymmetryType = src.iSymmetryType;
		iDivisionNumber_90deg = src.iDivisionNumber_90deg;
		arConUnitList.RemoveAll();
		arConUnitList.Copy(src.arConUnitList);
		arRbarUnitList.RemoveAll();
		arRbarUnitList.Copy(src.arRbarUnitList);
		return *this; 
	}	
};

struct _UMD_RC_GSEC_VERTEX
{
	double dX;  
  double dY;  
	void Init()	{	dX = dY = 0.0; }
	_UMD_RC_GSEC_VERTEX() { }  	
	_UMD_RC_GSEC_VERTEX(double x, double y) { dX = x; dY = y; }
  _UMD_RC_GSEC_VERTEX& Set(double x, double y) { dX = x; dY = y; return *this; }
  _UMD_RC_GSEC_VERTEX(const _UMD_RC_GSEC_VERTEX& src) { *this = src; }
  _UMD_RC_GSEC_VERTEX& operator= (const _UMD_RC_GSEC_VERTEX& src)
  {
    dX    = src.dX;    
		dY    = src.dY;    
    return *this;
  }
};

struct _UMD_RC_GSEC_POLYGON
{
	CArray<_UMD_RC_GSEC_VERTEX, _UMD_RC_GSEC_VERTEX&> aVertex;
	void Init()	{	aVertex.RemoveAll(); }
	_UMD_RC_GSEC_POLYGON() { }  	
	_UMD_RC_GSEC_POLYGON(const _UMD_RC_GSEC_POLYGON &src) { *this = src; }
  _UMD_RC_GSEC_POLYGON& operator= (const _UMD_RC_GSEC_POLYGON &src)
  {
		aVertex.RemoveAll();
    aVertex.Copy(src.aVertex);
    return *this;
  }
	void GetBoundary(double& dXmax, double& dXmin, double& dYmax, double& dYmin)
	{
		dXmax=dXmin=dYmax=dYmin=0.0;
		int nSize = aVertex.GetSize();
		if(nSize>0)
		{
			dXmax=dXmin = aVertex[0].dX;
			dYmax=dYmin = aVertex[0].dY;
		}
		for(int i=1 ; i<nSize ; i++)
		{
			dXmax = max(dXmax, aVertex[i].dX);
			dXmin = min(dXmin, aVertex[i].dX);
			dYmax = max(dYmax, aVertex[i].dY);
			dYmin = min(dYmin, aVertex[i].dY);
		}
	}
};


struct _UMD_RC_CON_POLY_DATA
{
	// �ܸ��� �ܰ�����
	_UMD_RC_GSEC_POLYGON OutPoly;
  // �ܸ��� ���ζ����� �迭
	CArray<_UMD_RC_GSEC_POLYGON,_UMD_RC_GSEC_POLYGON&> arInnPoly;
	// �ʱ�ȭ
  void Init()
	{    
		OutPoly.Init();
		arInnPoly.RemoveAll();
	}
	_UMD_RC_CON_POLY_DATA() { }  	
	_UMD_RC_CON_POLY_DATA(const _UMD_RC_CON_POLY_DATA &src) { *this = src; }
  _UMD_RC_CON_POLY_DATA& operator= (const _UMD_RC_CON_POLY_DATA &src)
	{
		OutPoly = src.OutPoly;
		arInnPoly.RemoveAll();
		arInnPoly.Copy(src.arInnPoly);
		return *this; 
	}
};

struct _UMD_RC_CONC_UNIT
{
	int iGrupNo;
	double dyz[2];				// Center(y,z).
	double dArea;
	double dyzPnt[4][2];	// Corner(y,z).
	void Init()
	{
		iGrupNo = 0;
		dyz[0] = dyz[1] = 0.0;
		dArea = 0.0;
		for(int i=0; i<4; ++i)
		{
			dyzPnt[i][0]=0.0;
			dyzPnt[i][1]=0.0;
		}
	}
	
	_UMD_RC_CONC_UNIT()	{}
	_UMD_RC_CONC_UNIT(_UMD_RC_CONC_UNIT& ConcUnit)	{*this = ConcUnit;}
	_UMD_RC_CONC_UNIT& operator= (const _UMD_RC_CONC_UNIT &src)
	{
		iGrupNo = src.iGrupNo;
		dyz[0] = src.dyz[0];
		dyz[1] = src.dyz[1];
		dArea = src.dArea;
		for(int i=0; i<4; ++i)
		{
			dyzPnt[i][0] = src.dyzPnt[i][0];
			dyzPnt[i][1] = src.dyzPnt[i][1];
		}
		return *this; 
	}
};

struct _UMD_RC_CONC
{
	CMap<int,int&,_UMD_RC_CONC_UNIT,_UMD_RC_CONC_UNIT&>	arConcUnit;
	void Init()	
	{
		arConcUnit.RemoveAll();
	}
	_UMD_RC_CONC()	{}
	_UMD_RC_CONC& operator= (const _UMD_RC_CONC &src) 
	{
		arConcUnit.RemoveAll();
		int iIndex = 0;
		POSITION pos = src.arConcUnit.GetStartPosition();
		while(pos)
		{
			_UMD_RC_CONC_UNIT rData;
			src.arConcUnit.GetNextAssoc(pos, iIndex, rData);
			arConcUnit.SetAt(iIndex, rData);
		}
		return *this; 
	}	
};

struct _UMD_RC_GEN_RBAR_UNIT
{
	double dyz[2];	// 0:y, 1:z.
	// Only Checking.
	double dArea;
	double dDia;
	void Init()
	{
		dyz[0]	 = 0.0;
		dyz[1]	 = 0.0;
		dArea    = 0.0;
		dDia     = 0.0;
	}
	_UMD_RC_GEN_RBAR_UNIT()	{Init();}
	_UMD_RC_GEN_RBAR_UNIT(double dy, double dz, double dSrcDia, double dSrcArea)
	{
		dyz[0] = dy;
		dyz[1] = dz;
		dDia   = dSrcDia;
		dArea  = dSrcArea;
	}
	_UMD_RC_GEN_RBAR_UNIT(const _UMD_RC_GEN_RBAR_UNIT &src)	{*this = src;}
	_UMD_RC_GEN_RBAR_UNIT& operator = (const _UMD_RC_GEN_RBAR_UNIT &src) 
	{ 
		dyz[0] = src.dyz[0];
		dyz[1] = src.dyz[1];
		dArea  = src.dArea;
		dDia   = src.dDia;
		return *this; 
	}
};

struct _UMD_RC_COL_MAINRBAR
{
	//DB�ܸ�
	double	dDc[5];	           //	�ܸ� �ܰ����� �Ÿ�
	int	    iBarNum[3][5];	   //	ö�� ����       (Check : ���(0 �̻�) )
	double  dBarDiaNa1[3][5];  //	Size1 ö�� ����
	double  dBarDiaNa2[3][5];  //	Size2 ö�� ����
	double  dBarAreaNa1[3][5]; //	Size1 ö��1�� �ܸ���
	double  dBarAreaNa2[3][5]; //	Size2 ö��1�� �ܸ���
  CString strBarNa1[3][5];   // Size1 ö�� �̸�
  CString strBarNa2[3][5];   // Size2 ö�� �̸�
	//���Ǵܸ�
	CArray<_UMD_RC_GEN_RBAR_UNIT, _UMD_RC_GEN_RBAR_UNIT&> arGenRbar;

	void Init() 
	{
		for(int i=0; i<5; ++i)
		{
			dDc[i] = 0.0;
			for(int j=0; j<3; ++j)
			{
				iBarNum[j][i]     = 0;
				dBarDiaNa1[j][i]  = 0.0;
				dBarDiaNa2[j][i]  = 0.0;
				dBarAreaNa1[j][i] = 0.0;
				dBarAreaNa2[j][i] = 0.0;
				strBarNa1[j][i]   = _T("");
				strBarNa2[j][i]   = _T("");
			}
		}
	}
	_UMD_RC_COL_MAINRBAR() { Init(); }
	_UMD_RC_COL_MAINRBAR(const _UMD_RC_COL_MAINRBAR &src) { *this = src; }
	_UMD_RC_COL_MAINRBAR& operator= (const _UMD_RC_COL_MAINRBAR &src) 
	{
		for(int i=0; i<5; ++i)
		{
			dDc[i] = src.dDc[i];
			for(int j=0; j<3; ++j)
			{
				iBarNum[j][i]     = src.iBarNum[j][i];
				dBarDiaNa1[j][i]  = src.dBarDiaNa1[j][i];
				dBarDiaNa2[j][i]  = src.dBarDiaNa2[j][i];
				dBarAreaNa1[j][i] = src.dBarAreaNa1[j][i];
				dBarAreaNa2[j][i] = src.dBarAreaNa2[j][i];
				strBarNa1[j][i]   = src.strBarNa1[j][i];
				strBarNa2[j][i]   = src.strBarNa2[j][i];
			}
		}
		return *this; 
	}
};

#define _UMD_RC_GSEC_POLYGON_LIST CArray<_UMD_RC_GSEC_POLYGON, _UMD_RC_GSEC_POLYGON&>


class DGL_3dp
{
protected:
  
public:
  double m_p[3];
  
  DGL_3dp()
  {
    m_p[0] = m_p[1] = m_p[2] = 0;
  }
  
  ~DGL_3dp()
  {
		
  }
	
  DGL_3dp(DGL_3dp& src)
  {
    m_p[0] = src.m_p[0];
    m_p[1] = src.m_p[1];
    m_p[2] = src.m_p[2];
    
  }
  
  double x(){ return  m_p[0]; }
  double y(){ return  m_p[1]; }
  double z(){ return  m_p[2]; }
  
  
  operator void  * ()
  {
    return m_p;
  }
	
  operator double* ()
  {
    return m_p;
  }
	
  DGL_3dp& operator = (const DGL_3dp& src)
  {
    m_p[0] = src.m_p[0];
    m_p[1] = src.m_p[1];
    m_p[2] = src.m_p[2];
    return *this;
  }
	
  void Set(double x, double y, double z)
  {
    m_p[0] = x; 
    m_p[1] = y;
    m_p[2] = z;
  }
	
  void Get(double& x, double& y, double & z)
  {
    x = m_p[0]; 
    y = m_p[1];
    z = m_p[2];
  }
	
  void Get(double rP[3])
  {
    rP[0] = m_p[0]; 
    rP[1] = m_p[1];
    rP[2] = m_p[2];
  }
	
  void Set3d(double p[3])
  {
    m_p[0] = p[0]; 
    m_p[1] = p[1]; 
    m_p[2] = p[2]; 
  }
};

class DGL_Line3d
{
public:
  DGL_Line3d()
  {
    ;
  }
  ~DGL_Line3d()
  {
		
  }
	
  DGL_3dp m_P1;
  DGL_3dp m_P2;
	
  DGL_Line3d& operator = (const DGL_Line3d& src)
  {
    m_P1 = src.m_P1;
    m_P2 = src.m_P2;
    return *this;
  }
};


class DgnPoint
{
public:
  int id;
  double density;
  double uv[2];
  double xyz[3];
public:
  DgnPoint() : id(0), density(0.0) {}
  DgnPoint(const DgnPoint& point)
  {
    id      = point.id;
    density = point.density;
    uv[0]   = point.uv[0];
    uv[1]   = point.uv[1];
    xyz[0]  = point.xyz[0];
    xyz[1]  = point.xyz[1];
    xyz[2]  = point.xyz[2];
  }
  DgnPoint& operator=(const DgnPoint& point)
  {
    if (&point == this) return *this;
    id      = point.id;
    density = point.density;
    uv[0]   = point.uv[0];
    uv[1]   = point.uv[1];
    xyz[0]  = point.xyz[0];
    xyz[1]  = point.xyz[1];
    xyz[2]  = point.xyz[2];
    return *this;
  }
  ~DgnPoint() {}
};


#endif // !defined(__DGNCALC_KR_RCSC_COLUMNSTRUCT_H__)