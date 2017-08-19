// PolyMaker.h: interface for the CPolyMaker class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_POLYMAKER_H__EA020F37_5794_40E8_9256_2ED1F21EE6C8__INCLUDED_)
#define AFX_POLYMAKER_H__EA020F37_5794_40E8_9256_2ED1F21EE6C8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <afxtempl.h>
#include "StructDef.h"

struct MVector3d
{
	double vec[3];
	void Init()
	{
		vec[0] = 0.0;
		vec[1] = 0.0;
		vec[2] = 0.0;
	}
	MVector3d() { Init(); }
	MVector3d(const double &val)
	{
		vec[0] = val;
		vec[1] = val;		
		vec[2] = val;		
	}
	MVector3d& operator= (const MVector3d &src)
	{
		vec[0] = src.vec[0];
		vec[1] = src.vec[1];
		vec[2] = src.vec[2];
		return *this;
	}
};



class CPolyMaker  
{
public:
	CPolyMaker();
	virtual ~CPolyMaker();

protected:

  CArray<DGL_3dp,DGL_3dp&>  m_arOPolyV  ;
  CArray<CArray<DGL_3dp,DGL_3dp&>*,CArray<DGL_3dp,DGL_3dp&>*> m_arHolePolys;
  CArray<DGL_3dp,DGL_3dp&>  m_CuttingLine;
  
  CArray<CArray<DGL_3dp,DGL_3dp&>*,CArray<DGL_3dp,DGL_3dp&>*> m_ResultPoly;
  CArray<CArray<DGL_3dp,DGL_3dp&>*,CArray<DGL_3dp,DGL_3dp&>*> m_ResultHolePoly;
  CArray<DGL_Line3d,DGL_Line3d& >                             m_arCuttingEdge ;
	
  CArray<CArray<DGL_3dp,DGL_3dp&>*,CArray<DGL_3dp,DGL_3dp&>*> m_arCutPolys1_Slice;
  CArray<CArray<DGL_3dp,DGL_3dp&>*,CArray<DGL_3dp,DGL_3dp&>*> m_arCutPolys2_Slice;
	
public:
  void SetOutterPoly   (CArray<DGL_3dp,DGL_3dp&>&OVerts); //�ܰ� ������ ���� ���� 
  void AddInnerPoly    (CArray<DGL_3dp,DGL_3dp&>&IVerts); //���� ������ ���� ���� 
  void SetCuttingLine  (DGL_3dp& P1, DGL_3dp& P2);        //���ܼ� ���� 
  
  BOOL MakeResult      ();                                //��ܸ� ���� ���� (������:������ ������ ���ܼ����� �ɸ��� �ʰ� ������ũ���� 2�� �̻� ������ ���� ��� �߸��� ����� �����Ҽ� ����)
  
  //------------------------------------------------------------------------------------------
  int  GetResultPolyCount ();  // ������ ������ ���� ��ȯ                             
  BOOL GetResultPoly      (int nIndex,CArray<DGL_3dp, DGL_3dp&> &RPoly);
  //------------------------------------------------------------------------------------------
  int  GetResultHolePolyCount( );  // ������ �󱸰����� Hole�� ���� ��ȯ 
  BOOL GetResultHolePoly  (int nIndex,CArray<DGL_3dp, DGL_3dp&> &RPoly);
  //----------------------------------------------------------------------------------------------
  int  GetCuttingEdgeCount();  // ������ ��ȿ�� ���ܼ� ���� ��ȯ 
  //------------------------------------------------------------------------------------------
  BOOL GetCuttingEdge     (int nIndex,DGL_3dp& P1, DGL_3dp& P2);
  //------------------------------------------------------------------------------------------
public:
  BOOL GetBoundingBox(DGL_3dp& MinP, DGL_3dp& MaxP);
  void RemoveAllInnerPolys();
  void RemoveAllOutterPoly();
  void ResetAll();
  void RemoveAllCutPolys();
  
public: // Sliced Polygon���� ���
  void RemoveAllResult();
	
  void SetSlicingLines   (CArray<DGL_Line3d,DGL_Line3d&>& arSlicingLines);
  BOOL MakeSlicedPolygons();      
  BOOL MakeSlicedPolygons2();      
  
  void AddCuttingPoly     (CArray<DGL_3dp,DGL_3dp&>&CVerts);
  BOOL MakeCuttingResult_Intersect();
  BOOL MakeCuttingResult_Extract  ();
  BOOL MakeCuttingResult_Merge    ();
  
  void GetOutterPoly        (CArray<DGL_3dp,DGL_3dp&>&OVerts); 
  int  GetInnerPolyCount    (); 
  void GetInnerPoly         (int nIndex, CArray<DGL_3dp,DGL_3dp&>&IVerts);
protected:
  BOOL _MakeCuttingResult(int nMode); // nMode : (1) Intersect (2)Extract (3) Merge
	
  BOOL PolygonShrink(CArray<DGL_3dp,DGL_3dp&>& arPolyVertsSrc,CArray<DGL_3dp,DGL_3dp&>& arPolyVertsShrinked,
                             double SFactorOrDist,BOOL bIsFactor);

};

#endif // !defined(AFX_POLYMAKER_H__EA020F37_5794_40E8_9256_2ED1F21EE6C8__INCLUDED_)
