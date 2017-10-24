#include "stdafx.h"

#include "PolyMaker.h"
#define __INC_GRENDERENGINE_H__
  #include "GRenderAll.h"
#include "GClipPolybyVatti.h"

CPolyMaker::CPolyMaker()
{

}

CPolyMaker::~CPolyMaker()
{
	ResetAll();  
}
//寇胞 弃府帮 屈惑 汲沥 
void CPolyMaker::SetOutterPoly(CArray<DGL_3dp,DGL_3dp&>&OVerts) 
{
  m_arOPolyV.RemoveAll();
  m_arOPolyV.SetSize(OVerts.GetSize());
  m_arOPolyV.Copy(OVerts);
}
//郴何 弃府帮 屈惑 汲沥
void CPolyMaker::AddInnerPoly(CArray<DGL_3dp,DGL_3dp&>&IVerts)  
{
  CArray<DGL_3dp,DGL_3dp&> * pVerts;
  pVerts = new CArray<DGL_3dp,DGL_3dp&>;
  pVerts->SetSize(IVerts.GetSize());
  pVerts->Copy(IVerts);
  m_arHolePolys.Add(pVerts);
}
//例窜急 沥狼 
void CPolyMaker::SetCuttingLine(DGL_3dp& P1, DGL_3dp& P2)        
{
  m_CuttingLine.RemoveAll();
  this->m_CuttingLine.Add(P1);
  this->m_CuttingLine.Add(P2);
}
//企窜搁 沥焊 积己 
BOOL CPolyMaker::MakeResult()                                
{
	RemoveAllResult();

  GClipPolyBy_Vatti_algorithm Alg_Vatti;

  //GClosedPolygonMaker* pGCM = GClosedPolygonMaker::GetInstance();
  int nOPolyV = m_arOPolyV.GetSize();
  int nIPolys = m_arHolePolys.GetSize();
  int nCLineV = m_CuttingLine.GetSize();

  if(nOPolyV <3 ) return FALSE;
  if(nCLineV <2 ) return FALSE;

  Alg_Vatti.AddSubjectPoly(m_arOPolyV,FALSE);
  
  for(int i = 0; i < nIPolys; i++)
  {
    if(3 > m_arHolePolys[i]->GetSize()) continue; // 公瓤茄 Polygon 力寇 
    Alg_Vatti.AddSubjectPoly(*m_arHolePolys[i],TRUE);
  }

  GVector CLineVector, YVector,ZVector;
  GPoint3D P1, P2,P3,P4;
  
  P1.Set(m_CuttingLine[0].x(),m_CuttingLine[0].y(),0);
  P2.Set(m_CuttingLine[1].x(),m_CuttingLine[1].y(),0);

  DGL_3dp MinP,MaxP;
  GLine3D TempLine;
  GetBoundingBox(MinP,MaxP);

  TempLine.P1.Set(MinP.x(),MinP.y(),0);
  TempLine.P2.Set(MaxP.x(),MaxP.y(),0);

  double LOffset, YOffset;
  LOffset = TempLine.Length();
  YOffset = LOffset*2.;
  
  ZVector.Set(0,0,1);
  CLineVector.Set(P1,P2); CLineVector.MakeUnit();
  YVector = ZVector * CLineVector;

  CLineVector.ScalingMySelf(LOffset);

  P1.x = P1.x - CLineVector.x;
  P1.y = P1.y - CLineVector.y;
  P2.x = P2.x + CLineVector.x;
  P2.y = P2.y + CLineVector.y;
    
  YVector.ScalingMySelf(YOffset);

  P3.x = P2.x + YVector.x;
  P3.y = P2.y + YVector.y;
  P4.x = P1.x + YVector.x;
  P4.y = P1.y + YVector.y;
  
  CArray<GPoint3D, GPoint3D&> arClipPoly;

  arClipPoly.Add(P1); arClipPoly.Add(P2);  
  arClipPoly.Add(P3); arClipPoly.Add(P4);


  Alg_Vatti.AddClipPoly(arClipPoly);

  //--------------------- 企窜搁 积己 ----------------------------
  Alg_Vatti.MakeResult_intersect();
  
  //--------------------- 积己等 企窜搁 历厘 ---------------------
  CArray<DGL_3dp,DGL_3dp&>  *pVTemp;
  CArray<int,int> arIPolys;
  
  Alg_Vatti.GetResultPolyIndex(FALSE,arIPolys);
  int nClosedPoly = arIPolys.GetSize();

  for(int i = 0; i < nClosedPoly; i++)
  {
    pVTemp = new CArray<DGL_3dp,DGL_3dp&>;
    Alg_Vatti.GetResultPoly(arIPolys[i],*pVTemp);
    m_ResultPoly.Add(pVTemp); 
  }

  //--------------------- 积己等 Hole 历厘 -----------------------
  Alg_Vatti.GetResultPolyIndex(TRUE,arIPolys);
  int nHolePoly = arIPolys.GetSize();

  for(int i = 0; i < nHolePoly; i++)
  {
    pVTemp = new CArray<DGL_3dp,DGL_3dp&>;
    Alg_Vatti.GetResultPoly(arIPolys[i],*pVTemp);
    m_ResultHolePoly.Add(pVTemp); 
  }

  //--------------------- Cutting Line 历厘 ----------------------
  GLine3D CuttingLine;
  CuttingLine.Set(P1,P2);
  CArray<CuttingEdge,CuttingEdge&>arEdge;
  if(!Alg_Vatti.GetCuttingEdge(arEdge,CuttingLine)) return FALSE;

  int narEdge;
  narEdge = arEdge.GetSize();
  DGL_Line3d TempL;
  for(int i = 0; i < narEdge; i++)
  {
    TempL.m_P1.Set(arEdge[i].m_Pos1.x,arEdge[i].m_Pos1.y,0);
    TempL.m_P2.Set(arEdge[i].m_Pos2.x,arEdge[i].m_Pos2.y,0);
    m_arCuttingEdge.Add(TempL);
  }
  
  return TRUE;
}

//------------------------------------------------------------------------------------------
// 积己等 弃府帮 俺荐 馆券  
int CPolyMaker::GetResultPolyCount()                             
{
  return m_ResultPoly.GetSize();
}

BOOL CPolyMaker::GetResultPoly(int nIndex,CArray<DGL_3dp, DGL_3dp&> &RPoly)
{
  if(nIndex > m_ResultPoly.GetSize()-1) return FALSE;
  RPoly.RemoveAll();
  RPoly.Copy(*(m_ResultPoly[nIndex]));
  return TRUE;
}

//------------------------------------------------------------------------------------------
// 积己等 企备埃俊辑 Hole狼 俺荐 馆券 
int  CPolyMaker::GetResultHolePolyCount( )  
{
  return this->m_ResultHolePoly.GetSize();
}

BOOL CPolyMaker::GetResultHolePoly  (int nIndex,CArray<DGL_3dp, DGL_3dp&> &RPoly)
{
  if(nIndex > m_ResultHolePoly.GetSize()-1) return FALSE;
  RPoly.RemoveAll();
  RPoly.Copy(*(m_ResultHolePoly[nIndex]));
  return TRUE; 
}

//------------------------------------------------------------------------------------------
// 积己等 蜡瓤茄 例窜急 俺荐 馆券 
int  CPolyMaker::GetCuttingEdgeCount()  
{
  return m_arCuttingEdge.GetSize();
}

//------------------------------------------------------------------------------------------
BOOL CPolyMaker::GetCuttingEdge(int nIndex,DGL_3dp& P1, DGL_3dp& P2)
{
  if(nIndex > m_arCuttingEdge.GetSize()-1) return FALSE;

  P1 = m_arCuttingEdge[nIndex].m_P1;
  P2 = m_arCuttingEdge[nIndex].m_P2;  

  return TRUE;
}

//------------------------------------------------------------------------------------------

BOOL CPolyMaker::GetBoundingBox(DGL_3dp& MinP, DGL_3dp& MaxP)
{
  int nOPolyV = m_arOPolyV.GetSize();
  int nIPolys = m_arHolePolys.GetSize();

  if(nOPolyV == 0) return FALSE;

  double MinX, MinY, MaxX ,MaxY;

  MinX = MinY = DBL_MAX;
  MaxX = MaxY = -DBL_MAX;
 

  for(int i = 0 ; i < nOPolyV; i++)
  {
    MinX = min(m_arOPolyV[i].x(),MinX); 
    MinY = min(m_arOPolyV[i].y(),MinY); 
    MaxX = max(m_arOPolyV[i].x(),MaxX); 
    MaxY = max(m_arOPolyV[i].y(),MaxY); 
  }

  for(int i = 0; i < nIPolys; i++)
  {
    int nIVert = m_arHolePolys[i]->GetSize();
    for( int j = 0; j < nIVert; j++)
    {
      MinX = min((*m_arHolePolys[i])[j].x(),MinX); 
      MinY = min((*m_arHolePolys[i])[j].y(),MinY); 
      MaxX = max((*m_arHolePolys[i])[j].x(),MaxX); 
      MaxY = max((*m_arHolePolys[i])[j].y(),MaxY);  
    }
  }

  MinP.Set(MinX,MinY,0);
  MaxP.Set(MaxX,MaxY,0);
  
  return TRUE;
}


void CPolyMaker::ResetAll()
{
  RemoveAllInnerPolys();
  m_arOPolyV.RemoveAll();
  m_CuttingLine.RemoveAll();
  RemoveAllResult();
  RemoveAllCutPolys();
}

void CPolyMaker::RemoveAllOutterPoly()
{
  m_arOPolyV.RemoveAll();
}

void CPolyMaker::RemoveAllInnerPolys()
{
  int nH = m_arHolePolys.GetSize();
  for( int i = 0; i < nH ; i++)
  {
    delete m_arHolePolys[i];
  }
  m_arHolePolys.RemoveAll();
}

void CPolyMaker::RemoveAllCutPolys()
{
  int nP = m_arCutPolys1_Slice.GetSize();
  for(int i = 0; i < nP; i++)
  {
    delete m_arCutPolys1_Slice[i];
  }
  
  m_arCutPolys1_Slice.RemoveAll();

  nP = m_arCutPolys2_Slice.GetSize();

  for(int i = 0; i < nP ; i++)
  {
    delete m_arCutPolys2_Slice[i];
  }

  m_arCutPolys2_Slice.RemoveAll();
}

void CPolyMaker::RemoveAllResult()
{
  int nP = m_ResultPoly.GetSize();
  for(int i = 0; i < nP ; i++)
  {
    delete m_ResultPoly[i];
  }
  m_ResultPoly.RemoveAll();
  
  nP = m_ResultHolePoly.GetSize();
  for(int i = 0; i < nP ; i++)
  {
    delete m_ResultHolePoly[i];
  }
  m_ResultHolePoly.RemoveAll();

  m_arCuttingEdge.RemoveAll();
}

void CPolyMaker::SetSlicingLines(CArray<DGL_Line3d,DGL_Line3d&>& arSlicingLines) 
{
  RemoveAllCutPolys();

  CArray<DGL_3dp,DGL_3dp&> * pPolys;

  int nSlicingLines = arSlicingLines.GetSize();
  for(int i = 0; i < nSlicingLines-1; i++)
  {
    pPolys = new CArray<DGL_3dp,DGL_3dp&>;
    
    pPolys->Add(arSlicingLines[i].m_P1);
    pPolys->Add(arSlicingLines[i].m_P2);
    pPolys->Add(arSlicingLines[i+1].m_P2);
    pPolys->Add(arSlicingLines[i+1].m_P1);

    if(i%2 == 0) // Even
      this->m_arCutPolys1_Slice.Add(pPolys);
    else         // Odd
      this->m_arCutPolys2_Slice.Add(pPolys);
  }
}

BOOL CPolyMaker::MakeSlicedPolygons()
{
  RemoveAllResult();

  GClipPolyBy_Vatti_algorithm Alg_VattiOdd ;
  GClipPolyBy_Vatti_algorithm Alg_VattiEven;

  int nOPolyV = m_arOPolyV.GetSize();
  int nIPolys = m_arHolePolys.GetSize();
  
  //----------------------------------------------------------------  !!!
  // 
  if(nOPolyV <3 ) return FALSE;
  
  Alg_VattiOdd  .AddSubjectPoly(m_arOPolyV,FALSE);
  Alg_VattiEven .AddSubjectPoly(m_arOPolyV,FALSE);
  
  for(int i = 0; i < nIPolys; i++)
  {
    if(3 > m_arHolePolys[i]->GetSize()) continue; // 公瓤茄 Polygon 力寇 
    Alg_VattiOdd.AddSubjectPoly(*m_arHolePolys[i],TRUE);
    Alg_VattiEven.AddSubjectPoly(*m_arHolePolys[i],TRUE);
  }

  CArray<DGL_3dp,DGL_3dp&>  *pVTemp;
  CArray<int,int> arIPolys;
  int nClosedPoly,nHolePoly,nSlicingPoly;

  //--------------------------------------------------------------
  //
  //  Odd
  //  
  nSlicingPoly = m_arCutPolys1_Slice.GetSize();
  for(int i = 0 ; i < nSlicingPoly; i++)
  {
    Alg_VattiOdd.AddClipPoly(*m_arCutPolys1_Slice[i]);
  }
  
  Alg_VattiOdd.MakeResult_intersect();

  Alg_VattiOdd.GetResultPolyIndex(FALSE,arIPolys);
  nClosedPoly = arIPolys.GetSize();
  
  //------------------ 积己等 Slice历厘 --------------------------
  for(int i = 0; i < nClosedPoly; i++)
  {
    pVTemp = new CArray<DGL_3dp,DGL_3dp&>;
    Alg_VattiOdd.GetResultPoly(arIPolys[i],*pVTemp);
    m_ResultPoly.Add(pVTemp); 
  }

  //--------------------- 积己等 Hole 历厘 -----------------------
  Alg_VattiOdd.GetResultPolyIndex(TRUE,arIPolys);
  nHolePoly = arIPolys.GetSize();

  for(int i = 0; i < nHolePoly; i++)
  {
    pVTemp = new CArray<DGL_3dp,DGL_3dp&>;
    Alg_VattiOdd.GetResultPoly(arIPolys[i],*pVTemp);
    m_ResultHolePoly.Add(pVTemp); 
  }

  //--------------------------------------------------------------
  //
  //  Even
  //  
  nSlicingPoly = m_arCutPolys2_Slice.GetSize();
  for(int i = 0; i < nSlicingPoly; i++)
  {
    Alg_VattiEven.AddClipPoly(*m_arCutPolys2_Slice[i]);
  }

  Alg_VattiEven.MakeResult_intersect();

  Alg_VattiEven.GetResultPolyIndex(FALSE,arIPolys);
  nClosedPoly = arIPolys.GetSize();
  //------------------ 积己等 Slice历厘 --------------------------
  for(int i = 0; i < nClosedPoly; i++)
  {
    pVTemp = new CArray<DGL_3dp,DGL_3dp&>;
    Alg_VattiEven.GetResultPoly(arIPolys[i],*pVTemp);
    m_ResultPoly.Add(pVTemp); 
  }

  //--------------------- 积己等 Hole 历厘 -----------------------
  Alg_VattiEven.GetResultPolyIndex(TRUE,arIPolys);
  nHolePoly = arIPolys.GetSize();
  for(int i = 0; i < nHolePoly; i++)
  {
    pVTemp = new CArray<DGL_3dp,DGL_3dp&>;
    Alg_VattiEven.GetResultPoly(arIPolys[i],*pVTemp);
    m_ResultHolePoly.Add(pVTemp); 
  }
  
  return TRUE;
}
  

BOOL CPolyMaker::MakeSlicedPolygons2()
{
  RemoveAllResult();

  GClipPolyBy_Vatti_algorithm Alg_VattiOdd ;
  GClipPolyBy_Vatti_algorithm Alg_VattiEven;

  int nOPolyV = m_arOPolyV.GetSize();
  int nIPolys = m_arHolePolys.GetSize();
  
  //----------------------------------------------------------------  !!!
  // 
  if(nOPolyV <3 ) return FALSE;
  
  Alg_VattiOdd  .AddSubjectPoly(m_arOPolyV,FALSE);
  Alg_VattiEven .AddSubjectPoly(m_arOPolyV,FALSE);
  
  for(int i = 0; i < nIPolys; i++)
  {
    if(3 > m_arHolePolys[i]->GetSize()) continue; // 公瓤茄 Polygon 力寇 
    Alg_VattiOdd.AddSubjectPoly(*m_arHolePolys[i],TRUE);
    Alg_VattiEven.AddSubjectPoly(*m_arHolePolys[i],TRUE);
  }

  CArray<DGL_3dp,DGL_3dp&>  *pVTemp;
  CArray<int,int> arIPolys;
  int nClosedPoly,nHolePoly,nSlicingPoly;

  //--------------------------------------------------------------
  //
  //  Odd
  //  
  nSlicingPoly = m_arCutPolys1_Slice.GetSize();
  for(int i = 0 ; i < nSlicingPoly; i++)
  {
    Alg_VattiOdd.ResetClipPolys();
    Alg_VattiOdd.ResetResultPolys();

    Alg_VattiOdd.AddClipPoly(*m_arCutPolys1_Slice[i]);
  
    Alg_VattiOdd.MakeResult_intersect();

    Alg_VattiOdd.GetResultPolyIndex(FALSE,arIPolys);
    nClosedPoly = arIPolys.GetSize();
    //------------------ 积己等 Slice历厘 --------------------------
    for(int i = 0; i < nClosedPoly; i++)
    {
      pVTemp = new CArray<DGL_3dp,DGL_3dp&>;
      Alg_VattiOdd.GetResultPoly(arIPolys[i],*pVTemp);
      m_ResultPoly.Add(pVTemp); 
    }

    //--------------------- 积己等 Hole 历厘 -----------------------
    Alg_VattiOdd.GetResultPolyIndex(TRUE,arIPolys);
    nHolePoly = arIPolys.GetSize();
    for(int i = 0; i < nHolePoly; i++)
    {
      pVTemp = new CArray<DGL_3dp,DGL_3dp&>;
      Alg_VattiOdd.GetResultPoly(arIPolys[i],*pVTemp);
      m_ResultHolePoly.Add(pVTemp); 
    }
  }

  //--------------------------------------------------------------
  //
  //  Even
  //  
  nSlicingPoly = m_arCutPolys2_Slice.GetSize();
  for(int i = 0; i < nSlicingPoly; i++)
  {
    Alg_VattiEven.ResetClipPolys();
    Alg_VattiEven.ResetResultPolys();

    Alg_VattiEven.AddClipPoly(*m_arCutPolys2_Slice[i]);

    Alg_VattiEven.MakeResult_intersect();

    Alg_VattiEven.GetResultPolyIndex(FALSE,arIPolys);
    nClosedPoly = arIPolys.GetSize();
    //------------------ 积己等 Slice历厘 --------------------------
    for(int j = 0; j < nClosedPoly; j++)
    {
      pVTemp = new CArray<DGL_3dp,DGL_3dp&>;
      Alg_VattiEven.GetResultPoly(arIPolys[j],*pVTemp);
      m_ResultPoly.Add(pVTemp); 
    }

    //--------------------- 积己等 Hole 历厘 -----------------------
    Alg_VattiEven.GetResultPolyIndex(TRUE,arIPolys);
    nHolePoly = arIPolys.GetSize();
    for(int j = 0; j < nHolePoly; j++)
    {
      pVTemp = new CArray<DGL_3dp,DGL_3dp&>;
      Alg_VattiEven.GetResultPoly(arIPolys[j],*pVTemp);
      m_ResultHolePoly.Add(pVTemp); 
    }
  }
  return TRUE;
}

void CPolyMaker::AddCuttingPoly (CArray<DGL_3dp,DGL_3dp&>&CVerts)
{
  CArray<DGL_3dp,DGL_3dp&>* pPolys = new CArray<DGL_3dp,DGL_3dp&>;
  pPolys->Copy(CVerts);
  m_arCutPolys1_Slice.Add(pPolys);
}

BOOL CPolyMaker::MakeCuttingResult_Intersect()
{
  return _MakeCuttingResult(1);
}

BOOL CPolyMaker::MakeCuttingResult_Extract  ()
{
  return _MakeCuttingResult(2);
}

BOOL CPolyMaker::MakeCuttingResult_Merge    ()
{
  return _MakeCuttingResult(3);
}


void CPolyMaker::GetOutterPoly        (CArray<DGL_3dp,DGL_3dp&>&OVerts)
{
	OVerts.Copy(m_arOPolyV);
}

int CPolyMaker::GetInnerPolyCount    ()
{
	int nIPolys = m_arHolePolys.GetSize(); 
	return nIPolys;
}

void CPolyMaker::GetInnerPoly         (int nIndex, CArray<DGL_3dp,DGL_3dp&>&IVerts)
{
	int nIPolys = m_arHolePolys.GetSize(); 
	if(nIndex >= nIPolys)
	{
		ASSERT(0);
		return;
	}  
	
	IVerts.Copy(*(m_arHolePolys.ElementAt(nIndex)));
}

// nMode : (1) Intersect (2)Extract (3) Merge
BOOL CPolyMaker::_MakeCuttingResult(int nMode)
{
  RemoveAllResult();

  GClipPolyBy_Vatti_algorithm Alg_Vatti;
  
  int nOPolyV = m_arOPolyV.GetSize();
  int nIPolys = m_arHolePolys.GetSize();
  
  //----------------------------------------------------------------  !!!
  // 
  if(nOPolyV <3 ) return FALSE;
  
  Alg_Vatti  .AddSubjectPoly(m_arOPolyV,FALSE);
  
  for(int i = 0; i < nIPolys; i++)
  {
    if(3 > m_arHolePolys[i]->GetSize()) continue; // 公瓤茄 Polygon 力寇 
    Alg_Vatti.AddSubjectPoly(*m_arHolePolys[i],TRUE);
  }

  CArray<DGL_3dp,DGL_3dp&>  *pVTemp;
  CArray<int,int> arIPolys;
  int nClosedPoly,nHolePoly,nSlicingPoly;

  //--------------------------------------------------------------
  //
  //  Cutter汲沥 
  //  
  nSlicingPoly = m_arCutPolys1_Slice.GetSize();
  for(int i = 0 ; i < nSlicingPoly; i++)
  {
    Alg_Vatti.AddClipPoly(*m_arCutPolys1_Slice[i]);
  }
  
  if(1 == nMode)
    Alg_Vatti.MakeResult_intersect();
  else if(2 == nMode)
    Alg_Vatti.MakeResult_extract();
  else if(3 == nMode)
    Alg_Vatti.MakeResult_merge();

  Alg_Vatti.GetResultPolyIndex(FALSE,arIPolys);
  nClosedPoly = arIPolys.GetSize();
  
  //------------------ 积己等 Slice历厘 --------------------------
  for(int i = 0; i < nClosedPoly; i++)
  {
    pVTemp = new CArray<DGL_3dp,DGL_3dp&>;
    Alg_Vatti.GetResultPoly(arIPolys[i],*pVTemp);
    m_ResultPoly.Add(pVTemp); 
  }

  //--------------------- 积己等 Hole 历厘 -----------------------
  Alg_Vatti.GetResultPolyIndex(TRUE,arIPolys);
  nHolePoly = arIPolys.GetSize();

  for(int i = 0; i < nHolePoly; i++)
  {
    pVTemp = new CArray<DGL_3dp,DGL_3dp&>;
    Alg_Vatti.GetResultPoly(arIPolys[i],*pVTemp);
    m_ResultHolePoly.Add(pVTemp); 
  }
  return TRUE;  
}


BOOL CPolyMaker::PolygonShrink(CArray<DGL_3dp,DGL_3dp&>& arPolyVertsSrc,CArray<DGL_3dp,DGL_3dp&>& arPolyVertsShrinked,
                                        double SFactorOrDist,BOOL bIsFactor)
{
  CArray<GPoint3D,GPoint3D&> PolyVertsSrc;
  CArray<GPoint3D,GPoint3D&> PolyVertsShrinked;
  GGeometryEngine GM;
  
  int nPVerts = arPolyVertsSrc.GetSize();
  PolyVertsSrc.SetSize(nPVerts);

  for(int i = 0; i < nPVerts ; i++)
  {
    DGL_3dp& RPos = arPolyVertsSrc.ElementAt(i);
    PolyVertsSrc.ElementAt(i).Set(RPos.x(),RPos.y(),RPos.z());
  }

  BOOL RetVal = GM.PolygonShrink_loveme(PolyVertsSrc,PolyVertsShrinked,SFactorOrDist,bIsFactor);

  if(RetVal)
  {
    nPVerts = PolyVertsShrinked.GetSize();
    arPolyVertsShrinked.SetSize(nPVerts);

    for(int i = 0 ; i < nPVerts; i++)
    {
      GPoint3D& pos3G = PolyVertsShrinked.ElementAt(i);
      arPolyVertsShrinked.ElementAt(i).Set(pos3G.x,pos3G.y,pos3G.z);
    }
  }
  return RetVal;
}
