#include "stdafx.h"
#define __INC_GLINE3D_H__
#include "GRenderAll.h"

#ifdef  _DEBUG
#define new DEBUG_NEW
#undef  THIS_FILE
static  char THIS_FILE[] = __FILE__;
#endif

/**********************************************************************
@@	GLine3D
@@	Viewing 변환에서 사용될 Line Object
@@	Grid 생성시 
*/
GLine3D::GLine3D()
{
  m_YIndex = -1 ; // SDS Mesh Line 관리를 위해 추가됨...
  m_XIndex = -1 ; // SDS Mesh Line 관리를 위해 추가됨...
}

GLine3D::~GLine3D()
{

}


GLine3D::GLine3D(GPoint3D&p1, GPoint3D&p2)
{
 	P1 = p1;
  P2 = p2;
  m_YIndex = -1 ; // SDS Mesh Line 관리를 위해 추가됨...
  m_XIndex = -1 ; // SDS Mesh Line 관리를 위해 추가됨...
}

GLine3D::GLine3D(GLine3D& S)
{
  P1 = S.P1;
  P2 = S.P2;
  
  m_XIndex = S.m_XIndex ; // SDS Mesh Line 관리를 위해 추가됨...
  m_YIndex = S.m_YIndex ; // SDS Mesh Line 관리를 위해 추가됨...
}

GLine3D& GLine3D::operator = (const GLine3D& Source)
{
  P1 = Source.P1;
  P2 = Source.P2;
  /*
  memcpy(&P1.x,&Source.P1.x,sizeof(LDOUBLE));
	memcpy(&P1.y,&Source.P1.y,sizeof(LDOUBLE));
	memcpy(&P1.z,&Source.P1.z,sizeof(LDOUBLE));
	
	memcpy(&P2.x,&Source.P2.x,sizeof(LDOUBLE));
	memcpy(&P2.y,&Source.P2.y,sizeof(LDOUBLE));
	memcpy(&P2.z,&Source.P2.z,sizeof(LDOUBLE));
	*/
  m_XIndex = Source.m_XIndex;
  m_YIndex = Source.m_YIndex;
	return (*this);
}

void GLine3D::Set(const GPoint3D& p1,const GPoint3D& p2)
{
	P1 = p1;
	P2 = p2;
}

void GLine3D::Set(const GLine2D& Line2D)
{
	P1.Set(Line2D.P1);
	P2.Set(Line2D.P2);
}

CString GLine3D::GetString(CString Format)
{
	CString StrTemp;
	StrTemp.Format(Format,
		P1.x,P1.y,P1.z,P2.x,P2.y,P2.z);
	
	return StrTemp;
}

void GLine3D::SetZero()
{
	P1.Set(0.0,0.0,0.0);
	P2.Set(0.0,0.0,0.0);
}

/***********************************************************
@@	현재 Line을 구성하는 점들이 유효한 Line을 구성하는 조건
@@	인지를 검사한다. 즉 Line을 구성하는 두점이 같은 점이면 
@@	Line은 유효하지 않은 Line으로 판정한다. 
**/
LBOOL GLine3D::IsValid()
{
	if(P1 == P2 ) return FALSE;
	return TRUE;
}

LBOOL GLine3D::IsValid2D()
{
  if(P1.IsEqual2D(P2)) return FALSE;
  return TRUE;
}

/******************************************************************
@@	현재 Line을 지정된 Delta값 만큼 이동시킨다. 
**/
GLine3D GLine3D::Offset(LDOUBLE DeltaX,LDOUBLE DeltaY, LDOUBLE DeltaZ)
{
	GLine3D TempLine;

	TempLine.P1 = P1.Offset(DeltaX, DeltaY, DeltaZ);
	TempLine.P2 = P2.Offset(DeltaX, DeltaY, DeltaZ);

	return TempLine;
}

LDOUBLE GLine3D::Length()
{

	GPoint3D DP;

	DP = P2 - P1;

	return sqrt(DP.x *DP.x + DP.y*DP.y + DP.z*DP.z);
}

LDOUBLE GLine3D::Length2D()
{
  LDOUBLE Dx,Dy;
  Dx  = P2.x - P1.x;
  Dy  = P2.y - P1.y;

  return sqrt(Dx*Dx + Dy*Dy);
}

LDOUBLE GLine3D::GetXDist()
{
  return fabs(P2.x - P1.x);
}

LDOUBLE GLine3D::GetYDist()
{
  return fabs(P2.y - P1.y);
}

LDOUBLE GLine3D::GetZDist()
{
  return fabs(P2.z - P1.z);
}

BOOL GLine3D::UpdatebyMinMax(GPoint3D& Pos)
{
  if(P1.x < P2.x)
  {
    P1.x = min(Pos.x, P1.x);
    P2.x = max(Pos.x, P2.x);
  }
  else // P1.x >= P2.x
  {
    P1.x = max(Pos.x, P1.x);
    P2.x = min(Pos.x, P2.x);
  }

  if(P1.y < P2.y)
  {
    P1.y = min(Pos.y, P1.y);
    P2.y = max(Pos.y, P2.y);
  }
  else // P1.x >= P2.x
  {
    P1.y = max(Pos.y, P1.y);
    P2.y = min(Pos.y, P2.y);
  }

  return TRUE;
}

LDOUBLE GLine3D::GetMinX(LDOUBLE CurMinX)
{
	LDOUBLE MinX;
	MinX = CurMinX;

	MinX = min(P1.x,MinX);
	MinX = min(P2.x,MinX);

	return MinX;
}

LDOUBLE GLine3D::GetMinY(LDOUBLE CurMinY)
{
	LDOUBLE MinY;
	MinY = CurMinY;

	MinY = min(P1.y,MinY);
	MinY = min(P2.y,MinY);

	return MinY;
}

LDOUBLE GLine3D::GetMaxX(LDOUBLE CurMaxX)
{
	LDOUBLE MaxX;
	MaxX = CurMaxX;

	MaxX = max(P1.x,MaxX);
	MaxX = max(P2.x,MaxX);

	return MaxX;
}

LDOUBLE GLine3D::GetMaxY(LDOUBLE CurMaxY)
{
	LDOUBLE MaxY;
	MaxY = CurMaxY;

	MaxY = max(P1.y,MaxY);
	MaxY = max(P2.y,MaxY);

	return MaxY;
}

void GLine3D::Mult(LDOUBLE Factor)
{
  P1 = P1*Factor;
  P2 = P2*Factor;
}

GLine3D& operator * (GLine3D& Src, LDOUBLE MultNum)
{
	static GLine3D Temp;
	Temp.P1 = Src.P1 * MultNum;
	Temp.P2 = Src.P2 * MultNum;

  return Temp;
}

GLine3D& operator * (LDOUBLE MultNum,GLine3D& Src)
{
	static GLine3D Temp;
	
  Temp.P1 = Src.P1 * MultNum;
	Temp.P2 = Src.P2 * MultNum;

	return Temp;
}

void GLine3D::GetCenterPoint(GPoint3D& CenterP)
{
  CenterP.x = (P1.x + P2.x) / 2.;
  CenterP.y = (P1.y + P2.y) / 2.;
  CenterP.z = (P1.z + P2.z) / 2.;
}