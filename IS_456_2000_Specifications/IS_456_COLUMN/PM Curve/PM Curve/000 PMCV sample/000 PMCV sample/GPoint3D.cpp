#include "stdafx.h"
#include <math.h>
#define __INC_GPOINT3D_H__
#include "GRenderAll.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

GPoint3D::GPoint3D()
{
	x = 0.;
  y = 0.;
  z = 0.;
}

void GPoint3D::Init()
{
	x = 0.;
  y = 0.;
  z = 0.;
}

GPoint3D::GPoint3D(LDOUBLE X, LDOUBLE Y, LDOUBLE Z)
{
 	x = X;
  y = Y;
  z = Z;
}

GPoint3D::GPoint3D(GPoint3D& S)
{
	x = S.x;
  y = S.y;
  z = S.z;
}

GPoint3D::GPoint3D(C3DPoint& P)
{
  x = P.x;
  y = P.y;
  z = P.z;
}

GPoint3D::~GPoint3D()
{

}

void GPoint3D::CopyPoint3D(GPoint3D Source)
{
	x = Source.x;
  y = Source.y;
  z = Source.z;
}

GPoint3D& GPoint3D::operator = (const GPoint3D& Source )
{
	x = Source.x;
  y = Source.y;
  z = Source.z;
	return (*this);
}

GPoint3D& GPoint3D::operator = (const C3DPoint& Source)
{
	x = Source.x;
  y = Source.y;
  z = Source.z;
	return (*this);
}


int GPoint3D::operator == (GPoint3D& point )
{
	if(fabs(x - point.x) >= V_LGLDBL_EPSILON) return 0;
	if(fabs(y - point.y) >= V_LGLDBL_EPSILON) return 0;
	if(fabs(z - point.z) >= V_LGLDBL_EPSILON) return 0;
	
	return 1;
}

int GPoint3D::IsSamePoint(GPoint3D& point , double Tol)
{
  if(fabs(x - point.x) >= fabs(Tol)) return 0;
	if(fabs(y - point.y) >= fabs(Tol)) return 0;
	if(fabs(z - point.z) >= fabs(Tol)) return 0;
  return 1;
}

int GPoint3D::IsEqual2D(GPoint3D& Point)
{
  if(fabs(x - Point.x) >= V_LGLDBL_EPSILON) return 0;
	if(fabs(y - Point.y) >= V_LGLDBL_EPSILON) return 0;
	
	return 1;
}


int GPoint3D::operator != (GPoint3D& point)
{
	if(	fabs(x - point.x) >= V_LGLDBL_EPSILON ) return 1; 
	if(	fabs(y - point.y) >= V_LGLDBL_EPSILON ) return 1;	
	if(	fabs(z - point.z) >= V_LGLDBL_EPSILON ) return 1;
	
	return 0;
}

GPoint3D& GPoint3D::operator + (GPoint3D& point)
{
	static GPoint3D Temp;
	Temp.x = x + point.x;
	Temp.y = y + point.y;
	Temp.z = z + point.z;
	return Temp;
}

GPoint3D& GPoint3D::operator - (GPoint3D& point)
{
	static GPoint3D Temp;
	Temp.x = x - point.x;
	Temp.y = y - point.y;
	Temp.z = z - point.z;
	return Temp;
}

/**************************************************
@@	Friend Member Function
*/

GPoint3D& operator / (GPoint3D& Src,LDOUBLE DivideNum)
{
	static GPoint3D Temp;
	Temp.x = Src.x / DivideNum;
	Temp.y = Src.y / DivideNum;
	Temp.z = Src.z / DivideNum;

	return Temp; 
}


GPoint3D& operator * (GPoint3D& Src, LDOUBLE MultNum)
{
	static GPoint3D Temp;
	Temp.x = Src.x * MultNum;
	Temp.y = Src.y * MultNum;
	Temp.z = Src.z * MultNum;
	
	return Temp;
}

GPoint3D& operator * (LDOUBLE MultNum,GPoint3D& Src)
{
	static GPoint3D Temp;
	Temp.x = Src.x * MultNum;
	Temp.y = Src.y * MultNum;
	Temp.z = Src.z * MultNum;

	return Temp;
}


GPoint3D& GPoint3D::operator + ( C3DPoint& point)
{
	static GPoint3D Temp;
	Temp.x = x + point.x;
	Temp.y = y + point.y;
	Temp.z = z + point.z;
	return Temp;
}

void GPoint3D::Set(LDOUBLE X, LDOUBLE Y, LDOUBLE Z)
{
  x = X;
  y = Y;
  z = Z;
}

void GPoint3D::Set(const GPoint2D& point2d)
{
	x = point2d.x;
	y = point2d.y;
	z = 0.0;
}

void GPoint3D::Set(const CPoint& cpoint)
{
	x = cpoint.x;
	y = cpoint.y;
	z = 0.0;
}

void GPoint3D::Set(C3DPoint& c3dpoint)
{
  x = c3dpoint.x;
  y = c3dpoint.y;
  z = c3dpoint.z;
}

CString GPoint3D::GetString(CString strFormat)
{
	CString strReturn;
    strReturn.Format(strFormat,x,y,z);
	return strReturn;
}

GPoint3D& GPoint3D::Offset(LDOUBLE DeltaX, LDOUBLE DeltaY, LDOUBLE DeltaZ)
{
	//static GPoint3D TempPoint;
	//TempPoint.x = x + DeltaX;
	//TempPoint.y = y + DeltaY;
	//TempPoint.z = z + DeltaZ;
	//return TempPoint;
	
	x = x + DeltaX;
	y = y + DeltaY;
	z = z + DeltaZ;

	return (*this);
	
}

// 설정된 x,y,z 값을 Vector 성분으로 간주 Vector Cross product연산을 수행한다. 
GPoint3D GPoint3D::Cross(GPoint3D vec)
{
	GPoint3D Temp;
     
  Temp.x = y*vec.z-z*vec.y;
  Temp.y = z*vec.x-x*vec.z;
  Temp.z = x*vec.y-y*vec.x;

  return Temp;
} 
// 설정된 x,y,z 값을 Vector 성분으로 간주 Vector Dot product 연산을 수행한다. 
LDOUBLE GPoint3D::Dot(GPoint3D vec)		// dot product
{
	GPoint3D Temp;

  Temp.x = x*vec.x;
  Temp.y = y*vec.y;
  Temp.z = z*vec.z;
	
	LDOUBLE ret = Temp.x+Temp.y+Temp.z;
  return  ret;
}