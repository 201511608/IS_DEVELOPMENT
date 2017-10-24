#include "stdafx.h"
#include <float.h>
#include <math.h>
#include "LGLtype.h"

#include "C3DPoint.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

C3DPoint::C3DPoint()
{
	x = 0.; 
  y = 0.;
  z = 0.;
}

void C3DPoint::Init()
{
	x = 0.; 
  y = 0.;
  z = 0.;
}

C3DPoint::C3DPoint(double X, double Y, double Z)
{
 	x = X;
  y = Y;
  z = Z;
}

C3DPoint::C3DPoint(C3DPoint& S)
{
	x = S.x;
  y = S.y;
  z = S.z;
}

C3DPoint::~C3DPoint()
{

}

C3DPoint& C3DPoint::operator = (const C3DPoint& Source )
{
	x = Source.x;
  y = Source.y;
  z = Source.z;
	
  return (*this);
}

int C3DPoint::operator == (C3DPoint& point )
{
	if(fabs(x - point.x) >= V_LGLDBL_EPSILON) return 0;
	if(fabs(y - point.y) >= V_LGLDBL_EPSILON) return 0;
	if(fabs(z - point.z) >= V_LGLDBL_EPSILON) return 0;
	
	return 1;
}

int C3DPoint::IsEqual2D(C3DPoint& Point)
{
  if(fabs(x - Point.x) >= V_LGLDBL_EPSILON) return 0;
	if(fabs(y - Point.y) >= V_LGLDBL_EPSILON) return 0;
	
	return 1;
}

int C3DPoint::operator != (C3DPoint& point)
{
	if(	fabs(x - point.x) >= V_LGLDBL_EPSILON ) return 1; 
	if(	fabs(y - point.y) >= V_LGLDBL_EPSILON ) return 1;	
	if(	fabs(z - point.z) >= V_LGLDBL_EPSILON ) return 1;
	
	return 0;
}

C3DPoint& C3DPoint::operator + (C3DPoint& point)
{
	static C3DPoint Temp;
	Temp.x = x + point.x;
	Temp.y = y + point.y;
	Temp.z = z + point.z;
	return Temp;
}

C3DPoint& C3DPoint::operator - (C3DPoint& point)
{
	static C3DPoint Temp;
	Temp.x = x - point.x;
	Temp.y = y - point.y;
	Temp.z = z - point.z;
	return Temp;
}

/**************************************************
@@	Friend Member Function
*/

C3DPoint& operator / (C3DPoint Src,double DivideNum)
{
	static C3DPoint Temp;
	Temp.x = Src.x / DivideNum;
	Temp.y = Src.y / DivideNum;
	Temp.z = Src.z / DivideNum;

	return Temp; 
}


C3DPoint& operator * (C3DPoint Src, double MultNum)
{
	static C3DPoint Temp;
	Temp.x = Src.x * MultNum;
	Temp.y = Src.y * MultNum;
	Temp.z = Src.z * MultNum;
	
	return Temp;
}

C3DPoint& operator * (double MultNum,C3DPoint Src)
{
	static C3DPoint Temp;
	Temp.x = Src.x * MultNum;
	Temp.y = Src.y * MultNum;
	Temp.z = Src.z * MultNum;

	return Temp;
}

void C3DPoint::Set(double X, double Y, double Z)
{
  x = X;
  y = Y;
  z = Z;
}

void C3DPoint::Set(const CPoint& cpoint)
{
	x = cpoint.x;
	y = cpoint.y;
	z = 0.0;
}

CString C3DPoint::GetString(CString strFormat)
{
	CString strReturn;
    strReturn.Format(strFormat,x,y,z);
	return strReturn;
}

C3DPoint& C3DPoint::Offset(double DeltaX, double DeltaY, double DeltaZ)
{
	x = x + DeltaX;
	y = y + DeltaY;
	z = z + DeltaZ;
	return (*this);
}
