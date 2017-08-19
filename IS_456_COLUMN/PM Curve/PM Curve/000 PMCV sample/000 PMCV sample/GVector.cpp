#include "stdafx.h"
#define __INC_GVECTOR_H__
#include "GRenderAll.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

GVector::GVector()     //dummy
{

}

GVector::GVector( LDOUBLE xi,LDOUBLE yi,LDOUBLE zi)
{
	
  x = xi;
  y = yi;
  z = zi;
	
	//memcpy(&x, &xi,sizeof(LDOUBLE));
	//memcpy(&y, &yi,sizeof(LDOUBLE));
	//memcpy(&z, &zi,sizeof(LDOUBLE));
}

GVector::GVector( GPoint3D StartPoint,GPoint3D EndPoint)
{
	x = EndPoint.x - StartPoint.x;
	y = EndPoint.y - StartPoint.y;
	z = EndPoint.z - StartPoint.z;
}

GVector::GVector( GPoint3D& VectorPoint)
{
	x = VectorPoint.x;
	y = VectorPoint.y;
	z = VectorPoint.z;
}

void GVector::ZeroVector()
{
	x = 0.;
  y = 0.;
  z = 0.;
  //memset(&x,0,sizeof(LDOUBLE));
	//memset(&y,0,sizeof(LDOUBLE));
	//memset(&z,0,sizeof(LDOUBLE));
}

int GVector::operator == (const GVector point )
{
	if(fabs(x - point.x) >= V_LGLDBL_EPSILON) return 0;
	if(fabs(y - point.y) >= V_LGLDBL_EPSILON) return 0;
	if(fabs(z - point.z) >= V_LGLDBL_EPSILON) return 0;
	
	return 1;
}

int GVector::operator != (const GVector point)
{
	if(fabs(x - point.x) >= V_LGLDBL_EPSILON) return 1;
	if(fabs(y - point.y) >= V_LGLDBL_EPSILON) return 1;
	if(fabs(z - point.z) >= V_LGLDBL_EPSILON) return 1;
	
	return 0;
}

GVector GVector::operator+(GVector vec)
{
	GVector Temp;

  Temp.x = x+vec.x;
  Temp.y = y+vec.y;
  Temp.z = z+vec.z;

  return Temp;
}

GVector GVector::operator-(GVector vec)
{
	GVector Temp;

  Temp.x = x-vec.x;
  Temp.y = y-vec.y;
  Temp.z = z-vec.z;

  return Temp;
}

GVector GVector::operator = (GVector src)
{
	x = src.x;
	y = src.y;
	z = src.z;
	
  return *this;
}

LBOOL GVector::IsAcuteAngle(GVector& OVec)
{
  LDOUBLE DotD = Dot(OVec);
  if(DotD > 0.) return TRUE;
  return FALSE;
}

LBOOL GVector::IsObtuseAngle(GVector& OVec)
{
  LDOUBLE DotD = Dot(OVec);
  if(DotD > 0.) return FALSE;
  return TRUE;
}


GVector GVector::operator*(GVector vec)  // cross operator
{
	GVector Temp;
     
  Temp.x = y*vec.z-z*vec.y;
  Temp.y = z*vec.x-x*vec.z;
  Temp.z = x*vec.y-y*vec.x;

  return Temp;
}

GVector GVector::Cross(GVector vec)
{
	GVector Temp;
     
  Temp.x = y*vec.z-z*vec.y;
  Temp.y = z*vec.x-x*vec.z;
  Temp.z = x*vec.y-y*vec.x;

  return Temp;
} 

LDOUBLE GVector::Dot(GVector vec)		// dot product
{
	GVector Temp;

  Temp.x = x*vec.x;
  Temp.y = y*vec.y;
  Temp.z = z*vec.z;
	
	LDOUBLE ret = Temp.x+Temp.y+Temp.z;
  return  ret;
}

// Radian 
BOOL GVector::Angle3D(GVector AVector, double & Angle)
{
  double VDot = Dot(AVector); 
  double V1Leng, V2Leng;
  double Temp,CosV;
  V1Leng = Length();
  V2Leng = AVector.Length();
  
  if(GMU::IsZero(V1Leng)) return FALSE;
  if(GMU::IsZero(V2Leng)) return FALSE;
  
  Temp = V1Leng* V2Leng;

  CosV = VDot / Temp;

  // _matherror를 overide하면 될텐데 .. Release Mode는 설정된것으로 
  // 이야기 들었다.DebugMode는 그냥 두자 
#ifdef _DEBUG
  if(CosV < -1)
  {
    //ASSERT(0);
    CosV = -1;
  }
  
  if(CosV > 1 )
  {
    //ASSERT(0);
    CosV = 1 ;
  }
#endif
  Angle = acos(CosV);

  return TRUE;
}

BOOL GVector::Cosine3D(GVector AVector, double & CosVal)
{
  double VDot = Dot(AVector); 
  double V1Leng, V2Leng;
  double Temp;
  V1Leng = Length();
  V2Leng = AVector.Length();
  
  if(GMU::IsZero(V1Leng)) return FALSE;
  if(GMU::IsZero(V2Leng)) return FALSE;
  
  Temp = V1Leng* V2Leng;

  CosVal = VDot / Temp;

  return TRUE;
}

LDOUBLE GVector::CCWAngle2D() // must be (z == 0.)
{
  double VLength = Length();
  if(GMU::IsZero(VLength)) return 0.;

  if(GMU::IsZero(x)) 
  {
    if(y > 0.) return  M_PI   /2;  // 90  Degree
    else       return (M_PI*3)/2;  // 270 Degree
  }

  if(GMU::IsZero(y))
  {
    if(x > 0.) return 0.;         // 0   Degree
    else       return M_PI;       // 180 Degree
  }
  
  MakeUnit();

  if( x > 0.)
  {
    if(y > 0.) // 1사분면 
    {
      return acos(x);  
    }
    else       // 4사분면 
    {
      double Theta = acos(x);
      return ((M_PI*2) - Theta);
    }
  }
  else
  {
    if(y > 0.) // 2사분면 
    {
      return acos(x);
    }
    else       // 3사분면 
    {
      double Theta = acos(x);
      return ((M_PI) - Theta); 
    }
  }

  return 0.;
}

void GVector::Set( LDOUBLE xi,LDOUBLE yi,LDOUBLE zi)
{
	memcpy(&x, &xi,sizeof(LDOUBLE));
	memcpy(&y, &yi,sizeof(LDOUBLE));
	memcpy(&z, &zi,sizeof(LDOUBLE));
}

void GVector::Set( GPoint3D StartPoint,GPoint3D EndPoint)
{
	x = EndPoint.x - StartPoint.x;
	y = EndPoint.y - StartPoint.y;
	z = EndPoint.z - StartPoint.z;
}

void GVector::Set( GPoint3D VectorPoint)
{
	x = VectorPoint.x;
	y = VectorPoint.y;
	z = VectorPoint.z;
}

GVector GVector::Unit()
{
	LDOUBLE length,TAbs;
	GVector a(0,0,0);

  TAbs = Abs();
  
  if ( TAbs!=0.0 )
  {
    length = 1.0 / TAbs;
    a.x = x*length;
    a.y = y*length;
    a.z = z*length;
  }
  
  return a;
}

GVector GVector::MakeUnit()
{
	LDOUBLE length,TAbs;
	
	TAbs = Abs();
  if ( TAbs!=0.0 )
  {
    length = 1.0 / TAbs;
		x = x*length;
		y = y*length;
		z = z*length;
	}
	return *this;
}

GVector GVector::Scaling( LDOUBLE r)
{
	GVector Temp;

  Temp.x = x*r;
  Temp.y = y*r;
  Temp.z = z*r;

  return Temp;
}

void GVector::ScalingMySelf(LDOUBLE r)
{
  x = x*r;
  y = y*r;
  z = z*r;
}


GVector GVector::Trans( GVector t)
{
	GVector Temp;

	Temp.x = x+t.x;
	Temp.y = y+t.y;
	Temp.z = z+t.z;

	return Temp;
}

GVector GVector::Trans( LDOUBLE xi,LDOUBLE yi,LDOUBLE zi)
{
  GVector Temp;

  Temp.x = x+xi;
  Temp.y = y+yi;
  Temp.z = z+zi;

  return Temp;
}

LDOUBLE GVector::Abs()
{
  LDOUBLE ret = x*x+y*y+z*z;
  ret = sqrt(ret);

  return ret;
}

LDOUBLE GVector::Length()
{
	LDOUBLE ret = x*x+y*y+z*z;
	ret = sqrt(ret);

  return ret;
}


CString GVector::GetString(CString strFormat)
{
	CString strReturn;
  strReturn.Format(strFormat,x,y,z);
	
  return strReturn;
}

LBOOL GVector::MakeNormalVector(GPoint3D& P1, GPoint3D& P2, GPoint3D& P3)
{
  GVector TXVector,TYVector;
  TXVector.Set(P1,P2);TXVector.MakeUnit();
  TYVector.Set(P1,P3);TYVector.MakeUnit();
  *this = TXVector * TYVector;
    
  return IsValid();
}

LBOOL GVector::MakeNormalVector(GPoint3D P1, GPoint3D P2, GPoint3D P3, GVector DirectionV)
{
  if(MakeNormalVector(P1,P2,P3))
  {
    if(Dot(DirectionV) < 0. )
      Inverse();
    return TRUE;
  }
  return FALSE;
}

LBOOL GVector::MakeNormalVector(GPoint3D P1, GPoint3D P2, GPoint3D P3, GPoint3D DirectionP)
{
  GVector DVector;
  
  if(MakeNormalVector(P1,P2,P3))
  {
    DVector.Set(P1,DirectionP);
    if(Dot(DVector) < 0. )
      Inverse();
    return TRUE;
  }
  return FALSE;  
}

// Zero Vector이면 FALSE....
LBOOL GVector::IsValid()
{
  if(GMU::IsZero(Length())) return FALSE;
  return TRUE;
}


void GVector::Inverse()
{
  x = -x;
  y = -y;
  z = -z;
}

