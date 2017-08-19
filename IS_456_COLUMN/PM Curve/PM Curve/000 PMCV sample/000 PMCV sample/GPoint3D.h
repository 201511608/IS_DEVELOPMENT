#ifndef __GPOINT3D_H__
#define __GPOINT3D_H__

#include "C3DPoint.h"

class GPoint2D;
class GPoint3D
{
public:
	LDOUBLE x,y,z;	
public:
	GPoint3D();
	GPoint3D(LDOUBLE X, LDOUBLE Y, LDOUBLE Z = 0);
	GPoint3D(GPoint3D &S); // Copy Constructor
	GPoint3D(C3DPoint &P);
  virtual ~GPoint3D();
public:
	int IsSamePoint(GPoint3D& point , double Tol);
	int IsEqual2D(GPoint3D& Point);
  int operator == (GPoint3D& point);
	int operator != (GPoint3D& point);

	void CopyPoint3D(GPoint3D Source);
	GPoint3D& operator = (const GPoint3D& Source);
  GPoint3D& operator = (const C3DPoint& Source);
	GPoint3D& operator + (GPoint3D& point);
	GPoint3D& operator - (GPoint3D& point);
	friend GPoint3D& operator / (GPoint3D& Src,LDOUBLE DivideNum);
	friend GPoint3D& operator * (GPoint3D& Src,LDOUBLE MultNum);
	friend GPoint3D& operator * (LDOUBLE MultNum, GPoint3D& Src);
  //friend GPoint3D& operator + (GPoint3D& Src,C3DPoint& point);
  GPoint3D& operator + (C3DPoint& point);
	void Set(LDOUBLE X, LDOUBLE Y, LDOUBLE Z = 0);
	void Set(const GPoint2D& point2d);
	void Set(const CPoint& cpoint);
  void Set(C3DPoint& c3dpoint);
  GPoint3D Cross(GPoint3D vec);
  LDOUBLE  Dot(GPoint3D vec)	;	// dot product
	GPoint3D& Offset(LDOUBLE DeltaX, LDOUBLE DeltaY, LDOUBLE DeltaZ);
	CString GetString(CString strFormat="(%f,%f,%f)");
	void Init();
};

#endif