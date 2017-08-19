#ifndef __C3DPOINT_H__
#define __C3DPOINT_H__

class C3DPoint
{
public:
	double x,y,z;	
public:
	C3DPoint();
	C3DPoint(double X, double Y, double Z = 0);
	C3DPoint(C3DPoint &S); // Copy Constructor
	virtual ~C3DPoint();
public:
	int IsEqual2D(C3DPoint& Point);
  int operator == (C3DPoint& point);
	int operator != (C3DPoint& point);
	
	C3DPoint& operator = (const C3DPoint& Source);
	C3DPoint& operator + (C3DPoint& point);
	C3DPoint& operator - (C3DPoint& point);
	friend C3DPoint& operator / (C3DPoint Src,double DivideNum);
	friend C3DPoint& operator * (C3DPoint Src,double MultNum);
	friend C3DPoint& operator * (double MultNum, C3DPoint Src);
	void Set(double X, double Y, double Z = 0);
	void Set(const CPoint& cpoint);
	C3DPoint& Offset(double DeltaX, double DeltaY, double DeltaZ);
	CString GetString(CString strFormat="(%f,%f,%f)");
	void Init();
};

#endif