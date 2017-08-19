#ifndef __GPOINT2D_H__
#define __GPOINT2D_H__

class GPoint3D;
class GPoint2D
{
public:
	LINT x,y;	
public:
	GPoint2D();
	GPoint2D(LINT X,LINT Y);
	GPoint2D(GPoint2D &S); // Copy Constructor
	virtual ~GPoint2D();
public:
	void DrawEllipse(CDC* pDC,int Width,int Height);
	int operator == (const GPoint2D point);
	int operator != (const GPoint2D point);

	GPoint2D& operator = (const GPoint2D& Source);
	GPoint2D	operator  + (const GPoint2D point);
	GPoint2D operator  - (const GPoint2D point);

	friend GPoint2D operator / (GPoint2D Src,LINT DivideNum);
	friend GPoint2D operator * (GPoint2D Src,LINT MultNum);
	friend GPoint2D operator * (LINT MultNum, GPoint2D Src);
	
	void Set(GPoint3D& point3d);
	void Set(LINT X, LINT Y);
	CString GetString(CString strFormat="(%d,%d)");
	void Init();
	void Draw(CDC *pdc,COLORREF Color);
  void DrawLarge(CDC *pdc,COLORREF Color);
  void DrawCross(CDC* pDC, int Length);
};

//extern LINT RoundUpInt(LDOUBLE dVal);
//extern LINT RoundOffInt(LDOUBLE dVal);

#endif