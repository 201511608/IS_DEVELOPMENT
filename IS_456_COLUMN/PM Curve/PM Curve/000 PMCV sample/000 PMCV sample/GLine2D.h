#ifndef __GLINE2D_H__
#define __GLINE2D_H__

class GLine3D;
class GLine2D
{
public:
	GLine2D();
	virtual ~GLine2D();
	GPoint2D P1;
	GPoint2D P2;
	GLine2D operator = (GLine2D Source);
	void Set(GPoint2D p1,GPoint2D p2);
	void Set(GLine3D& Line3D);
	CString GetString();
	void SetZero();
	LBOOL IsValid();
	/*************************************************************
	@@	Line이 현재 좌표계의 X축과 이루는 각도를 반환한다.
	@@  반환값은 Radian값이다.	
	@@	X축에 수직일 경우에는 Y좌표를 비교해서 각도의 부호를 결정
	@@	한다. 
	**/
	LDOUBLE GetAngle();
	
	LINT Length();
	void DrawStyled(CDC *pdc,COLORREF color,int nStyle);
  void Draw(CDC *pdc);
	void Draw(CDC *pdc,COLORREF color);
  void Draw(CDC* pdc, COLORREF color, int LineWidth);
	/************************************************************
	@@	Line의 Point들에서 MinX,MinY,MaxX,MaxY를 산출한다.
	@@	Argument는 현재의 MinX,MinY,MaxX,MaxY를 의미한다. 
	**/
	LINT GetMinX(LINT CurMinX);
	LINT GetMinY(LINT CurMinY);
	LINT GetMaxX(LINT CurMaxX);
	LINT GetMaxY(LINT CurMaxY);

	int operator == ( GLine2D line);
	int operator != ( GLine2D line);
};

#endif