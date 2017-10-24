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
	@@	Line�� ���� ��ǥ���� X��� �̷�� ������ ��ȯ�Ѵ�.
	@@  ��ȯ���� Radian���̴�.	
	@@	X�࿡ ������ ��쿡�� Y��ǥ�� ���ؼ� ������ ��ȣ�� ����
	@@	�Ѵ�. 
	**/
	LDOUBLE GetAngle();
	
	LINT Length();
	void DrawStyled(CDC *pdc,COLORREF color,int nStyle);
  void Draw(CDC *pdc);
	void Draw(CDC *pdc,COLORREF color);
  void Draw(CDC* pdc, COLORREF color, int LineWidth);
	/************************************************************
	@@	Line�� Point�鿡�� MinX,MinY,MaxX,MaxY�� �����Ѵ�.
	@@	Argument�� ������ MinX,MinY,MaxX,MaxY�� �ǹ��Ѵ�. 
	**/
	LINT GetMinX(LINT CurMinX);
	LINT GetMinY(LINT CurMinY);
	LINT GetMaxX(LINT CurMaxX);
	LINT GetMaxY(LINT CurMaxY);

	int operator == ( GLine2D line);
	int operator != ( GLine2D line);
};

#endif