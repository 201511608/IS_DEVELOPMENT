#ifndef __GLINE3D_H__
#define __GLINE3D_H__

/**************************************************
@@	Viewing��ȯ�� ���� �ӽ÷� ���Ǵ� Line ����ü 
@@	Modeling Data�ʹ� ������ �ڷ���. 
**/
class GLine2D;
class GLine3D
{
public:
	int m_YIndex; // SDS Mesh Line ������ ���� �߰���...
  int m_XIndex; // SDS Mesh Line ������ ���� �߰���...
  GPoint3D P1;
	GPoint3D P2;
public:
	GLine3D();
  GLine3D(GPoint3D&p1, GPoint3D&p2);
  GLine3D(GLine3D& S);

	virtual ~GLine3D();
	GLine3D& operator = (const GLine3D& Source);
	void Set(const GPoint3D& p1,const GPoint3D& p2);
	void Set(const GLine2D& Line2D);
	
	CString GetString(CString Format = "(%f, %f, %f), (%f, %f, %f)");
	
	void SetZero();
	LBOOL IsValid();
  LBOOL IsValid2D();

	GLine3D Offset(LDOUBLE DeltaX, LDOUBLE DeltaY, LDOUBLE DeltaZ);
	void Mult(LDOUBLE Factor);
  BOOL UpdatebyMinMax(GPoint3D& Pos);
  /************************************************************
	@@	Line�� ���̸� ��ȯ�Ѵ�.
	**/
	LDOUBLE Length();
  LDOUBLE Length2D();
	/************************************************************
	@@	Line�� Point�鿡�� MinX,MinY,MaxX,MaxY�� �����Ѵ�.
	@@	Argument�� ������ MinX,MinY,MaxX,MaxY�� �ǹ��Ѵ�. 
	**/
	LDOUBLE GetMinX(LDOUBLE CurMinX);
	LDOUBLE GetMinY(LDOUBLE CurMinY);
	LDOUBLE GetMaxX(LDOUBLE CurMaxX);
	LDOUBLE GetMaxY(LDOUBLE CurMaxY);

  LDOUBLE GetXDist();
  LDOUBLE GetYDist();
  LDOUBLE GetZDist();
  void GetCenterPoint(GPoint3D& CenterP);

  friend GLine3D& operator * (GLine3D& Src,LDOUBLE MultNum);
	friend GLine3D& operator * (LDOUBLE MultNum, GLine3D& Src);
};


#endif