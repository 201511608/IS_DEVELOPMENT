#ifndef __GRECT_H__
#define __GRECT_H__

/***********************************************
@@	LDOUBLE�� Rectangle ����ü 
@@	������ �� ��ǥ�迡���� Rectangle ���� 
**/
class GRect
{
public:
	GPoint3D& GetRB();
	GPoint3D& GetRT();
	GPoint3D& GetLT();
	GPoint3D& GetLB();
	LBOOL IsIntersect(GRect& Irect);
	
	LDOUBLE left;
	LDOUBLE bottom;
	LDOUBLE right;
	LDOUBLE top;
  LDOUBLE depth; // Z Depth
	
	GRect();
	GRect(LDOUBLE Left, LDOUBLE Bottom, LDOUBLE Right, LDOUBLE Top, LDOUBLE Depth = 0.0);
	GRect(GRect& SRect);
  virtual ~GRect();

  void	Set(LDOUBLE Left, LDOUBLE Bottom, LDOUBLE Right, LDOUBLE Top, LDOUBLE Depth = 0.0 );
	void  Set(GPoint3D & P1, GPoint3D & P2);
  void  Set(GLine3D &line3d);
  GRect	operator = (const GRect& Source);
	int		operator == (const GRect rect);
	int		operator != (const GRect rect);
	
	/*******************************************************************
	@@	�־��� GRect Object�� ���� GRect Object�� ������ ���ԵǸ� TRUE��
	@@	��ȯ�Ѵ�. 
	**/
	LBOOL IsInner(GRect IRect);
	/*******************************************************************
	@@	�־��� Point�� X Y ��ǥ�� ���� GRect Object�� ������ ���ԵǸ� TRUE��
	@@	��ȯ�Ѵ�. 
	**/
	LBOOL IsInner(GPoint3D Point);
	LBOOL IsInner(GPoint2D Point);
	/************************************************************************************************
	@@	������ Delta�� ��ŭ �̵� ��Ų��. 
	**/
	void Move(LDOUBLE DeltaX, LDOUBLE DeltaY);
	/***************************************************************************
	@@	������ Delta�� ��ŭ  ��Ǯ���� 
	@@	Argument�� -���� ���޵Ǹ� �� �ݴ��� ó���� ��������. 
	**/
	void	InflateRect(LDOUBLE DeltaLeft,LDOUBLE DeltaBottom, LDOUBLE DeltaRight, LDOUBLE DeltaTop );
	/***************************************************************************
	@@	������ Delta�� ��ŭ  ���δ�. 
	@@	Argument�� -���� ���޵Ǹ� �� �ݴ��� ó���� ��������. 
	**/
	void	DeflateRect(LDOUBLE DeltaLeft, LDOUBLE DeltaBottom, LDOUBLE DeltaRight, LDOUBLE DeltaTop);
	LDOUBLE Width();
	LDOUBLE Height();
	/**********************************************************************
	@@	�־��� Rectangle�� ���̰� �� ���� 1�̵ǵ��� �ϴ� ���� ��ȯ�Ѵ�. 
	@@	�� Normalizing Factor�� ��ȯ�Ѵ�.
	**/
	LDOUBLE GetNormalizingFactor();
	/**********************************************************************
	@@	�־��� Color�� Rectangle�� �׸���. 
	**/
	void Draw(CDC *pdc,COLORREF color);
	/**********************************************************************
	@@	left > right �Ǵ� bottom > top ���ǵ��� �Ϲ����� Rectangle�� ������ 
	@@	�������� ���ϴ� Data�� �����Ǿ����� �̷��� Data�� ǥ��ȭ ��Ų��.
	**/
	void NormalizeRect();
	
	CString GetString(CString Format ="(L: %f B: %f R: %f T: %f)");

	/*****************************************************************************
	@@	���õ� Rectangle��  ��ġ�� Line�� ���� TRUE�� ��ȯ�Ѵ�.
	**/
	//LBOOL IsIntersect(GOB_Line Line);

	/******************************************************************************
	@@	Line�� ������ ���Եɶ� TRUE�� ��ȯ�Ѵ�.
	@@	���õǴ� Line�� 2D Line Object���� �����Ѵ�. 
	**/
	//LBOOL IsWithin(GOB_Line Line);
	
	/*********************************************
	@@ GGeometryEngine���� �ű� ���� 
	**/
	//LBOOL IsWithin(GraphicObject* pObject);
	//LBOOL IsIntersect(GraphicObject* pObject);

  GPoint3D& GetCenterP3D(LDOUBLE ZCenter  = 0.);
};

#endif