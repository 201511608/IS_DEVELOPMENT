#ifndef __GRECT_H__
#define __GRECT_H__

/***********************************************
@@	LDOUBLE형 Rectangle 구조체 
@@	오른손 실 좌표계에서의 Rectangle 정의 
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
	@@	주어진 GRect Object가 현재 GRect Object에 완전히 포함되면 TRUE를
	@@	반환한다. 
	**/
	LBOOL IsInner(GRect IRect);
	/*******************************************************************
	@@	주어진 Point의 X Y 좌표가 현재 GRect Object에 완전히 포함되면 TRUE를
	@@	반환한다. 
	**/
	LBOOL IsInner(GPoint3D Point);
	LBOOL IsInner(GPoint2D Point);
	/************************************************************************************************
	@@	지정된 Delta값 만큼 이동 시킨다. 
	**/
	void Move(LDOUBLE DeltaX, LDOUBLE DeltaY);
	/***************************************************************************
	@@	지정된 Delta값 만큼  부풀린다 
	@@	Argument로 -값이 전달되면 그 반대의 처리가 행해진다. 
	**/
	void	InflateRect(LDOUBLE DeltaLeft,LDOUBLE DeltaBottom, LDOUBLE DeltaRight, LDOUBLE DeltaTop );
	/***************************************************************************
	@@	지정된 Delta값 만큼  줄인다. 
	@@	Argument로 -값이 전달되면 그 반대의 처리가 행해진다. 
	**/
	void	DeflateRect(LDOUBLE DeltaLeft, LDOUBLE DeltaBottom, LDOUBLE DeltaRight, LDOUBLE DeltaTop);
	LDOUBLE Width();
	LDOUBLE Height();
	/**********************************************************************
	@@	주어진 Rectangle을 길이가 긴 변이 1이되도록 하는 값을 반환한다. 
	@@	즉 Normalizing Factor를 반환한다.
	**/
	LDOUBLE GetNormalizingFactor();
	/**********************************************************************
	@@	주어진 Color로 Rectangle을 그린다. 
	**/
	void Draw(CDC *pdc,COLORREF color);
	/**********************************************************************
	@@	left > right 또는 bottom > top 조건등의 일반적인 Rectangle의 조건을 
	@@	만족하지 못하는 Data가 설정되었을때 이러한 Data를 표준화 시킨다.
	**/
	void NormalizeRect();
	
	CString GetString(CString Format ="(L: %f B: %f R: %f T: %f)");

	/*****************************************************************************
	@@	제시된 Rectangle에  걸치는 Line에 대해 TRUE를 반환한다.
	**/
	//LBOOL IsIntersect(GOB_Line Line);

	/******************************************************************************
	@@	Line이 완전히 포함될때 TRUE를 반환한다.
	@@	제시되는 Line은 2D Line Object임을 가정한다. 
	**/
	//LBOOL IsWithin(GOB_Line Line);
	
	/*********************************************
	@@ GGeometryEngine으로 옮길 예정 
	**/
	//LBOOL IsWithin(GraphicObject* pObject);
	//LBOOL IsIntersect(GraphicObject* pObject);

  GPoint3D& GetCenterP3D(LDOUBLE ZCenter  = 0.);
};

#endif