#include "stdafx.h"
#define __INC_GRECT_H__
#include "GRenderAll.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
/**********************************************************************
@@	GRect Class Implementation 
@@	Viewing Window�� ���� �Ǽ��� Rectangle�� �ʿ��� ���� ���ȴ�. 
*/
GRect::GRect()
{
	left   = 0.0;
	top    = 0.0;
	right  = 0.0;
	bottom = 0.0;
  depth  = 0.0;
}

GRect::GRect(LDOUBLE Left, LDOUBLE Bottom, LDOUBLE Right, LDOUBLE Top, LDOUBLE Depth/* = 0.0*/)
{
	left   = Left;
	bottom = Bottom;
	right  = Right;
	top    = Top;
  depth  = Depth;
}

// Copy Constructor
GRect::GRect(GRect& SRect)
{
  left   = SRect.left;
  bottom = SRect.bottom;
  right  = SRect.right;
  top    = SRect.top;
  depth  = SRect.depth;
}


GRect::~GRect()
{

}

void GRect::Set(LDOUBLE Left, LDOUBLE Bottom, LDOUBLE Right, LDOUBLE Top, LDOUBLE Depth/* = 0.0*/)
{
	left   = Left;
	bottom = Bottom;
	right  = Right;
	top    = Top;
  depth  = Depth;
}

void GRect::Set(GPoint3D & P1, GPoint3D & P2)
{
  left   = P1.x;
  bottom = P1.y;
  right  = P2.x;
  top    = P2.y;
  depth  = (P1.z + P2.z) /2.0;
  NormalizeRect();
}

void GRect::Set(GLine3D& line3d)
{
  left   = line3d.P1.x;
  bottom = line3d.P1.y;
  right  = line3d.P2.x;
	top    = line3d.P2.y;
  depth  = (line3d.P1.z + line3d.P2.z) /2.0;
  NormalizeRect();
}

void GRect::Move(LDOUBLE DeltaX, LDOUBLE DeltaY)
{
	left   = left   + DeltaX;
	bottom = bottom + DeltaY;
	right  = right  + DeltaX;
	top    = top    + DeltaY;  
}

GRect GRect::operator = (const GRect& Source)
{
	left   = Source.left;
	bottom = Source.bottom;
	right  = Source.right;
	top    = Source.top;
  depth  = Source.depth;
	return (*this);
}

/****************************************************************************************************
@@	������ Delta�� ��ŭ Rectangle�� ��Ǯ����. 
**/

void GRect::InflateRect(LDOUBLE DeltaLeft, LDOUBLE DeltaBottom, LDOUBLE DeltaRight, LDOUBLE DeltaTop)
{
	// ������ �� ��ǥ�迡���� ���� 
	left	=	left   - DeltaLeft;
	bottom	=	bottom - DeltaBottom;
	right	=	right  + DeltaRight;
	top		=	top    + DeltaTop;
}

/****************************************************************************************************
@@	������ Delta�� ��ŭ Rectangle�� ����Ѵ�. 
**/

void GRect::DeflateRect(LDOUBLE DeltaLeft,LDOUBLE DeltaBottom, LDOUBLE DeltaRight, LDOUBLE DeltaTop )
{
	// ������ �� ��ǥ�迡���� ���� 
	left	      =	  left   + DeltaLeft;
	bottom	    =	  bottom + DeltaBottom;
	right		    =	  right  - DeltaRight;
	top			    =	  top    - DeltaTop;
}

LDOUBLE GRect::Width()
{
	return fabs(right - left);
}

LDOUBLE GRect::Height()
{
	return fabs(top - bottom);
}

/*******************************************************************
@@	�־��� GRect Object�� ���� GRect Object�� ������ ���ԵǸ� TRUE��
@@	��ȯ�Ѵ�. 
**/

LBOOL GRect::IsInner(GRect IRect)
{
	NormalizeRect();
	IRect.NormalizeRect();
	if(left <= IRect.left && bottom <= IRect.bottom && 
	   top  >= IRect.top  && right  >= IRect.right ) return TRUE;

	return FALSE;
}

LBOOL GRect::IsIntersect(GRect &Irect)
{
  //return !(left   >= Irect.right || right <= Irect.left||
	//	       bottom >= Irect.top   || top   <= Irect.bottom);
  return !(left   > Irect.right || right < Irect.left||
		       bottom > Irect.top   || top   < Irect.bottom);
}

LBOOL GRect::IsInner(GPoint3D point)
{
	NormalizeRect();
	if(left <= point.x && bottom <= point.y &&
		 top  >= point.y && right  >= point.x   ) return TRUE;
	return FALSE;
}

LBOOL GRect::IsInner(GPoint2D point)
{
	NormalizeRect();
	if(left <= point.x && bottom <= point.y &&
		top >= point.y && right  >= point.x) return TRUE;
	return FALSE;
}

/*****************************************************************************
@@	���� �׽�Ʈ�� �����ϱ� ���� Flag�� �����Ѵ�. 
@@	Line Clipping�˰��򿡼� 
******************************************************************************
@@	B4 | B3 | B2 | B1 
@@	-----------------------------------------------------------------
@@	B1 : ������ ������ Window�� ���ʿ� ���̴� ���         : true
@@	B2 : ������ ������ Window�� �Ʒ��ʿ� ���̴� ���      : true
@@	B3 : ������ ������ Window�� �����ʿ� ���̴� ���      : true
@@	B4 : ������ ������ Window�� ���ʿ� ���̴� ���         : true 
@@  -----------------------------------------------------------------
@@	Bit Field�� �������� �մ��ϳ� ����� Hex Field�� �����Ǿ� �ִ�.
*/
/*
void _SetIntersectFlag(GRect &Window,LINT &P1Flag,LINT &P2Flag,GPoint3D P1,GPoint3D P2)
{
	memset(&P1Flag,0,sizeof(LINT)); // P1Flag �ʱ�ȭ 
	memset(&P2Flag,0,sizeof(LINT)); // P2Flag �ʱ�ȭ  
	
	//	�־��� �� P1�� P2�� ���� Clipping Flag�� �����Ѵ� 
	if(P1.x < Window.left)    
		P1Flag |= 0x1000;				//	B1 : P1�� Window�� ���ʿ� ���̴� ���  
	if(P1.x > Window.right)      
		P1Flag |= 0x0100;				//	B2 : P1�� Window�� �����ʿ� ���̴� ���
	if(P1.y > Window.top)		
		P1Flag |= 0x0001;				//	B3 : P1�� Window�� ���ʿ� ���̴� ���
	if(P1.y < Window.bottom)
		P1Flag |= 0x0010;				//	B4 : P1�� Window�� �Ʒ��ʿ� ���̴� ���  

	
	if(P2.x < Window.left)
		P2Flag |= 0x1000;				//	B1 : P2�� Window�� ���ʿ� ���̴� ���  
	if(P2.x > Window.right)
		P2Flag |= 0x0100;				//	B2 : P2�� Window�� �����ʿ� ���̴� ���
	if(P2.y > Window.top)
		P2Flag |= 0x0001;				//	B3 : P2�� Window�� ���ʿ� ���̴� ���
	if(P2.y < Window.bottom)
		P2Flag |= 0x0010;				//	B4 : P2�� Window�� �Ʒ��ʿ� ���̴� ���  
}


LBOOL GRect::IsIntersect(GraphicObject* pObject)
{
	LINT P1Flag,P2Flag;
	
	GPoint3D *pPoint1, *pPoint2;
	GPoint3D pt1, pt2;
	LINT nGeomType = pObject->WhoAmI();
	pPoint1 = pObject->GetPointForQT(0);
	pt1.x = pPoint1->x;
	pt1.y = pPoint1->y;
	NormalizeRect();
	switch(nGeomType)
	{
	case GOB_POINT:
		if(left <= pt1.x && bottom <= pt1.y &&
			top >= pt1.y && right  >= pt1.x) return TRUE;
		return FALSE;
	case GOB_LINE:
		pPoint2 = pObject->GetPointForQT(1);
		pt2.x = pPoint2->x;
		pt2.y = pPoint2->y;
		_SetIntersectFlag(*this, P1Flag, P2Flag, pt1, pt2);
		if( (P1Flag & P2Flag)  != 0x0000 ) return FALSE;
		return TRUE;
	}
	return TRUE;
}

LBOOL GRect::IsWithin(GraphicObject* pObject)
{
	LINT P1Flag,P2Flag;
	
	GPoint3D *pPoint1, *pPoint2;
	GPoint3D pt1, pt2;
	LINT nGeomType = pObject->WhoAmI();
	pPoint1 = pObject->GetPointForQT(0);
	pt1.x = pPoint1->x;
	pt1.y = pPoint1->y;
	NormalizeRect();
	switch(nGeomType)
	{
	case GOB_POINT:
		if(left <= pt1.x && bottom <= pt1.y &&
			top >= pt1.y && right  >= pt1.x) return TRUE;
		return FALSE;
	case GOB_LINE:
		pPoint2 = pObject->GetPointForQT(1);
		pt2.x = pPoint2->x;
		pt2.y = pPoint2->y;
		_SetIntersectFlag(*this,P1Flag,P2Flag,pt1,pt2);
		if(P1Flag == 0x0000 && P2Flag == 0x0000) return TRUE;
		return FALSE;
	}

	return FALSE;
}
*/
void GRect::NormalizeRect()
{
	LDOUBLE Temp;
	
	if(left > right)
	{
		Temp  = left;
		left  = right;
		right = Temp;
	}

	if(bottom > top)
	{
		Temp   = bottom;
		bottom = top;
		top    = Temp;
	}
}
//
void GRect::Draw(CDC *pdc,COLORREF color)
{
	//int idc;
	//idc = pdc->SaveDC();
	CPen *pOldPen;
	CPen pen;
	pen.CreatePen(PS_SOLID,1,color);
	pOldPen = pdc->SelectObject(&pen);
	
  pdc->MoveTo(GMU::RoundUpInt(left) ,GMU::RoundUpInt(bottom));
  pdc->LineTo(GMU::RoundUpInt(left) ,GMU::RoundUpInt(top));
  pdc->LineTo(GMU::RoundUpInt(right),GMU::RoundUpInt(top));
  pdc->LineTo(GMU::RoundUpInt(right),GMU::RoundUpInt(bottom));
  pdc->LineTo(GMU::RoundUpInt(left) ,GMU::RoundUpInt(bottom));
	
	pdc->SelectObject(pOldPen);
	pen.DeleteObject();
	//pdc->RestoreDC(idc);
}

/**************************************************************************
@@	�־��� Rectangle�� ���̰� �� ���� 1�̵ǵ��� �ϴ� ���� ��ȯ�Ѵ�. 
@@	�� Normalizing Factor�� ��ȯ�Ѵ�.
@@	Viewing ��ǥ����� ��ȯ�� Display�� Model Data�� Projection �۾����� 
@@	Model Data�� ���� Normalized Device Coordinate�� �����ϱ� ���� 
@@	�뵵�� ����� ����.
**/
LDOUBLE GRect::GetNormalizingFactor()
{
	LDOUBLE Width,Height,NorFactor;
	Width	= fabs( right - left );
	Height	= fabs( top - bottom ); 

	NorFactor = 1 / max(Width,Height);

	return NorFactor;
}

CString GRect::GetString(CString Format )
{
	CString StrTemp;
	StrTemp.Format(Format,left,bottom,right,top);
	return StrTemp;
}

int GRect::operator == (const GRect rect)
{
	NormalizeRect();

	if(fabs(left   - rect.left)			>= V_LGLDBL_EPSILON) return 0;
	if(fabs(right  - rect.right)		>= V_LGLDBL_EPSILON) return 0;
	if(fabs(top    - rect.top)			>= V_LGLDBL_EPSILON) return 0;
	if(fabs(bottom - rect.bottom)	  >= V_LGLDBL_EPSILON) return 0;
  if(fabs(depth  - rect.depth)    >= V_LGLDBL_EPSILON) return 0;
	
	return 1;
}

int GRect::operator != (const GRect rect)
{
	NormalizeRect();

	if(fabs(left   - rect.left)   >= V_LGLDBL_EPSILON) return 1;
	if(fabs(right  - rect.right)  >= V_LGLDBL_EPSILON) return 1;
	if(fabs(top    - rect.top)    >= V_LGLDBL_EPSILON) return 1;
	if(fabs(bottom - rect.bottom) >= V_LGLDBL_EPSILON) return 1;
  if(fabs(depth  - rect.depth)  >= V_LGLDBL_EPSILON) return 1;
	
	return 0;
}

GPoint3D __tempP3D;


GPoint3D& GRect::GetLB()
{
  __tempP3D.Set(left,bottom,depth);
  return __tempP3D;
}

GPoint3D& GRect::GetLT()
{
  __tempP3D.Set(left,top,depth);
  return __tempP3D;
}


GPoint3D& GRect::GetRT()
{
  __tempP3D.Set(right,top,depth);
  return __tempP3D;
}

GPoint3D& GRect::GetRB()
{
  __tempP3D.Set(right,bottom,depth);
  return __tempP3D;
}



GPoint3D& GRect::GetCenterP3D(LDOUBLE ZCenter /* = 0. */)
{
  __tempP3D.Set((left+right)/2.,(bottom+top)/2.,ZCenter);
  return __tempP3D;
}

//////////////////////////////////////////////
// NDC���� Window Rect�� depth�� View��ǥ�� ������ NDC��ȯ�� View Plane�� 
// Z ��ǥ�� �ǹ��Ѵ�. 

