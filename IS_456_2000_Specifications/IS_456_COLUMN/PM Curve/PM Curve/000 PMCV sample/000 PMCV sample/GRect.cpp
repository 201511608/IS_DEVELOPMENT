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
@@	Viewing Window등 각종 실수형 Rectangle이 필요한 곳에 사용된다. 
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
@@	지정된 Delta값 만큼 Rectangle을 부풀린다. 
**/

void GRect::InflateRect(LDOUBLE DeltaLeft, LDOUBLE DeltaBottom, LDOUBLE DeltaRight, LDOUBLE DeltaTop)
{
	// 오른손 실 좌표계에서의 정의 
	left	=	left   - DeltaLeft;
	bottom	=	bottom - DeltaBottom;
	right	=	right  + DeltaRight;
	top		=	top    + DeltaTop;
}

/****************************************************************************************************
@@	지정된 Delta값 만큼 Rectangle을 축소한다. 
**/

void GRect::DeflateRect(LDOUBLE DeltaLeft,LDOUBLE DeltaBottom, LDOUBLE DeltaRight, LDOUBLE DeltaTop )
{
	// 오른손 실 좌표계에서의 정의 
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
@@	주어진 GRect Object가 현재 GRect Object에 완전히 포함되면 TRUE를
@@	반환한다. 
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
@@	포함 테스트를 수행하기 위한 Flag를 설정한다. 
@@	Line Clipping알고리즘에서 
******************************************************************************
@@	B4 | B3 | B2 | B1 
@@	-----------------------------------------------------------------
@@	B1 : 선분의 끝점이 Window의 위쪽에 놓이는 경우         : true
@@	B2 : 선분의 끝점이 Window의 아래쪽에 놓이는 경우      : true
@@	B3 : 선분의 끝점이 Window의 오른쪽에 놓이는 경우      : true
@@	B4 : 선분의 끝점이 Window의 왼쪽에 놓이는 경우         : true 
@@  -----------------------------------------------------------------
@@	Bit Field로 구성함이 합당하나 현재는 Hex Field로 구성되어 있다.
*/
/*
void _SetIntersectFlag(GRect &Window,LINT &P1Flag,LINT &P2Flag,GPoint3D P1,GPoint3D P2)
{
	memset(&P1Flag,0,sizeof(LINT)); // P1Flag 초기화 
	memset(&P2Flag,0,sizeof(LINT)); // P2Flag 초기화  
	
	//	주어진 점 P1과 P2에 대한 Clipping Flag를 설정한다 
	if(P1.x < Window.left)    
		P1Flag |= 0x1000;				//	B1 : P1이 Window의 왼쪽에 놓이는 경우  
	if(P1.x > Window.right)      
		P1Flag |= 0x0100;				//	B2 : P1이 Window의 오른쪽에 놓이는 경우
	if(P1.y > Window.top)		
		P1Flag |= 0x0001;				//	B3 : P1이 Window의 위쪽에 놓이는 경우
	if(P1.y < Window.bottom)
		P1Flag |= 0x0010;				//	B4 : P1이 Window의 아래쪽에 놓이는 경우  

	
	if(P2.x < Window.left)
		P2Flag |= 0x1000;				//	B1 : P2가 Window의 왼쪽에 놓이는 경우  
	if(P2.x > Window.right)
		P2Flag |= 0x0100;				//	B2 : P2가 Window의 오른쪽에 놓이는 경우
	if(P2.y > Window.top)
		P2Flag |= 0x0001;				//	B3 : P2가 Window의 위쪽에 놓이는 경우
	if(P2.y < Window.bottom)
		P2Flag |= 0x0010;				//	B4 : P2가 Window의 아래쪽에 놓이는 경우  
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
@@	주어진 Rectangle을 길이가 긴 변이 1이되도록 하는 값을 반환한다. 
@@	즉 Normalizing Factor를 반환한다.
@@	Viewing 좌표계로의 변환후 Display될 Model Data의 Projection 작업에서 
@@	Model Data에 대한 Normalized Device Coordinate를 생성하기 위한 
@@	용도로 사용할 목적.
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
// NDC에서 Window Rect의 depth는 View좌표계 기준의 NDC변환된 View Plane의 
// Z 좌표를 의미한다. 

