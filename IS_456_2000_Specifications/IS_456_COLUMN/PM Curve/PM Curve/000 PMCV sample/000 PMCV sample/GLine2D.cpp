#include "stdafx.h"
#define __INC_GLINE2D_H__
#include "GRenderAll.h"

#ifdef  _DEBUG
#define new DEBUG_NEW
#undef  THIS_FILE
static  char THIS_FILE[] = __FILE__;
#endif
/************************************************
   *   GLine2D Class Implementation
   *
 * * *
  ***
   */
GLine2D::GLine2D()
{

}

GLine2D::~GLine2D()
{

}

GLine2D GLine2D::operator = (GLine2D Source)
{
  P1 = Source.P1;
  P2 = Source.P2;
  /*
  memcpy(&P1.x,&Source.P1.x,sizeof(LINT));
	memcpy(&P1.y,&Source.P1.y,sizeof(LINT));
	memcpy(&P2.x,&Source.P2.x,sizeof(LINT));
	memcpy(&P2.y,&Source.P2.y,sizeof(LINT));
	*/
	return *this;
}

void GLine2D::Set(GPoint2D p1,GPoint2D p2)
{
	P1 = p1;
	P2 = p2;
}

void GLine2D::Set(GLine3D& Line3D)
{
	P1.Set(Line3D.P1);
	P2.Set(Line3D.P2);
}

CString GLine2D::GetString()
{
	CString StrTemp;
	StrTemp.Format("(%d, %d),(%d, %d)",P1.x,P1.y,P2.x,P2.y);
	
	return StrTemp;
}

void GLine2D::SetZero()
{
	P1.Set(0,0);
	P2.Set(0,0);
}

/***********************************************************
@@	현재 Line을 구성하는 점들이 유효한 Line을 구성하는 조건
@@	인지를 검사한다. 즉 Line을 구성하는 두점이 같은 점이면 
@@	Line은 유효하지 않은 Line으로 판정한다. 
**/
LBOOL GLine2D::IsValid()
{
	if(P1 == P2 ) return FALSE;
	return TRUE;
}

/***********************************************************
@@	Line이 현재 좌표계의 X축과 이루는 각도를 반환한다.
@@  반환값은 Radian값이다. 
@@	X축에 수직일 경우에는 Y좌표를 비교해서 각도의 부호를 결정
@@	한다. 
**/
LDOUBLE GLine2D::GetAngle()
{
	LDOUBLE ValTan;
	LDOUBLE DivVal;
	DivVal = P2.x - P1.x;
	if(fabs(DivVal) <= V_LGLDBL_EPSILON)
	{
		/***********************************************
		@@ 원점에 가까운 점을 시작점으로 한다. 
		**/
		if(P2.x*P2.x + P2.y*P2.y > P1.x*P1.x + P1.y*P1.y)
		{
			return RAD(90.0);
		}
		else
		{
			return RAD(-90.0);
		}
	}
	ValTan = (P2.y - P1.y) / (P2.x - P1.x);
	return atan(ValTan);
}

LINT GLine2D::Length()
{
	LDOUBLE LengthD,DeltaX,DeltaY;
	DeltaX = (P2.x - P1.x);
	DeltaY = (P2.y - P1.y);

	LengthD = sqrt( DeltaX*DeltaX + DeltaY*DeltaY );

  return GMU::RoundUpInt(LengthD);

}

void GLine2D::Draw(CDC *pdc)
{
	pdc->MoveTo(P1.x,P1.y);
	pdc->LineTo(P2.x,P2.y);
}


void GLine2D::DrawStyled(CDC *pdc,COLORREF color,int nStyle)
{
	int idDC = pdc->SaveDC();
	//CPen *pOldPen;
	CPen pen;
	pen.CreatePen(nStyle,1,color);
	//pOldPen = pdc->SelectObject(&pen);
	pdc->SelectObject(&pen);
	Draw(pdc);
	//pdc->SelectObject(pOldPen);
	pdc->RestoreDC(idDC);
	pen.DeleteObject();
}


void GLine2D::Draw(CDC *pdc,COLORREF color)
{
	int idDC = pdc->SaveDC();
	//CPen *pOldPen;
	CPen pen;
	pen.CreatePen(PS_SOLID,1,color);
	//pOldPen = pdc->SelectObject(&pen);
	pdc->SelectObject(&pen);
	Draw(pdc);
	//pdc->SelectObject(pOldPen);
	pdc->RestoreDC(idDC);
	pen.DeleteObject();
}

void GLine2D::Draw(CDC* pdc, COLORREF color, int LineWidth)
{
	int idDC = pdc->SaveDC();
	//CPen *pOldPen;
	CPen pen;
	pen.CreatePen(PS_SOLID,LineWidth,color);
	//pOldPen = pdc->SelectObject(&pen);
	pdc->SelectObject(&pen);
	Draw(pdc);
	//pdc->SelectObject(pOldPen);
	pdc->RestoreDC(idDC);
	pen.DeleteObject();
}

LINT GLine2D::GetMinX(LINT CurMinX)
{
	LINT MinX;
	MinX = CurMinX;

	MinX = min(P1.x,MinX);
	MinX = min(P2.x,MinX);

	return MinX;
}

LINT GLine2D::GetMinY(LINT CurMinY)
{
	LINT MinY;
	MinY = CurMinY;

	MinY = min(P1.y,MinY);
	MinY = min(P2.y,MinY);

	return MinY;
}


LINT GLine2D::GetMaxX(LINT CurMaxX)
{
	LINT MaxX;
	MaxX = CurMaxX;

	MaxX = max(P1.x,MaxX);
	MaxX = max(P2.x,MaxX);

	return MaxX;
}

LINT GLine2D::GetMaxY(LINT CurMaxY)
{
	LINT MaxY;
	MaxY = CurMaxY;

	MaxY = max(P1.y,MaxY);
	MaxY = max(P2.y,MaxY);

	return MaxY;
}

int GLine2D::operator == (GLine2D line)
{
	if(P1 != line.P1 ) return 0;
	if(P2 != line.P2 ) return 0;

	return 1;
}

int GLine2D::operator != (GLine2D line)
{
	if(P1 != line.P1 ) return 1;
	if(P2 != line.P2 ) return 1;

	return 0;
}
