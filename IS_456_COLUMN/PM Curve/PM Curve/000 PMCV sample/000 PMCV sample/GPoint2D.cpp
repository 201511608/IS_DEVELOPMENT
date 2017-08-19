#include "stdafx.h"
#define __INC_GPOINT2D_H__
#include "GRenderAll.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
/************************************************
   *   GPoint2D Class Implementation
   *
 * * *
  ***
   */
GPoint2D::GPoint2D()
{
	x = 0;
  y = 0;
}

void GPoint2D::Init()
{
	x = 0;
  y = 0;
}

GPoint2D::GPoint2D(LINT X, LINT Y)
{
	x = X;
  y = Y;
}

// Copy Constructor
GPoint2D::GPoint2D(GPoint2D& S)
{
	x = S.x;
  y = S.y;
}

GPoint2D::~GPoint2D()
{

}

GPoint2D& GPoint2D::operator = (const GPoint2D& Source )
{
	x = Source.x;
  y = Source.y;
	return (*this);
}

int GPoint2D::operator == (const GPoint2D point )
{
	if(x != point.x) return 0;
	if(y != point.y) return 0;
	return 1;
}


int GPoint2D::operator != (const GPoint2D point)
{
	if(x != point.x) return 1;
	if(y != point.y) return 1;
	return 0;
}

GPoint2D GPoint2D::operator + (const GPoint2D point)
{
	GPoint2D Temp;
	Temp.x = x + point.x;
	Temp.y = y + point.y;
	return Temp;
}

GPoint2D GPoint2D::operator - (const GPoint2D point)
{
	GPoint2D Temp;
	Temp.x = x - point.x;
	Temp.y = y - point.y;
	return Temp;
}
/*
LINT RoundUpInt(LDOUBLE dVal)
{
	if(dVal > 0)
	   return (LINT)floor(dVal+V_LGL_ROUND_UP_INT);
	else
	   return (LINT)ceil(dVal-V_LGL_ROUND_UP_INT);
}

LINT RoundOffInt(LDOUBLE dVal)
{	
	if(dVal > 0)
		return (LINT)floor(dVal);
	else
		return (LINT)ceil(dVal);
}
*/
/**************************************************
@@	Friend Member Function
*/
void GPoint2D::Set(GPoint3D& point3d)
{
  x = GMU::RoundUpInt(point3d.x);
  y = GMU::RoundUpInt(point3d.y);   
}


GPoint2D operator / (GPoint2D Src,LINT DivideNum)
{
	GPoint2D Temp;
	Temp.x = Src.x / DivideNum;
	Temp.y = Src.y / DivideNum;

	return Temp; 
}


GPoint2D operator * (GPoint2D Src, LINT MultNum)
{
	GPoint2D Temp;
	Temp.x = Src.x * MultNum;
	Temp.y = Src.y * MultNum;
	
	return Temp;
}


GPoint2D operator * (LINT MultNum,GPoint2D Src)
{
	GPoint2D Temp;
	Temp.x = Src.x * MultNum;
	Temp.y = Src.y * MultNum;
	
	return Temp;
}

void GPoint2D::Set(LINT X, LINT Y)
{
  x = X;
  y = Y;
}


CString GPoint2D::GetString(CString strFormat)
{
	CString strReturn;
    strReturn.Format(strFormat,x,y);
	return strReturn;
}

void GPoint2D::Draw(CDC *pdc,COLORREF Color)
{
	pdc->SetPixelV(x,y,Color);
}

void GPoint2D::DrawLarge(CDC *pdc,COLORREF Color)
{
  /*
  BOOL bSucc = TRUE;
  
  bSucc &= pdc->SetPixelV(x,y+1,Color);
  bSucc &= pdc->SetPixelV(x,y-1,Color);
  bSucc &= pdc->SetPixelV(x-1,y,Color);
  bSucc &= pdc->SetPixelV(x+1,y,Color);
  bSucc &= pdc->SetPixelV(x+1,y+1,Color);
  bSucc &= pdc->SetPixelV(x+1,y-1,Color);
  bSucc &= pdc->SetPixelV(x-1,y+1,Color);
  bSucc &= pdc->SetPixelV(x-1,y-1,Color);

  if(bSucc == FALSE) GSaveHistoryNF("SetPixelV :: Fail!");
  */
  /*
  BOOL bSucc = TRUE;

  if( Color != pdc->SetPixel(x,y+1,Color  )) bSucc = FALSE;
  if( Color != pdc->SetPixel(x,y-1,Color  )) bSucc = FALSE;
  if( Color != pdc->SetPixel(x-1,y,Color  )) bSucc = FALSE;
  if( Color != pdc->SetPixel(x+1,y,Color  )) bSucc = FALSE;
  if( Color != pdc->SetPixel(x+1,y+1,Color)) bSucc = FALSE;
  if( Color != pdc->SetPixel(x+1,y-1,Color)) bSucc = FALSE;
  if( Color != pdc->SetPixel(x-1,y+1,Color)) bSucc = FALSE;
  if( Color != pdc->SetPixel(x-1,y-1,Color)) bSucc = FALSE;

  if(bSucc == FALSE) GSaveHistoryNF("SetPixel :: Fail!");
  */
  
  /*
  int nSaveDC = pdc->SaveDC();
  
  LOGBRUSH lb;
	lb.lbStyle = BS_HOLLOW;
	CBrush br;
	br.CreateBrushIndirect(&lb);
	

  CRect rect;
  CPen pen(PS_SOLID,1,Color);
  pdc->SelectObject(&br);
  pdc->SelectObject(&pen);
  rect.SetRect(x-2,y+2,x+2,y-2);
  pdc->Rectangle(&rect);
  
  pdc->RestoreDC(nSaveDC);
  */
  CRect rect;
  rect.SetRect(x-2,y+2,x+2,y-2);
  pdc->Rectangle(rect);
  //pdc->Ellipse(rect);
  
}

void GPoint2D::DrawEllipse(CDC *pDC, int Width, int Height)
{
  /*
  LOGPEN LogPen;
  pDC->GetCurrentPen()->GetLogPen(&LogPen);
  LogPen.lopnColor;
    
  pDC->SetPixelV(x  ,y,LogPen.lopnColor);
  pDC->SetPixelV(x-1,y,LogPen.lopnColor);
  pDC->SetPixelV(x+1,y,LogPen.lopnColor);
  pDC->SetPixelV(x,y-1,LogPen.lopnColor);
  pDC->SetPixelV(x,y+1,LogPen.lopnColor);

  // Grad Color
  COLORREF PCo = LogPen.lopnColor;
  COLORREF GradColor;

  DWORD Red   = GetRValue(PCo) + 80;
	DWORD Blue  = GetBValue(PCo) + 80;
	DWORD Green = GetGValue(PCo) + 80;
	
  if(Red>255)   Red  = 255;
	if(Green>255) Green= 255;
	if(Blue>255)  Blue = 255;

  pDC->SetPixelV(x-1,y-1,GradColor);
  pDC->SetPixelV(x-1,y+1,GradColor);
  pDC->SetPixelV(x+1,y+1,GradColor);
  pDC->SetPixelV(x+1,y-1,GradColor);
  */



  
  
  
  
  
  //int DeltaX = GMU::RoundUpInt(Width/2.0);
  //int DeltaY = GMU::RoundUpInt(Height/2.0);
  int DeltaX = Width  / 2;
  int DeltaY = Height / 2;

  CRect ERect;
  ERect.SetRect(x-DeltaX,y+DeltaY,x+DeltaX,y-DeltaY);
  pDC->Ellipse(ERect);
  
  //pDC->Rectangle(&ERect);
}


void GPoint2D::DrawCross(CDC* pDC, int Length)
{
  int DeltaX = Length / 2;
  int DeltaY = Length / 2;

  pDC->MoveTo(x-DeltaX, y-DeltaY);
  pDC->LineTo(x+DeltaX, y+DeltaY);
  pDC->MoveTo(x-DeltaX, y+DeltaY);
  pDC->LineTo(x+DeltaX, y-DeltaY);
}
