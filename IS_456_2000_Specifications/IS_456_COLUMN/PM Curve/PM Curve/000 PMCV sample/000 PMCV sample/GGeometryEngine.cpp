#include "stdafx.h"
#define __INC_GRENDERENGINE_H__
#include "GRenderAll.h"
#include <math.h>


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

GGeometryEngine::GGeometryEngine()
{
}

GGeometryEngine::~GGeometryEngine()
{
}

LBOOL GGeometryEngine::_IsWithin2D(GRect rect, GPoint3D * pPoint)
{
	if(rect.left <= pPoint->x && rect.bottom <= pPoint->y &&
	   rect.top >= pPoint->y && rect.right  >= pPoint->x) return TRUE;
	return FALSE;	
}



BOOL GGeometryEngine::NewellsMethod(CArray<GPoint3D,GPoint3D&>&VArr,GVector&NorVector)
{
  LDOUBLE  nor_X=0,nor_Y=0,nor_Z=0;
  
  int nV = VArr.GetSize(), j;
  for (int i=0;i<nV;i++) 
  {
    if (i== nV-1 ) j=0;
    else j=i+1;
    nor_X=nor_X+(VArr[i].y-VArr[j].y)*(VArr[i].z+VArr[j].z);
    nor_Y=nor_Y+(VArr[i].z-VArr[j].z)*(VArr[i].x+VArr[j].x);
    nor_Z=nor_Z+(VArr[i].x-VArr[j].x)*(VArr[i].y+VArr[j].y);
  }

  NorVector.Set(nor_X,nor_Y,nor_Z);
  NorVector.MakeUnit();
  return TRUE;
}

BOOL GGeometryEngine::LeftInnerIntersectPoint(GPoint3D P1,GPoint3D PC, GPoint3D P2,GVector NorVector, double Dist,GPoint3D& IPoint)
{
  // P1 -> PC -> P2
  GVector P1PC,PCP2  ; // XDir 로 간주한다. 
  GVector P1PCY,PCP2Y;

  P1PC.Set(P1,PC); P1PC.MakeUnit();
  PCP2.Set(PC,P2); PCP2.MakeUnit();
  
  // J -> K -> I 
  P1PCY = NorVector * P1PC;
  PCP2Y = NorVector * PCP2;

  P1PCY.MakeUnit(); 
  PCP2Y.MakeUnit(); 

  P1PCY.ScalingMySelf(Dist);
  PCP2Y.ScalingMySelf(Dist);

  GLine3D Line1,Line2;
  Line1.Set(P1,PC);
  Line2.Set(PC,P2);
  Line1.Offset(P1PCY.x,P1PCY.y,P1PCY.z);
  Line2.Offset(PCP2Y.x,PCP2Y.y,PCP2Y.z);

  GPoint3D IPos;
  
  if(GetIntersectLine3DPoint(Line1,Line2,IPos))
  {
    IPoint = IPos;
  }
  else
  {
    BOOL bZero1C = GMU::IsZero(P1PC.Length());
    BOOL bZeroC2 = GMU::IsZero(PCP2.Length());
    if(bZero1C && bZeroC2) 
      return FALSE;

    if(bZero1C)
    {
      IPoint = Line2.P1;
      return TRUE;
    }
    
    if(bZeroC2) 
    {
      IPoint = Line1.P2;
      return TRUE;
    }
    
    if(P1PC == PCP2) 
    {
      IPoint = Line1.P2;
      return TRUE;
    }
    
    PCP2.Inverse();
    if(P1PC == PCP2) 
    {
      IPos = Line1.P1;
      return FALSE;
    }

    IPoint = Line1.P2;
  }
  return TRUE;
}

// bIsFactor (TRUE) SFactorOrDist --> Shrink Factor (FALSE) SFactorOrDist --> Distance
BOOL GGeometryEngine::PolygonShrink_loveme(CArray<GPoint3D,GPoint3D&>&VArr,CArray<GPoint3D,GPoint3D&>&ShrinkedVArr,
                            double SFactorOrDist,BOOL bIsFactor) 
{
  ShrinkedVArr.RemoveAll();
  LDOUBLE TotalLength = 0;
  GVector NorVector,TShrVector;

  GPoint3D P1,PC,P2,TPoint;
  GLine3D  TLength;

  int nV;

  nV = VArr.GetSize();
  
  if(nV < 3) return FALSE;

  CArray<GVector,GVector&> ShrVectArr;
  
  NewellsMethod(VArr,NorVector);
  
  double OffsetDist;

  if(bIsFactor)
  {
    GLine3D LengthL;
    double TotalLength = 0;
    for(int i = 0; i < nV ; i++)
    {
      LengthL.P1 = VArr[i];  
      if(i == nV-1)
        LengthL.P2 = VArr[0];
      else
        LengthL.P2 = VArr[i+1];
        
      TotalLength += LengthL.Length();
    }
    OffsetDist = TotalLength*SFactorOrDist;
  }
  else
  {
    OffsetDist = SFactorOrDist;
  }
   
  
  GPoint3D IPos;
  GLine3D Line1,Line2;
  for(int i = 0 ; i < nV ; i++)
  {
    if(i == 0)
    {
      P1 = VArr[nV-1]; PC = VArr[0]; P2 = VArr[1];
    }
    else if(i == nV-1)
    {
      P1 = VArr[nV-2]; PC = VArr[nV-1]; P2 = VArr[0];
    }
    else
    {
      P1 = VArr[i-1] ; PC = VArr[i]   ; P2 = VArr[i+1];
    }
      
    if(LeftInnerIntersectPoint(P1,PC,P2,NorVector,OffsetDist,IPos))
      ShrinkedVArr.Add(IPos);
  }
  return TRUE;
}

// 무한직선으로 간주 -------------------------------------------------------------------
BOOL GGeometryEngine::GetIntersectLine3DPoint(GLine3D & L1,GLine3D & L2,GPoint3D & IPos)
{
  GVector V1,V2;
	
  //GPoint3D L1P1,L1P2,L2P1,L2P2;
  
  V1.Set(L1.P1,L1.P2);
  V2.Set(L2.P1,L2.P2);
  
  //V1.MakeUnit();
  //V2.MakeUnit();
  
  if(!V1.IsValid() || !V2.IsValid()) return FALSE;
	
  GVector TVec = V1.Cross(V2);
  TVec.MakeUnit();
  if(!TVec.IsValid()) return FALSE;;
  
  LDOUBLE D = sqrt((V1.y*V2.z - V2.y*V1.z)*(V1.y*V2.z - V2.y*V1.z) +
                   (V1.z*V2.x - V2.z*V1.x)*(V1.z*V2.x - V2.z*V1.x) +
                   (V1.x*V2.y - V2.x*V1.y)*(V1.x*V2.y - V2.x*V1.y));
  if(GMU::IsZero(D)) return FALSE;                   

  LDOUBLE Dist = fabs(((L2.P1.x - L1.P1.x)*(V1.y*V2.z-V2.y*V1.z)-
                       (L2.P1.y - L1.P1.y)*(V1.x*V2.z-V2.x*V1.z)+
                       (L2.P1.z - L1.P1.z)*(V1.x*V2.y-V2.x*V1.y))/D);
  if(!GMU::IsZero(Dist)) return FALSE;
	  
	LDOUBLE  Q1 = V2.x * V1.y - V1.x * V2.y;
	LDOUBLE  Q2 = V2.x * V1.z - V1.x * V2.z;
	LDOUBLE  Q3 = V2.y * V1.z - V1.y * V2.z;
	LDOUBLE  t;

  
  if(!GMU::IsZero(Q1))
	{
		t = (V1.x*(L2.P1.y - L1.P1.y) - V1.y *(L2.P1.x - L1.P1.x)) / Q1;
	}
  else if(!GMU::IsZero(Q2))
	{
    t = (V1.x*(L2.P1.z - L1.P1.z) - V1.z *(L2.P1.x - L1.P1.x)) / Q2;
	}
  else if(!GMU::IsZero(Q3))
	{
    t = (V1.y*(L2.P1.z - L1.P1.z) - V1.z *(L2.P1.y - L1.P1.y)) / Q3;
	}
	else
		return FALSE;
  
	IPos.x = L2.P1.x + V2.x * t;
	IPos.y = L2.P1.y + V2.y * t;
	IPos.z = L2.P1.z + V2.z * t;

	return TRUE;
}