#include "stdafx.h"
#include "GRenderAll.h"
//#include "GMathUtil.h"
/*
LDOUBLE RoundUpDouble(LDOUBLE dVal)
{
	if(dVal > 0)
		return floor(dVal+V_LGL_ROUND_UP_DOUBLE);
	else
	   return ceil(dVal-V_LGL_ROUND_UP_DOUBLE);
}

LINT	RoundUpInt(LDOUBLE dVal)
{
  if(dVal > 0)
	   return (LINT)floor(dVal+V_LGL_ROUND_UP_INT);
	else
	   return (LINT)ceil(dVal-V_LGL_ROUND_UP_INT);
}

LINT  RoundOffInt(LDOUBLE dVal)
{
	if(dVal > 0)
		return (LINT)floor(dVal);
	else
		return (LINT)ceil(dVal);
}

LBOOL IsSameDouble(LDOUBLE &a, LDOUBLE &b)
{
  if( fabs(a-b) <= V_LGLDBL_EPSILON) return TRUE;
  return FALSE;
}

int CompareInt( const void* arg1, const void* arg2 )
{
  return *(int*)arg1 - *(int*)arg2;
}

int CompareDoubleV(double arg1, double arg2)
{
  if( fabs(arg1-arg2) <= V_LGLDBL_EPSILON) return 0;
  if(arg1 < arg2) return -1;
  return 1;
}

int CompareDouble (const void* arg1, const void* arg2)
{
  if(IsSameDouble(*(double*)arg1,*(double*)arg2)) return 0;
  if(*(double*)arg1 < *(double*)arg2) return -1;
  return 1;
}

int CompareStr(const void* arg1, const void* arg2)
{
   //Compare all of both strings: 
   return _stricmp( (char* ) arg1,  ( char* ) arg2 );
}
*/  
///////////////////////////////////////////////////////
// 주어진 Tolerance로 주어진 LDOUBLE 값을 변환한다. 

char* _strPrecF[]={
  "%.0f",
  "%.1f",
  "%.2f",
  "%.3f",
  "%.4f",
  "%.5f",
  "%.6f",
  "%.7f",
  "%.8f",
  "%.9f",
  "%.10f",
  "%.11f",
  "%.12f",
  "%.13f",
  "%.14f",
  "%.15f"
};

LDOUBLE GMU::RoundUpDouble(LDOUBLE dVal, LINT nPrec)
{

 LDOUBLE rdVal;
 CString strF;
 strF.Format(_strPrecF[nPrec],dVal);
 sscanf_s(strF.GetBuffer(strF.GetLength()+1),_strPrecF[nPrec],&rdVal); 
 strF.ReleaseBuffer();

 return rdVal;
}

LDOUBLE GMU::RoundUpDoubleR(LDOUBLE dVal, LINT nPrec)
{
  double TVal = pow(10.0, (double)nPrec); 
  dVal = dVal* TVal;

  if(dVal > 0.)
		dVal =  floor(dVal+V_LGL_ROUND_UP_DOUBLE);
	else
	  dVal = ceil(dVal-V_LGL_ROUND_UP_DOUBLE);
  
  dVal = dVal / TVal;

  return dVal;

}


LDOUBLE GMU::RoundUpDouble(LDOUBLE dVal)
{
  //사용되지 않음.
	if(dVal > 0.)
		return floor(dVal+V_LGL_ROUND_UP_DOUBLE);
	else
	  return ceil(dVal-V_LGL_ROUND_UP_DOUBLE);
}

LINT	GMU::RoundUpInt(LDOUBLE dVal)  
{
  return RoundOffInt(dVal);

  //if(dVal > 0.)
	//   return (LINT)floor(dVal+V_LGL_ROUND_UP_INT);
	//else
	//   return (LINT)ceil(dVal-V_LGL_ROUND_UP_INT);
  //------------------------------------------------
  // 속도개선 : PIG
  //if(dVal > 0.)
	//   return (LINT)(dVal+V_LGL_ROUND_UP_INT);
	//else
	//   return (LINT)(dVal-V_LGL_ROUND_UP_INT);
  //<<-- GL GDI 혼용 에서 문제 발생 
  // 사용하지 말것. 
}

LINT  GMU::RoundOffInt(LDOUBLE dVal)  
{
	
  if(dVal > 0.)
		return (LINT)floor(dVal);
	else
		return (LINT)ceil(dVal);
  //---------------------------------------------------
  // 속도개선 : PIG
  //return (LINT)dVal; <<-- GL GDI 혼용 에서 문제 발생 
  //                   <<-- 사용하지 말것. 
}

LBOOL GMU::IsZero(LDOUBLE Val)
{
  if( Val <= V_LGLDBL_EPSILON && Val >= -V_LGLDBL_EPSILON) return TRUE;
  return FALSE;
}

void GMU::MakeNotZeroMinus(LDOUBLE &V1,LDOUBLE &V2,LDOUBLE &V3)
{
  if(IsZero(V1))
    V1 -= V_LGLDBL_EPSILON;
  if(IsZero(V2))
    V2 -= V_LGLDBL_EPSILON;
  if(IsZero(V3))
    V3 -= V_LGLDBL_EPSILON;
}

void GMU::MakeNotZeroPlus(LDOUBLE &V1,LDOUBLE &V2,LDOUBLE &V3)
{
  if(IsZero(V1))
    V1 += V_LGLDBL_EPSILON;
  if(IsZero(V2))
    V2 += V_LGLDBL_EPSILON;
  if(IsZero(V3))
    V3 += V_LGLDBL_EPSILON;
}
//LBOOL GMU::IsZeroFor3DLineIntersect

LBOOL GMU::IsZeroForDisplay(LDOUBLE Val)
{
  if( Val <= V_LGLDBL_EPSILON_FOR_DISPLAY && Val >= -V_LGLDBL_EPSILON_FOR_DISPLAY ) return TRUE;
  return FALSE;
}

LBOOL GMU::IsSameDouble(LDOUBLE &a, LDOUBLE &b)
{
  if(IsZero(a) && IsZero(b)) return TRUE;
  if(a * b >= 0.0)
  {
    if(fabs(a-b) <= V_LGLDBL_EPSILON) return TRUE;
  }
  return FALSE;
}

LBOOL GMU::IsSameDoubleV(LDOUBLE a, LDOUBLE b)
{
  if(IsZero(a) && IsZero(b)) return TRUE;
  if(a * b >= 0.0)
  {
    if(fabs(a-b) <= V_LGLDBL_EPSILON) return TRUE;
  }
  return FALSE;
}

int  GMU::CompareInt( const void* arg1, const void* arg2 )
{
  return *(int*)arg1 - *(int*)arg2;
}

int  GMU::CompareDoubleV(double arg1, double arg2)
{
  //if( fabs(arg1-arg2) <= V_LGLDBL_EPSILON) return 0;
  if(IsSameDoubleV(arg1,arg2)) return 0;
  if(arg1 < arg2) return -1;
  return 1;
}

/*
int __CompareDoubleV(double& a, double & b)
{
  if(IsSameDoubleV(arg1,arg2)) return 0;
  if(arg1 < arg2) return -1;
  return 1;
}
// == 
LBOOL GMU::EQ(double& arg1,double& arg2)
{
  if(arg1 * arg2 >= 0.0)
  {
    if(fabs(arg1-arg2) <= V_LGLDBL_EPSILON) return TRUE;
  }
  return FALSE;
}
// >= 
LBOOL GMU::GE(double& arg1,double& arg2)
{
  

}
// >
LBOOL GMU::GT(double& arg1,double& arg2)
{

}
// <= 
LBOOL GMU::LE(double& arg1,double& arg2)
{

}
// <
LBOOL GMU::LT(double& arg1,double& arg2)
{

}
*/

int  GMU::CompareDouble (const void* arg1, const void* arg2)
{
  if(IsSameDouble(*(double*)arg1,*(double*)arg2)) return 0;
  if(*(double*)arg1 < *(double*)arg2) return -1;
  return 1;
}

int GMU::CompareStr(const void* arg1, const void* arg2)
{
   //Compare all of both strings: 
   return _stricmp( (char* ) arg1,  ( char* ) arg2 );
}

