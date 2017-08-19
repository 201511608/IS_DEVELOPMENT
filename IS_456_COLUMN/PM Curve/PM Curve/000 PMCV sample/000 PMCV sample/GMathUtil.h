#ifndef __GMATHUTIL_H__
#define __GMATHUTIL_H__
#include "LGLType.h"

class GMU
{
public:
	static LBOOL   IsZeroForDisplay(LDOUBLE Val);
  //static LBOOL   IsZero(LDOUBLE Val);
  static LBOOL   IsZero(LDOUBLE Val);
  static void    MakeNotZeroMinus(LDOUBLE &V1,LDOUBLE &V2,LDOUBLE &V3);
  static void    MakeNotZeroPlus(LDOUBLE &V1,LDOUBLE &V2,LDOUBLE &V3);
  static LDOUBLE RoundUpDouble(LDOUBLE dVal, LINT nPrec);
  static LDOUBLE RoundUpDouble(LDOUBLE dVal);
  static LINT	   RoundUpInt(LDOUBLE dVal);
  static LINT    RoundOffInt(LDOUBLE dVal);
  static LBOOL   IsSameDouble(LDOUBLE &a, LDOUBLE &b);
  static LBOOL   IsSameDoubleV(LDOUBLE a, LDOUBLE b);
  static LDOUBLE RoundUpDoubleR(LDOUBLE dVal, LINT nPrec);

  static int CompareInt       (const void* arg1, const void* arg2);
  static int CompareDouble (const void* arg1, const void* arg2);
  static int CompareStr       (const void* arg1, const void* arg2);

  static int CompareDoubleV (double arg1, double arg2);
};

/*
extern LDOUBLE RoundUpDouble(LDOUBLE dVal);
extern LINT	   RoundUpInt(LDOUBLE dVal);
extern LINT    RoundOffInt(LDOUBLE dVal);
extern LBOOL   IsSameDouble(LDOUBLE &a, LDOUBLE &b);

extern int CompareInt       (const void* arg1, const void* arg2);
extern int CompareDouble (const void* arg1, const void* arg2);
extern int CompareStr       (const void* arg1, const void* arg2);

extern int CompareDoubleV (double arg1, double arg2);
*/

#endif