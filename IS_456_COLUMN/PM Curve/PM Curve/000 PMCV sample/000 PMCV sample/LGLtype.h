#ifndef __LGLTYPE_H__
#define __LGLTYPE_H__

/******  DEFINE IN BORLANDC++ 3.0 MATH.H  *******/
#define M_E         2.71828182845904523536
#define M_LOG2E     1.44269504088896340736
#define M_LOG10E    0.434294481903251827651
#define M_LN2       0.693147180559945309417
#define M_LN10      2.30258509299404568402
#define M_PI        3.14159265358979323846
#define M_PI_2      1.57079632679489661923
#define M_PI_4      0.785398163397448309616
#define M_1_PI      0.318309886183790671538
#define M_2_PI      0.636619772367581343076
#define M_1_SQRTPI  0.564189583547756286948
#define M_2_SQRTPI  1.12837916709551257390
#define M_SQRT2     1.41421356237309504880
#define M_SQRT_2    0.707106781186547524401

/****** Generate by L.C.G ***********************************/
#define M_DEG_TO_RAD  0.0174532925199432950   // ( M_PI  / 180.0 )
#define M_RAD_TO_DEG  57.2957795130823230	  // ( 180.0 / M_PI	 )

/****** Define by L.C.G **************************/  
#define RAD(Degree) ((Degree) * M_DEG_TO_RAD)
#define DEG(Rad)    ((Rad   ) * M_RAD_TO_DEG) 

/**************************************************
@@	From VC++ FLOAT.H
@@	Smallest value such that 1.0+LDBL_EPSILON != 1.0  
@@  Same as define macro DBL_EPSILON in FLOAT.H
**/
//#define LGLDBL_EPSILON              DBL_EPSILON * 100
//#define LGLDBL_EPSILON              1.0E-9
#define LGLDBL_EPSILON_FOR_DISPLAY    LGLDBL_EPSILON * 100
//#define LGLDBL_MBR_PRECISION        1.0E-9
#define LGLDBL_EPSILON                2.2204460492503131e-016
#define LGLDBL_MBR_PRECISION          2.2204460492503131e-016

#define LGL_LINE_INTERSECT_EPSILON    LGLDBL_EPSILON 
#define LGL_ROUND_UP_INT              0.5 // LDOUBLE에서 LINT로의 형변환을 위해서 더해지는 값
#define LGL_ROUND_UP_DOUBLE           0.5 // LDOUBLE값의 반올림을 위해 더해지는 값
#define LGLDBL_AREA_MIN              -DBL_MAX / 10.0
#define LGLDBL_AREA_MAX               DBL_MAX / 10.0


#define V_LGLDBL_EPSILON              1.0E-9
#define V_LGLDBL_EPSILON_FOR_DISPLAY  LGLDBL_EPSILON * 100
#define V_LGLDBL_MBR_PRECISION        1.0E-9
//#define LGL_LINE_INTERSECT_EPSILON   10E-9 // Line Intersect 조건 검사때의 정도 
#define V_LGL_LINE_INTERSECT_EPSILON  LGLDBL_EPSILON 
#define V_LGL_ROUND_UP_INT            0.5 // LDOUBLE에서 LINT로의 형변환을 위해서 더해지는 값
#define V_LGL_ROUND_UP_DOUBLE         0.5 // LDOUBLE값의 반올림을 위해 더해지는 값
#define V_LGLDBL_AREA_MIN             -DBL_MAX / 10.0
#define V_LGLDBL_AREA_MAX             DBL_MAX / 10.0
/*
extern double  V_LGLDBL_EPSILON              ;//1.0E-9
extern double  V_LGLDBL_EPSILON_FOR_DISPLAY  ;//LGLDBL_EPSILON * 100
extern double  V_LGLDBL_MBR_PRECISION        ;//1.0E-9
//#define LGL_LINE_INTERSECT_EPSILON   10E-9 // Line Intersect 조건 검사때의 정도 
extern double  V_LGL_LINE_INTERSECT_EPSILON  ;//LGLDBL_EPSILON 
extern double  V_LGL_ROUND_UP_INT            ;//0.5 // LDOUBLE에서 LINT로의 형변환을 위해서 더해지는 값
extern double  V_LGL_ROUND_UP_DOUBLE         ;//0.5 // LDOUBLE값의 반올림을 위해 더해지는 값
extern double  V_LGLDBL_AREA_MIN             ;//-DBL_MAX / 10.0
extern double  V_LGLDBL_AREA_MAX             ;//DBL_MAX / 10.0
*/

#define LFLOAT                float
#define LDOUBLE               double
#define LLONG_DOUBLE          long double 

#define LSHORT                short int
#define LINT                  int
#define LUINT                 unsigned int 
#define LLONG                 long int
#define LUNSIGNED_LONG        unsigned long int 

#define LBYTE                 BYTE
#define LCHAR                 char 
#define LUNSIGNED_CHAR        unsigned char
#define LBOOL                 BOOL
#define LPOSITION             POSITION

#define LLOGPEN               LOGPEN
#define LLOGBRUSH             LOGBRUSH
#define LLOGFONT              LOGFONT
#ifndef REM
	#define REM(param)
#endif
/*
double  V_LGLDBL_EPSILON              = 1.0E-9;
double  V_LGLDBL_EPSILON_FOR_DISPLAY  = LGLDBL_EPSILON * 100;
double  V_LGLDBL_MBR_PRECISION        = 1.0E-6;
//#define LGL_LINE_INTERSECT_EPSILON   10E-9 // Line Intersect 조건 검사때의 정도 
double  V_LGL_LINE_INTERSECT_EPSILON  = LGLDBL_EPSILON ;
double  V_LGL_ROUND_UP_INT            = 0.5 ;// LDOUBLE에서 LINT로의 형변환을 위해서 더해지는 값
double  V_LGL_ROUND_UP_DOUBLE         = 0.5 ;// LDOUBLE값의 반올림을 위해 더해지는 값
double  V_LGLDBL_AREA_MIN             = -DBL_MAX / 10.0;
double  V_LGLDBL_AREA_MAX             =  DBL_MAX / 10.0;
*/

#endif
