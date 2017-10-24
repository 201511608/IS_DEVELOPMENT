#ifndef __GVECTOR_H__
#define __GVECTOR_H__

class GVector 
{
public:
	void ScalingMySelf(LDOUBLE r);
	void Inverse();
	LBOOL IsValid();
	LBOOL MakeNormalVector(GPoint3D P1, GPoint3D P2, GPoint3D P3, GPoint3D DirectionP);
	LBOOL MakeNormalVector(GPoint3D P1, GPoint3D P2, GPoint3D P3, GVector DirectionV);
	LBOOL MakeNormalVector(GPoint3D& P1,GPoint3D& P2,GPoint3D& P3);
	LDOUBLE x,y,z;
      
	GVector();     //dummy
  GVector(LDOUBLE X,LDOUBLE Y,LDOUBLE Z);
  GVector(GPoint3D StartPoint,GPoint3D EndPoint);
  GVector(GPoint3D &VectorPoint);
	
	LBOOL  IsAcuteAngle(GVector& OVec);
  LBOOL  IsObtuseAngle(GVector& OVec);
	int operator == (const GVector point );
	int operator != (const GVector point);
	GVector operator+(GVector);
  GVector operator-(GVector);
  GVector operator*(GVector);   //cross
	GVector Cross(GVector vec);
	LDOUBLE Dot(GVector);
  LDOUBLE CCWAngle2D();
  BOOL    Angle3D(GVector AVector, double & Angle);
	BOOL    Cosine3D(GVector AVector, double & CosVal);
	GVector operator=(GVector);
  /****************************************************
	@@	해당 Vector에대한 Unit Vector를 반환한다. 
	@@	원본이 Unit Vector로 바뀐다 
	*/
	GVector MakeUnit();
	/****************************************************
	@@	해당 Vector에대한 Unit Vector를 반환한다. 
	@@	원본은 바뀌지 않는다. 
	*/
	void Set( LDOUBLE xi,LDOUBLE yi,LDOUBLE zi);
	void Set( GPoint3D StartPoint,GPoint3D EndPoint);
	void Set( GPoint3D VectorPoint);
	void ZeroVector();
	
	GVector Unit();
  LDOUBLE  Abs();
	LDOUBLE  Length();
	GVector        Scaling( LDOUBLE S);
  GVector        Trans( GVector );
  GVector        Trans( LDOUBLE X,LDOUBLE Y,LDOUBLE Z);

	CString GetString(CString strFormat = "[ %f ][ %f ][ %f ]");
};
#endif