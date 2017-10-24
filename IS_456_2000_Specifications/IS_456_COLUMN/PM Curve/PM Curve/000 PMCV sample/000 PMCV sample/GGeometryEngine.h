#ifndef __GGEOMETRYENGINE_H__
#define __GGEOMETRYENGINE_H__

class GGeometryEngine
{  
public:  
  GGeometryEngine();
	virtual ~GGeometryEngine();
public:	
	LBOOL _IsWithin2D(GRect rect, GPoint3D * pPoint);

	BOOL     NewellsMethod(CArray<GPoint3D,GPoint3D&>&VArr,GVector&NorVector);
	BOOL     LeftInnerIntersectPoint(GPoint3D P1,GPoint3D PC, GPoint3D P2,GVector NorVector,double Dist,GPoint3D& IPoint);
	BOOL     PolygonShrink_loveme(CArray<GPoint3D,GPoint3D&>&VArr,CArray<GPoint3D,GPoint3D&>&ShrinkedVArr,
                            double SFactorOrDist,BOOL bIsFactor); 
	// 무한직선으로 간주 
  BOOL     GetIntersectLine3DPoint(GLine3D & L1,GLine3D & L2,GPoint3D & IPos);
};
#endif

