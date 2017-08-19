// MathFunc.h: interface for the CMathFunc class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(__MATHFUNC_H__)
#define __MATHFUNC_H__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <math.h>
#include <afxtempl.h>

typedef double _matrix[4][4];
typedef double _vector[4];

///////////////////////////////////////
class CMathFunc
{
public:
  CMathFunc(){};
  virtual ~CMathFunc(){};
  static double m_pi;
  static double m_trrad;
  static double m_trang;
  static double m_NormalZero;

	static double GetAngle(double x1, double y1, double x2, double y2, BOOL bPosAngle = TRUE);

  static int math_ccw(double ax, double ay, double bx, double by, double cx, double cy);

  static void mathMxIdentity(_matrix m);
  static void mathMxMult(_matrix m1, _matrix m2, _matrix mout);
  static void mathMxCopy(_matrix min, _matrix mout);

  static void mathVecMult(_vector v, _matrix m, _vector vout);
  static double mathVecDot(_vector v1, _vector v2);
  static void mathVecCross(_vector v1, _vector v2, _vector vout);

  static void mathRotateX(double angle, double& rx, double& ry, double& rz);  //Global-X축에 대한 회전, angle:회전각[deg], rx,ry,rz:회전대상 좌표
  static void mathRotateY(double angle, double& rx, double& ry, double& rz);  //Global-Y축에 대한 회전, angle:회전각[deg], rx,ry,rz:회전대상 좌표
  static void mathRotateZ(double angle, double& rx, double& ry, double& rz);  //Global-Z축에 대한 회전, angle:회전각[deg], rx,ry,rz:회전대상 좌표
  static void mathRotate(double angle, double ux, double uy, double uz, double& rx, double& ry, double& rz);  //임의 축에 대한 회전, angle:회전각[deg], ux,uy,uz:회전축 벡터, rx,ry,rz:회전대상 좌표
  static void mathRotate(double angle, double px, double py, double pz, double ux, double uy, double uz,
               double& rx, double& ry, double& rz); //임의 축에 대한 회전, angle:회전각[deg], px,py,pz:회전축상의 한점, ux,uy,uz:회전축 벡터, rx,ry,rz:회전대상 좌표
  static double mathDot(double xyz1[3], double xyz2[3]);
  static double mathDot(double u1x, double u1y, double u2x, double u2y);
  static void mathCross(double xyz1[3], double xyz2[3], double xyzout[3]);
	static double mathCrossAngleNormalize(double Vector1[3], double Vector2[3]);
  static double mathCrossAngleNormalizeWithSign(double Vector1[3], double Vector2[3]);	
  static double mathCrossAngle(double Vector1[3], double Vector2[3]);                                 //두벡터 사이의 각도를 계산
  static double mathCrossAngle2D(double u1x, double u1y, double u2x, double u2y);                     //normalize  된 두 2차원벡터 사이의 각도를 계산
  static double mathCrossAngle2D(double Vector1[3], double Vector2[3]);                                 //두 2차원벡터 사이의 각도를 계산 / +-180도 범위로.
  static double mathLength(double dx, double dy, double dz=0.);                                       //1점과 원점과의 거리(2,3차원)
  static double mathLength(double dxi, double dyi, double dxj, double dyj);                           //2점간의 거리(2차원)
  static double mathLength(double dxi, double dyi, double dzi, double dxj, double dyj, double dzj);   //2점간의 거리(3차원)
  static double mathArea(double xyz1[3], double xyz2[3], double xyzout[3]);                           //삼각형 면적
  static double mathArea2d(double x1, double y1, double x2, double y2, double x3, double y3);         //삼각형 면적(2차원)
  static double mathArea2d(int nPoint, double x[], double y[]);
  static double mathVolume(int nVertex, double xyz[][3]);                                             // 4,6,8 개의 꼭지점을 가진 고체요소의 체적
  static void mathNormal(double Vector1[3], double Vector2[3], double VectorN[3]);                    //원점기준의 2선의 법선
  static void mathNormal(double Vector1[3], double Vector2[3], double Vector3[3], double VectorN[3]);//Vector1기준의 2선의 법선
  static BOOL mathNormalize(double xyz[3], double xyzn[3]);                                           //단위벡터로 변환
  static BOOL mathNormalize(double dx, double dy, double dz, double& dxn, double& dyn, double& dzn);  //단위벡터로 변환(3d)
  static BOOL mathNormalize(double dx, double dy, double& dxn, double& dyn);                          //단위벡터로 변환(2d)
  static void mathPLCrossPoint(double a, double b, double c, double l, double m, double n,
                 double x, double y, double z, double& rpx, double& rpy, double& rpz);  //선과 점의 수직교점, a,b,c:선상의 1점, l,m,n:선의 방향벡터, x,y,z:점의 좌표, rpx,rpy,rpz:교점의 좌표
  static void mathPLCrossPoint(double a, double b, double c, double l, double m, double n,
                 double x, double y, double z, double& rpx, double& rpy, double& rpz, double& t);  //선과 점의 수직교점, a,b,c:선상의 1점, l,m,n:선의 방향벡터, x,y,z:점의 좌표, rpx,rpy,rpz:교점의 좌표
  static void mathVec2Ang(double ux, double uy, double uz, double& rx, double& ry);  //벡터를 각도로 변환(Z축을 Y축으로 ry, X축으로 rx 회전시키면 ux,uy,uz위치가 됨)
  static BOOL mathPlaneEquation(double p1[3], double p2[3], double p3[3], double& a, double& b, double& c, double & h);
  static double mathDistanceToPoint(double p1[3], double p2[3]);
  static double mathDistanceToLine(double line_i[3], double line_j[3], double point[3]);// 점과 선과의 최단거리(선 범위내에서)
  static double mathDistanceToLine2D(double line_i_x, double line_i_y, double line_j_x, double line_j_y, double point_x, double point_y);// 점과 선과의 최단거리(선 범위내에서)
  static double mathDistanceToLine2DEx(double xi, double yi, double xj, double yj, double x, double y); // 점과 직선과의 최단거리
  static double mathDistanceToPlane(double a, double b, double c, double h, double x, double y, double z); // 점과 면의 수직거리
  static double mathDistanceToPlaneWithSign(double a, double b, double c, double h, double x, double y, double z); // 점과 면의 수직거리 부호 포함

  static void   mathCross_product_3d(const double xyz0[3], const double xyz1[3], const double xyz2[3], double result_vector[3]);

  /////////////////////////////////////////////////
  // Add by sshan(090708)
  // 직선에 점의 포함여부 반환(2D)  
  static BOOL mathIncludePointInLine(double xi, double yi, double xj, double yj, double x, double y, double Tol = 0.);
  /////////////////////////////////////////////////
  // Add by sshan(090708)
  // 평면에 점의 포함여부 반환
  // p1, p2, p3 : Plane을 형성할 포인터
  // point : 포함여부 판단할 점
  static BOOL mathIncludePointInPlane(double p1[3], double p2[3], double p3[3], double point[3], double Tol = 0.);

  // Add by ZINU.('02.09.26). 점과 직선사이의 수직 거리
  static double mathDistanceFromIntersectPointToLine(double line_i[3], double line_j[3], double point[3]);
  // 직선과 평면과의 교점계산
  static BOOL mathIntersectPointToPlane(double LineVector[3], double LinePoint[3], double PlaneNormal[3], double PlanePoint[3], double pInts[3]);
  // 직선과 평면과의 교점계산, dDistance:교점이 직선내에 있을경우에는 0. 아니면 교점과의 거리
  static BOOL mathIntersectPointToPlane(double LinePoint1[3], double LinePoint2[3], double PlaneNormal[3], double PlanePoint[3], double pInts[3], double& dDistance);
  static double mathDistanceLine(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3]); // 두개의 선분의 최단 거리
  static BOOL mathIntersectLine(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3], double Tolerance, double& dblDistance, double pInts[3]); // 두개의 선분의 교점
  static BOOL mathIntersectLine2(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3], double Tolerance, double& dblDistance, double pInts[3]); // 두개의 선분 내의 교점
  static BOOL mathIntersectLine2D(double p1i_x, double p1i_y, double p1j_x, double p1j_y, double p2i_x, double p2i_y, double p2j_x, double p2j_y, double& pInts_x, double& pInts_y, BOOL bCheckOnLine=FALSE); // 두개의 선분의 교점
  // 임의 좌표계로 좌표변환(ux,uy,uz는 단위 벡터, 단 원점이동은 없음)
  static void mathTranUCS(const double ux[3], const double uy[3], const double uz[3], const int nData, double coor[][3]);
  static void mathGCS2UCS(double& x, double& y, double& z, double ucs[3][3]);
  static void mathUCS2GCS(double& x, double& y, double& z, double ucs[3][3]);
  static void mathGCS2UCS(double& x, double& y, double ucs[2][2]);
  static void mathUCS2GCS(double& x, double& y, double ucs[2][2]);
  static BOOL mathPolyCentroid(int n, double x[], double y[], double& xCentroid, double& yCentroid, double& area);
  // 직교 좌표계의 한 점을 원통 좌표계로 바꾸는 함수(x, y, z) => (radius, angle, height)   
  static void mathConvertToCylinderCoordi(double node_xyz[3], double cyn_xyz[3], double org_xyz[3], double rot_xyz[3], double pol_xyz[3]);
                                              // 대상 좌표    // 변환된 원통 좌표    // 극점좌표    // 회전축 위의 좌표// 극축 위의 좌표        
  static bool project_on_cylinder_normal(const double node_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3]);
  static bool project_on_cylinder_vector(const double node_xyz[3], const double vector_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3]);
  static bool project_on_cone_normal(const double node_xyz[3], const double radius1, const double radius2, const double height, const double tolerance, double projected_node_xyz[3]);
  static bool project_on_cone_vector(const double node_xyz[3], const double vector_xyz[3], const double radius1, const double radius2, const double height, const double tolerance, double projected_node_xyz[3]);
  static bool project_on_sphere_normal(const double node_xyz[3], const double radius, const double tolerance, double projected_node_xyz[3]);
  static bool project_on_sphere_vector(const double node_xyz[3], const double vector_xyz[3], const double radius, const double tolerance, double projected_node_xyz[3]);
  static bool project_on_ellipsoid_normal(const double node_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3]);
  static bool project_on_ellipsoid_vector(const double node_xyz[3], const double vector_xyz[3], const double radius, const double height, const double tolerance, double projected_node_xyz[3]);
  static bool project_on_plane_normal(const double node_xyz[3], const double plane_xyz1[3], const double plane_xyz2[2], const double plane_xyz3[3], const double tolerance, double projected_node_xyz[3]);
  static bool project_on_plane_vector(const double node_xyz[3], const double vector_xyz[3], const double plane_xyz1[3], const double plane_xyz2[2], const double plane_xyz3[3], const double tolerance, double projected_node_xyz[3]);
  static bool project_on_plane_vector_for_CuttingPlane(const double node_xyz[3], const double vector_xyz[3], const double plane_xyz1[3], const double plane_xyz2[2], const double plane_xyz3[3],double projected_node_xyz[3]);
  static bool project_on_quad_curve_normal(const double curve_xyz_1[3], const double curve_xyz_2[3], const double curve_xyz_3[3], const double point[3], const double tolerance, double project_point[3]);
  static bool project_on_quad_curve_vector(const double curve_xyz_1[3], const double curve_xyz_2[3], const double curve_xyz_3[3], const double point[3], const double direction_vector[3], const double tolerance, double project_point[3]);
  static bool check_quad_curve(const double curve_xyz_1[3], const double curve_xyz_2[3], const double curve_xyz_3[3], const double ZERO);

  static double mathInterpolate(double dLeft, double dRight, double dDistRatioFromLeft);
  static double mathInterpolate(double dLeftBot, double dRightBot, double dLeftTop, double dRightTop, 
                                double dDistRatioFromLeft, double dDistRatioFromBot);

  // 실수값에 대한 반올림 (val:대상값, point:위치(소수점 아래)
  static double mathRoundOff(const double &val, const int &point);

  // SHIN (2006.2.3) 2차원 도형계산을 위해 추가함
  //normalize  된 두 2차원벡터 사이의 각도를 계산 +-180도 범위의 값을 가짐
  static double mathCrossAngle2DSign(double u1x, double u1y, double u2x, double u2y);
  // 두 2차원 vector의 외적을 구하여 넘겨줌 
  static double mathCross2D(double vector1[2], double vector2[2]);
  // 2차원 vector의 크기를 단위vector로 변환
  //   vector  : 변환할 벡터
  //   vectorn : Normalize된 벡터를 넘겨받을 변수
  static BOOL mathNormalize2D(double vector[2], double vectorn[2]); 
  // 2차원 좌표 p1과 p2의 값을 서로 바꾸어줌 
  static void mathSwap2D( double p1[2], double p2[2] );
  // 2차원 Line내에 해당좌표가 포함되는지 여부를 판단함
  //   bound1  : Line의 구성 좌표
  //   bound2  : Line의 구성 좌표
  //   targetPt: 포함여부를 확인할 좌표
  //   isOnLine: bound좌표와 targetPt좌표가 일치할경우 포함으로 볼지 여부
  //   return  : TRUE=포함   FALSE=미포함
  static BOOL mathIsPointOfLine2D( double bound1[2], double bound2[2], double targetPt[2], BOOL isOnLine );
  // 해당 좌표(p1)가 2차원 polyLine내에 들어가는지 여부를 판단함 
  //   p1      : 포함여부를 확인할 좌표
  //   nData   : polyLine의 좌표 갯수
  //   polyLine: 비교할 2차원 좌표의 배열
  //   bIncludeOutLine : 해당좌표가 Line 및 꼭지점에 일치할경우 포함으로 볼지 여부 
  //   return  : TRUE=포함   FALSE=미포함
  // ※참고문헌(기하 알고리즘 P1124)
  static BOOL mathIsInsidePoint2D(double p1[2], const int nData, double polyLine[][2], BOOL bIncludeOutLine );
  // 두 선분이 교차하는지 여부를 판단함
  //   p1Org : 첫번째 선분의 구성좌표
  //   p2Org : 첫번째 선분의 구성좌표
  //   p3Org : 두번째 선분의 구성좌표
  //   p4Org : 두번째 선분의 구성좌표
  //   return: 1=교차하는 경우  -1=교차하지 않는경우  0=한선분의 끝점이 다른선분에 포함되는 경우  
  // ※참고문헌(기하 알고리즘 P1118)
  static int  mathIntersect_ccw2D( double p1Org[2], double p2Org[2], double p3Org[2], double p4Org[2] );

  static BOOL mathIsInsidePoint2D_Tol(double p1[2], const int nData, double polyLine[][2], double dTol, BOOL bIncludeOutLine=TRUE);

  // 직선과 직선의 교점 구하기
  //   line1 : Line1을 구성하는 2개의 좌표점 (line1[1][0] 2번째 좌표의 x값, line1[1][1] 2번째 좌표의 y값)
  //   line2 : Line2을 구성하는 2개의 좌표점
  //   cross : 교차점의 좌표
  //   return: 0:평행  1:교차
  static int mathLineLineCross2D(double line1[2][2], double line2[2][2], double cross[2]);
  // 직선과 선분의 교점 구하기
  //   line1 : Line을 구성하는 2개의 좌표점 (line1[1][0] 2번째 좌표의 x값, line1[1][1] 2번째 좌표의 y값)
  //   sttP  : 선분을 구성하는 좌표점
  //   endP  : 선분을 구성하는 좌표점
  //   cross : 교차점의 좌표
  //   return: 0:교차점이 선분내에 없을때  1:교차점이 선분내에 있을때
  static int mathLineSegCross2D(double line[2][2], double sttP[2], double endP[2], double cross[2]);

	// polyLine내부에 포함되는 직선의 길이를 계산하여 넘겨줌
	//   line1 : Line을 구성하는 2개의 좌표점 (line1[1][0] 2번째 좌표의 x값, line1[1][1] 2번째 좌표의 y값)
	//   nData   : polyLine의 좌표 갯수
  //   polyLine: 비교할 2차원 좌표의 배열
  //   return  : polyLine내부에 포함되는 직선의 길이
	static double mathsInsideLength(double line[2][2], const int nData, double polyLine[][2]);

	// ※주의 : 아직 검증하지 않았음 SHIN.07.11.12
	// 2차원 polyline을 dOffset거리만큼 도형의 바깥쪽, 안쪽으로 크기를 변경함
	//   dOffset : Offset거리(-:도형의 안쪽으로 Offset,  +:도형의 바깥쪽으로 Offset)
	//   nData   : polyLine의 좌표 갯수
  // < polyLine: 입력된 polyLine을 Offset시켜서 넘겨줌
	//   nRotType: polyLine의 회전방향(0:Auto,  1:반시계방향,  2:시계방향)
	static BOOL mathOffsetOfPolyline(double dOffset, const int nData, double polyLine[][2], int nRotType = 0);

	// 3점으로 부터 원의 정보를 계산
	//   P1, P2, P3 : 호를 구성하는 3점(1,2,3순서)
	// < CenterP    : 원의 중심점
	// < Radius     : 반지름
	//   RETURN     : -1:3점이 모두 같은 경우, 0:3점이 직선인경우, 1:계산가능한 경우
	static int  mathCircleForm3Point(double p1[2], double p2[2], double p3[2], double CenterP[2], double& dRadius);
	// 3점으로 부터 호의 정보를 계산
	//   P1, P2, P3 : 호를 구성하는 3점(1,2,3순서)
	// < CenterP    : 원의 중심점
	// < Radius     : 반지름
	// < StartAngle : 시작각도
	// < SweepAngle : 내부각도(반시계방향 +) 
	//   RETURN     : -1:3점이 모두 같은 경우, 0:3점이 직선인경우, 1:계산가능한 경우
	static int  mathArcForm3Point(double p1[2], double p2[2], double p3[2], double CenterP[2], double& dRadius, double& dStartAngle, double& dSweepAngle);

	static BOOL   GetCenter_Bulge(double line[2][2], double bulge, double centerP[2]);
	static double GetBulge       (double line[2][2], double radius, bool inside);

  static double mathMin(double m1, double m2);
  static double mathMin(double m1, double m2, double m3);
  static double mathMin(double m1, double m2, double m3, double m4);	

  inline static double mathSin(double x)
  {
    double r=sin(x);
    if(fabs(r) < m_NormalZero)r=0.;
    if(fabs(1.-r) < m_NormalZero)r=1.;
    if(fabs(1.+r) < m_NormalZero)r=-1.;
    return r;
  };
  inline static double mathCos(double x)
  {
    double r=cos(x);
    if(fabs(r) < m_NormalZero)r=0.;
    if(fabs(1.-r) < m_NormalZero)r=1.;
    if(fabs(1.+r) < m_NormalZero)r=-1.;
    return r;
  };
  static double mathAsin(double x)
  {
    if(x < -1)x=-1.;
    if(x >  1)x=1.;
    return asin(x);
  };
  static double mathAcos(double x)
  {
    if(x < -1)x=-1.;
    if(x >  1)x=1.;
    return acos(x);
  };
  static double mathSqrt(double x)
  {
    if(x < 0)x=0.;
    return sqrt(x);
  };
  static int mathAbs(int x)
  {
    return abs(x);
  };
  static double mathAbs(double x)
  {
    return fabs(x);
  };

protected:
  static double Solid_4(double xyz[][3]);
  static double Solid_6(double xyz[][3]);
  static double Solid_8(double xyz[][3]);

public:
  // gauss 소거법 방법1.
  static bool gauss_elimination1(const int n, double **A, double *B, double *X);
  static int  invrs_matrix(const int n, double** pm1, double **r);
  static double det(const int n, double **pm1);
  static void mult_vector(const int n, double **matrix, double *v, double *r);
  // gauss 소거법 방법2.
  static bool gauss_elimination2(const int n, double **A, double *B, double *X);
	
public:
	static BOOL GetDivideNum(const double dVal, int &nDiv, int &nItr);
  static int RealToInt(const double dVal, double dTol=1e-6);
	
};

#endif // !defined(__MATHFUNC_H__)
