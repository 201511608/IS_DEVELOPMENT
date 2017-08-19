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

  static void mathRotateX(double angle, double& rx, double& ry, double& rz);  //Global-X�࿡ ���� ȸ��, angle:ȸ����[deg], rx,ry,rz:ȸ����� ��ǥ
  static void mathRotateY(double angle, double& rx, double& ry, double& rz);  //Global-Y�࿡ ���� ȸ��, angle:ȸ����[deg], rx,ry,rz:ȸ����� ��ǥ
  static void mathRotateZ(double angle, double& rx, double& ry, double& rz);  //Global-Z�࿡ ���� ȸ��, angle:ȸ����[deg], rx,ry,rz:ȸ����� ��ǥ
  static void mathRotate(double angle, double ux, double uy, double uz, double& rx, double& ry, double& rz);  //���� �࿡ ���� ȸ��, angle:ȸ����[deg], ux,uy,uz:ȸ���� ����, rx,ry,rz:ȸ����� ��ǥ
  static void mathRotate(double angle, double px, double py, double pz, double ux, double uy, double uz,
               double& rx, double& ry, double& rz); //���� �࿡ ���� ȸ��, angle:ȸ����[deg], px,py,pz:ȸ������� ����, ux,uy,uz:ȸ���� ����, rx,ry,rz:ȸ����� ��ǥ
  static double mathDot(double xyz1[3], double xyz2[3]);
  static double mathDot(double u1x, double u1y, double u2x, double u2y);
  static void mathCross(double xyz1[3], double xyz2[3], double xyzout[3]);
	static double mathCrossAngleNormalize(double Vector1[3], double Vector2[3]);
  static double mathCrossAngleNormalizeWithSign(double Vector1[3], double Vector2[3]);	
  static double mathCrossAngle(double Vector1[3], double Vector2[3]);                                 //�κ��� ������ ������ ���
  static double mathCrossAngle2D(double u1x, double u1y, double u2x, double u2y);                     //normalize  �� �� 2�������� ������ ������ ���
  static double mathCrossAngle2D(double Vector1[3], double Vector2[3]);                                 //�� 2�������� ������ ������ ��� / +-180�� ������.
  static double mathLength(double dx, double dy, double dz=0.);                                       //1���� �������� �Ÿ�(2,3����)
  static double mathLength(double dxi, double dyi, double dxj, double dyj);                           //2������ �Ÿ�(2����)
  static double mathLength(double dxi, double dyi, double dzi, double dxj, double dyj, double dzj);   //2������ �Ÿ�(3����)
  static double mathArea(double xyz1[3], double xyz2[3], double xyzout[3]);                           //�ﰢ�� ����
  static double mathArea2d(double x1, double y1, double x2, double y2, double x3, double y3);         //�ﰢ�� ����(2����)
  static double mathArea2d(int nPoint, double x[], double y[]);
  static double mathVolume(int nVertex, double xyz[][3]);                                             // 4,6,8 ���� �������� ���� ��ü����� ü��
  static void mathNormal(double Vector1[3], double Vector2[3], double VectorN[3]);                    //���������� 2���� ����
  static void mathNormal(double Vector1[3], double Vector2[3], double Vector3[3], double VectorN[3]);//Vector1������ 2���� ����
  static BOOL mathNormalize(double xyz[3], double xyzn[3]);                                           //�������ͷ� ��ȯ
  static BOOL mathNormalize(double dx, double dy, double dz, double& dxn, double& dyn, double& dzn);  //�������ͷ� ��ȯ(3d)
  static BOOL mathNormalize(double dx, double dy, double& dxn, double& dyn);                          //�������ͷ� ��ȯ(2d)
  static void mathPLCrossPoint(double a, double b, double c, double l, double m, double n,
                 double x, double y, double z, double& rpx, double& rpy, double& rpz);  //���� ���� ��������, a,b,c:������ 1��, l,m,n:���� ���⺤��, x,y,z:���� ��ǥ, rpx,rpy,rpz:������ ��ǥ
  static void mathPLCrossPoint(double a, double b, double c, double l, double m, double n,
                 double x, double y, double z, double& rpx, double& rpy, double& rpz, double& t);  //���� ���� ��������, a,b,c:������ 1��, l,m,n:���� ���⺤��, x,y,z:���� ��ǥ, rpx,rpy,rpz:������ ��ǥ
  static void mathVec2Ang(double ux, double uy, double uz, double& rx, double& ry);  //���͸� ������ ��ȯ(Z���� Y������ ry, X������ rx ȸ����Ű�� ux,uy,uz��ġ�� ��)
  static BOOL mathPlaneEquation(double p1[3], double p2[3], double p3[3], double& a, double& b, double& c, double & h);
  static double mathDistanceToPoint(double p1[3], double p2[3]);
  static double mathDistanceToLine(double line_i[3], double line_j[3], double point[3]);// ���� ������ �ִܰŸ�(�� ����������)
  static double mathDistanceToLine2D(double line_i_x, double line_i_y, double line_j_x, double line_j_y, double point_x, double point_y);// ���� ������ �ִܰŸ�(�� ����������)
  static double mathDistanceToLine2DEx(double xi, double yi, double xj, double yj, double x, double y); // ���� �������� �ִܰŸ�
  static double mathDistanceToPlane(double a, double b, double c, double h, double x, double y, double z); // ���� ���� �����Ÿ�
  static double mathDistanceToPlaneWithSign(double a, double b, double c, double h, double x, double y, double z); // ���� ���� �����Ÿ� ��ȣ ����

  static void   mathCross_product_3d(const double xyz0[3], const double xyz1[3], const double xyz2[3], double result_vector[3]);

  /////////////////////////////////////////////////
  // Add by sshan(090708)
  // ������ ���� ���Կ��� ��ȯ(2D)  
  static BOOL mathIncludePointInLine(double xi, double yi, double xj, double yj, double x, double y, double Tol = 0.);
  /////////////////////////////////////////////////
  // Add by sshan(090708)
  // ��鿡 ���� ���Կ��� ��ȯ
  // p1, p2, p3 : Plane�� ������ ������
  // point : ���Կ��� �Ǵ��� ��
  static BOOL mathIncludePointInPlane(double p1[3], double p2[3], double p3[3], double point[3], double Tol = 0.);

  // Add by ZINU.('02.09.26). ���� ���������� ���� �Ÿ�
  static double mathDistanceFromIntersectPointToLine(double line_i[3], double line_j[3], double point[3]);
  // ������ ������ �������
  static BOOL mathIntersectPointToPlane(double LineVector[3], double LinePoint[3], double PlaneNormal[3], double PlanePoint[3], double pInts[3]);
  // ������ ������ �������, dDistance:������ �������� ������쿡�� 0. �ƴϸ� �������� �Ÿ�
  static BOOL mathIntersectPointToPlane(double LinePoint1[3], double LinePoint2[3], double PlaneNormal[3], double PlanePoint[3], double pInts[3], double& dDistance);
  static double mathDistanceLine(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3]); // �ΰ��� ������ �ִ� �Ÿ�
  static BOOL mathIntersectLine(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3], double Tolerance, double& dblDistance, double pInts[3]); // �ΰ��� ������ ����
  static BOOL mathIntersectLine2(double pl1_i[3], double pl1_j[3], double pl2_i[3], double pl2_j[3], double Tolerance, double& dblDistance, double pInts[3]); // �ΰ��� ���� ���� ����
  static BOOL mathIntersectLine2D(double p1i_x, double p1i_y, double p1j_x, double p1j_y, double p2i_x, double p2i_y, double p2j_x, double p2j_y, double& pInts_x, double& pInts_y, BOOL bCheckOnLine=FALSE); // �ΰ��� ������ ����
  // ���� ��ǥ��� ��ǥ��ȯ(ux,uy,uz�� ���� ����, �� �����̵��� ����)
  static void mathTranUCS(const double ux[3], const double uy[3], const double uz[3], const int nData, double coor[][3]);
  static void mathGCS2UCS(double& x, double& y, double& z, double ucs[3][3]);
  static void mathUCS2GCS(double& x, double& y, double& z, double ucs[3][3]);
  static void mathGCS2UCS(double& x, double& y, double ucs[2][2]);
  static void mathUCS2GCS(double& x, double& y, double ucs[2][2]);
  static BOOL mathPolyCentroid(int n, double x[], double y[], double& xCentroid, double& yCentroid, double& area);
  // ���� ��ǥ���� �� ���� ���� ��ǥ��� �ٲٴ� �Լ�(x, y, z) => (radius, angle, height)   
  static void mathConvertToCylinderCoordi(double node_xyz[3], double cyn_xyz[3], double org_xyz[3], double rot_xyz[3], double pol_xyz[3]);
                                              // ��� ��ǥ    // ��ȯ�� ���� ��ǥ    // ������ǥ    // ȸ���� ���� ��ǥ// ���� ���� ��ǥ        
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

  // �Ǽ����� ���� �ݿø� (val:���, point:��ġ(�Ҽ��� �Ʒ�)
  static double mathRoundOff(const double &val, const int &point);

  // SHIN (2006.2.3) 2���� ��������� ���� �߰���
  //normalize  �� �� 2�������� ������ ������ ��� +-180�� ������ ���� ����
  static double mathCrossAngle2DSign(double u1x, double u1y, double u2x, double u2y);
  // �� 2���� vector�� ������ ���Ͽ� �Ѱ��� 
  static double mathCross2D(double vector1[2], double vector2[2]);
  // 2���� vector�� ũ�⸦ ����vector�� ��ȯ
  //   vector  : ��ȯ�� ����
  //   vectorn : Normalize�� ���͸� �Ѱܹ��� ����
  static BOOL mathNormalize2D(double vector[2], double vectorn[2]); 
  // 2���� ��ǥ p1�� p2�� ���� ���� �ٲپ��� 
  static void mathSwap2D( double p1[2], double p2[2] );
  // 2���� Line���� �ش���ǥ�� ���ԵǴ��� ���θ� �Ǵ���
  //   bound1  : Line�� ���� ��ǥ
  //   bound2  : Line�� ���� ��ǥ
  //   targetPt: ���Կ��θ� Ȯ���� ��ǥ
  //   isOnLine: bound��ǥ�� targetPt��ǥ�� ��ġ�Ұ�� �������� ���� ����
  //   return  : TRUE=����   FALSE=������
  static BOOL mathIsPointOfLine2D( double bound1[2], double bound2[2], double targetPt[2], BOOL isOnLine );
  // �ش� ��ǥ(p1)�� 2���� polyLine���� ������ ���θ� �Ǵ��� 
  //   p1      : ���Կ��θ� Ȯ���� ��ǥ
  //   nData   : polyLine�� ��ǥ ����
  //   polyLine: ���� 2���� ��ǥ�� �迭
  //   bIncludeOutLine : �ش���ǥ�� Line �� �������� ��ġ�Ұ�� �������� ���� ���� 
  //   return  : TRUE=����   FALSE=������
  // ��������(���� �˰��� P1124)
  static BOOL mathIsInsidePoint2D(double p1[2], const int nData, double polyLine[][2], BOOL bIncludeOutLine );
  // �� ������ �����ϴ��� ���θ� �Ǵ���
  //   p1Org : ù��° ������ ������ǥ
  //   p2Org : ù��° ������ ������ǥ
  //   p3Org : �ι�° ������ ������ǥ
  //   p4Org : �ι�° ������ ������ǥ
  //   return: 1=�����ϴ� ���  -1=�������� �ʴ°��  0=�Ѽ����� ������ �ٸ����п� ���ԵǴ� ���  
  // ��������(���� �˰��� P1118)
  static int  mathIntersect_ccw2D( double p1Org[2], double p2Org[2], double p3Org[2], double p4Org[2] );

  static BOOL mathIsInsidePoint2D_Tol(double p1[2], const int nData, double polyLine[][2], double dTol, BOOL bIncludeOutLine=TRUE);

  // ������ ������ ���� ���ϱ�
  //   line1 : Line1�� �����ϴ� 2���� ��ǥ�� (line1[1][0] 2��° ��ǥ�� x��, line1[1][1] 2��° ��ǥ�� y��)
  //   line2 : Line2�� �����ϴ� 2���� ��ǥ��
  //   cross : �������� ��ǥ
  //   return: 0:����  1:����
  static int mathLineLineCross2D(double line1[2][2], double line2[2][2], double cross[2]);
  // ������ ������ ���� ���ϱ�
  //   line1 : Line�� �����ϴ� 2���� ��ǥ�� (line1[1][0] 2��° ��ǥ�� x��, line1[1][1] 2��° ��ǥ�� y��)
  //   sttP  : ������ �����ϴ� ��ǥ��
  //   endP  : ������ �����ϴ� ��ǥ��
  //   cross : �������� ��ǥ
  //   return: 0:�������� ���г��� ������  1:�������� ���г��� ������
  static int mathLineSegCross2D(double line[2][2], double sttP[2], double endP[2], double cross[2]);

	// polyLine���ο� ���ԵǴ� ������ ���̸� ����Ͽ� �Ѱ���
	//   line1 : Line�� �����ϴ� 2���� ��ǥ�� (line1[1][0] 2��° ��ǥ�� x��, line1[1][1] 2��° ��ǥ�� y��)
	//   nData   : polyLine�� ��ǥ ����
  //   polyLine: ���� 2���� ��ǥ�� �迭
  //   return  : polyLine���ο� ���ԵǴ� ������ ����
	static double mathsInsideLength(double line[2][2], const int nData, double polyLine[][2]);

	// ������ : ���� �������� �ʾ��� SHIN.07.11.12
	// 2���� polyline�� dOffset�Ÿ���ŭ ������ �ٱ���, �������� ũ�⸦ ������
	//   dOffset : Offset�Ÿ�(-:������ �������� Offset,  +:������ �ٱ������� Offset)
	//   nData   : polyLine�� ��ǥ ����
  // < polyLine: �Էµ� polyLine�� Offset���Ѽ� �Ѱ���
	//   nRotType: polyLine�� ȸ������(0:Auto,  1:�ݽð����,  2:�ð����)
	static BOOL mathOffsetOfPolyline(double dOffset, const int nData, double polyLine[][2], int nRotType = 0);

	// 3������ ���� ���� ������ ���
	//   P1, P2, P3 : ȣ�� �����ϴ� 3��(1,2,3����)
	// < CenterP    : ���� �߽���
	// < Radius     : ������
	//   RETURN     : -1:3���� ��� ���� ���, 0:3���� �����ΰ��, 1:��갡���� ���
	static int  mathCircleForm3Point(double p1[2], double p2[2], double p3[2], double CenterP[2], double& dRadius);
	// 3������ ���� ȣ�� ������ ���
	//   P1, P2, P3 : ȣ�� �����ϴ� 3��(1,2,3����)
	// < CenterP    : ���� �߽���
	// < Radius     : ������
	// < StartAngle : ���۰���
	// < SweepAngle : ���ΰ���(�ݽð���� +) 
	//   RETURN     : -1:3���� ��� ���� ���, 0:3���� �����ΰ��, 1:��갡���� ���
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
  // gauss �ҰŹ� ���1.
  static bool gauss_elimination1(const int n, double **A, double *B, double *X);
  static int  invrs_matrix(const int n, double** pm1, double **r);
  static double det(const int n, double **pm1);
  static void mult_vector(const int n, double **matrix, double *v, double *r);
  // gauss �ҰŹ� ���2.
  static bool gauss_elimination2(const int n, double **A, double *B, double *X);
	
public:
	static BOOL GetDivideNum(const double dVal, int &nDiv, int &nItr);
  static int RealToInt(const double dVal, double dTol=1e-6);
	
};

#endif // !defined(__MATHFUNC_H__)
