// GClipPolybyVatti.h: interface for the GClipPolybyVatti class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(__GCLIPPOLYBYVATTI_H__)
#define __GCLIPPOLYBYVATTI_H__

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StructDef.h"
#include "GPoint3D.h"
#include "GLine3D.h"

//--------------------------------------------------------------------------------
/* Increase GPC_EPSILON to encourage merging of near coincident edges    */

//#define GPC_EPSILON (DBL_EPSILON)
#define GPC_EPSILON (1e-7)

#define GPC_VERSION "2.31"


/*
===========================================================================
                           Public Data Types
===========================================================================
*/

typedef enum                        /* Set operation type                */
{
  GPC_DIFF,                         /* Difference                        */
  GPC_INT,                          /* Intersection                      */
  GPC_XOR,                          /* Exclusive or                      */
  GPC_UNION                         /* Union                             */
} gpc_op;

typedef struct                      /* Polygon vertex structure          */
{
  double              x;            /* Vertex x component                */
  double              y;            /* vertex y component                */
} gpc_vertex;

class  gpc_vertex_list              /* Vertex list structure             */
{
public:
  int                 num_vertices; /* Number of vertices in list        */
  gpc_vertex         *vertex;       /* Vertex array pointer              */

  gpc_vertex_list()
  {
    num_vertices = 0;
    vertex       = 0;
  }

  ~gpc_vertex_list()
  {

  }

};

class gpc_polygon                   /* Polygon set structure             */ 
{
public:
  int                 num_contours; /* Number of contours in polygon     */
  int                *hole;         /* Hole / external contour flags     */
  gpc_vertex_list    *contour;      /* Contour array pointer             */

  gpc_polygon()
  {
    num_contours = 0;
    hole         = 0;
    contour      = 0;
  }
  ~gpc_polygon()
  {

  }
} ;

class  gpc_tristrip                 /* Tristrip set structure            */
{
public:
  int                 num_strips;   /* Number of tristrips               */
  gpc_vertex_list    *strip;        /* Tristrip array pointer            */

  gpc_tristrip()
  {
    num_strips = 0;
    strip      = 0;
  }
  
};


/*
===========================================================================
                       Public Function Prototypes
===========================================================================
*/

void gpc_read_polygon        (FILE            *infile_ptr, 
                              int              read_hole_flags,
                              gpc_polygon     *polygon);

void gpc_write_polygon       (FILE            *outfile_ptr,
                              int              write_hole_flags,
                              gpc_polygon     *polygon);

void gpc_add_contour         (gpc_polygon     *polygon,
                              gpc_vertex_list *contour,
                              int              hole);

void gpc_polygon_clip        (gpc_op           set_operation,
                              gpc_polygon     *subject_polygon,
                              gpc_polygon     *clip_polygon,
                              gpc_polygon     *result_polygon);

void gpc_tristrip_clip       (gpc_op           set_operation,
                              gpc_polygon     *subject_polygon,
                              gpc_polygon     *clip_polygon,
                              gpc_tristrip    *result_tristrip);

void gpc_polygon_to_tristrip (gpc_polygon     *polygon,
                              gpc_tristrip    *tristrip);

void gpc_free_polygon        (gpc_polygon     *polygon);

void gpc_free_tristrip       (gpc_tristrip    *tristrip);


class CuttingEdge
{
public:
  GPoint3D m_Pos1;
  GPoint3D m_Pos2;

  bool IsValid()
  {
    if(m_Pos1 == m_Pos2) return false;
    return true;
  }
};

class GClipPolyBy_Vatti_algorithm  
{
public:
  gpc_polygon m_SPolys;// Subject Polys 
  gpc_polygon m_CPolys;// Cutter Polys
  gpc_polygon m_RPolys;// Results
public:
  GClipPolyBy_Vatti_algorithm ();
  ~GClipPolyBy_Vatti_algorithm();
  void ResetAll             ();
  void ResetSubjectPolys    ();
  void ResetClipPolys       ();
  void ResetResultPolys     ();
  
  void AddSubjectPoly       (CArray<GPoint3D, GPoint3D&>& Poly,BOOL bIsHole);
  void AddSubjectPoly       (CArray<DGL_3dp,DGL_3dp&>& Poly,BOOL bIsHole);
  void AddClipPoly          (CArray<DGL_3dp,DGL_3dp&>& Poly);
  void AddClipPoly          (CArray<GPoint3D, GPoint3D&>& Poly);
  void AddSubjectPoly       (CArray<GPoint3D, GPoint3D> & Poly,BOOL bIsHole);
  void AddClipPoly          (CArray<GPoint3D, GPoint3D> & Poly);
  void MakeResult_intersect ();
  void MakeResult_merge     ();
  void MakeResult_extract   ();
  
  int  GetResultPolyIndex(BOOL bIsHole,CArray<int,int>&arIPoly);
  BOOL GetResultPoly     (int nIndex,CArray<GPoint3D, GPoint3D&> &Poly);
  BOOL GetResultPoly     (int nIndex,CArray<GPoint3D, GPoint3D> &Poly);
  BOOL GetResultPoly     (int nIndex,CArray<DGL_3dp,DGL_3dp&>& Poly);
    
  //-----------------------------------------------------------------------------------------------------
  // 주어진 Point가 Line 범위를 벗어나면 오류처리 한다. 
  bool DistFromLineToPoint(GPoint3D &LP1,GPoint3D &LP2,GPoint3D &FromP,double & Dist);
  bool IsWithin2DByAngle(gpc_vertex_list &pPoly, GPoint3D point,bool& bIsOnEdge);
  bool GetLineIntersectPoint2D(GLine3D &LToTest, GLine3D& TempL, GPoint3D &IntersectPoint);
  bool GetCuttingEdge(CArray<CuttingEdge,CuttingEdge&> &Edges,GLine3D& CuttingLine);
  
};

/*
#include "GViewWnd.h"

class GData_Geom;
class GViewWnd;
class VattisPolyClipEventHandler : public I_VWEventHandler
{
public:
  CArray<GPoint3D,GPoint3D>* m_pCurPolyV;
  CArray<CArray<GPoint3D,GPoint3D>*,CArray<GPoint3D,GPoint3D>*> m_arOPolys;
  CArray<CArray<GPoint3D,GPoint3D>*,CArray<GPoint3D,GPoint3D>*> m_arIPolys;
  CArray<CArray<GPoint3D,GPoint3D>*,CArray<GPoint3D,GPoint3D>*> m_arCPolys;
  CArray<CArray<GPoint3D,GPoint3D>*,CArray<GPoint3D,GPoint3D>*> m_arRPolys; // Result Polygon 
  CArray<CArray<GPoint3D,GPoint3D>*,CArray<GPoint3D,GPoint3D>*> m_arHPolys; // Result Hole 
  CArray<CuttingEdge,CuttingEdge&>                         m_arCuttingEdge ;


  int m_nOperation; // (0) Intersect (1) Merge (2) Extract
  //------------------------------------------------------------------------------------------
  // (0) Input OPoly (1) Input IPolys (3) Input CuttingLine  (4) Display Result 
  int m_nMode; 

  VattisPolyClipEventHandler();
  virtual ~VattisPolyClipEventHandler();
  
  void SetOperationIntersect(){ m_nOperation = 0;}
  void SetOperationMerge    (){ m_nOperation = 1;}
  void SetOperationExtract  (){ m_nOperation = 2;}
  BOOL IsOperationIntersect (){ return m_nOperation == 0;}
  BOOL IsOperationMerge     (){ return m_nOperation == 1;}
  BOOL IsOperationExtract   (){ return m_nOperation == 2;}

  void StartSPolyInput (BOOL bIsHole);
  void StartCPolyInput ();
  void RemoveCPolys    ();
  void RemoveSPolys    (BOOL bIsHole);
  void RemoveRPolys    ();
  void ResetAllData    ();
  void MakeResult      (GViewWnd* pVW);

  BOOL IsSamePoint(GPoint3D* pP1,GPoint3D* pP2, int Tol,GViewWnd* pVW);
  //------------------------------------------------------------------------------------------
  // Event Handler !!!
  virtual void On_LButtonDown(GPoint3D *pWcsPos,GViewWnd* pVW);
  //------------------------------------------------------------------------------------------
  // Paint Handler !!!
  virtual void On_Paint      (GViewWnd* pVW,CDC*pDC);
protected:
  bool GetBoundingBox(GPoint3D& MinP, GPoint3D& MaxP);
};

class GVattisPolyImporter
{
protected:
  VattisPolyClipEventHandler* m_pVattiEVHandler;
public:
  GVattisPolyImporter(VattisPolyClipEventHandler* pVEVHandler)
  {
    m_pVattiEVHandler = NULL;
  }
  
  BOOL ImportOutterPolys(GData_Geom* pDataGeom);
  BOOL ImportInnerPolys (GData_Geom* pDataGeom);
  BOOL ImportCutterPolys(GData_Geom* pDataGeom);
  BOOL ImportResultPolys(GData_Geom* pDataGeom);
  BOOL ImportAllPolys   (GData_Geom* pDataGeom);
};
*/
#endif // !defined(__GCLIPPOLYBYVATTI_H__)
