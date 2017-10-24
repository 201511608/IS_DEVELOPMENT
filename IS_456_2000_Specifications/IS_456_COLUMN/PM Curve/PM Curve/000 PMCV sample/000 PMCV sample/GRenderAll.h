#ifndef __GRENDERALL_H__
#define __GRENDERALL_H__

#ifndef __AFXTEMPL_H__
	#include <afxtempl.h>
#endif

#ifndef _INC_FLOAT
	#include <float.h>
#endif

#ifndef _INC_MATH
	#include <math.h>
#endif

#ifndef __LGLTYPE_H__
	#include "LGLType.h"
#endif

#ifdef __INC_GPOINT3D_H__
  #ifndef __C3DPOINT_H__
    #include "C3DPoint.h"
  #endif
  #ifndef __GPOINT3D_H__ 
		#include "GPoint3D.h"
	#endif
	#ifndef __GPOINT2D_H__	
		#include "GPoint2D.h"
	#endif
#endif

#ifdef __INC_GPOINT2D_H__
  #ifndef __C3DPOINT_H__
    #include "C3DPoint.h"
  #endif	
  #ifndef __GPOINT3D_H__ 
		#include "GPoint3D.h"
	#endif
	#ifndef __GPOINT2D_H__	
		#include "GPoint2D.h"
	#endif
#endif

#ifdef __INC_GLINE3D_H__
  #ifndef __C3DPOINT_H__
    #include "C3DPoint.h"
  #endif	
  #ifndef __GPOINT3D_H__ 
		#include "GPoint3D.h"
	#endif
	#ifndef __GPOINT2D_H__	
		#include "GPoint2D.h"
	#endif
	#ifndef __GLINE3D_H__ 
		#include "GLine3D.h"
	#endif
	#ifndef __GLINE2D_H__
		#include "GLine2D.h"
	#endif
#endif

#ifdef __INC_GLINE2D_H__
  #ifndef __C3DPOINT_H__
    #include "C3DPoint.h"
  #endif	
  #ifndef __GPOINT3D_H__ 
		#include "GPoint3D.h"
	#endif
	#ifndef __GPOINT2D_H__	
		#include "GPoint2D.h"
	#endif
	#ifndef __GLINE3D_H__ 
		#include "GLine3D.h"
	#endif
	#ifndef __GLINE2D_H__
		#include "GLine2D.h"
	#endif
#endif

#ifdef __INC_GRECT_H__
	#ifndef __C3DPOINT_H__
    #include "C3DPoint.h"
  #endif	
  #ifndef __GPOINT3D_H__ 
		#include "GPoint3D.h"
	#endif
	#ifndef __GPOINT2D_H__	
		#include "GPoint2D.h"
	#endif
  #ifndef __GLINE3D_H__ 
		#include "GLine3D.h"
	#endif
	#ifndef __GLINE2D_H__
		#include "GLine2D.h"
	#endif
	#ifndef __GRECT_H__
		#include "GRect.h"
	#endif
#endif

#ifdef __INC_GVECTOR_H__
  #ifndef __C3DPOINT_H__
    #include "C3DPoint.h"
  #endif	
  #ifndef __GPOINT3D_H__ 
		#include "GPoint3D.h"
	#endif
	#ifndef __GPOINT2D_H__	
		#include "GPoint2D.h"
	#endif
	#ifndef __GLINE3D_H__ 
		#include "GLine3D.h"
	#endif
	#ifndef __GVECTOR_H__
		#include "GVector.h"
	#endif
#endif

#ifndef __GMATHUTIL_H__
  #include "GMathUtil.h"
#endif

#ifdef __INC_GRENDERENGINE_H__
  #ifndef __C3DPOINT_H__
    #include "C3DPoint.h"
  #endif	
	#ifndef __GPOINT3D_H__ 
		#include "GPoint3D.h"
	#endif
	#ifndef __GPOINT2D_H__	
		#include "GPoint2D.h"
	#endif
	#ifndef __GLINE3D_H__ 
		#include "GLine3D.h"
	#endif
	#ifndef __GVECTOR_H__
		#include "GVector.h"
	#endif
	#ifndef __GRECT_H__
		#include "GRect.h"
	#endif
	#ifndef __GGEOMETRYENGINE_H__
		#include "GGeometryEngine.h"
	#endif	
#endif


#endif
