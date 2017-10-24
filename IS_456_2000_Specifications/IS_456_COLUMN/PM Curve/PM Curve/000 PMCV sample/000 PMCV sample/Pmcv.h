// Pmcv.h : main header file for the PMCV application
//

#if !defined(AFX_PMCV_H__45F39745_0B0A_4E58_B003_0AE9EC0966C6__INCLUDED_)
#define AFX_PMCV_H__45F39745_0B0A_4E58_B003_0AE9EC0966C6__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// main symbols

/////////////////////////////////////////////////////////////////////////////
// CPmcvApp:
// See Pmcv.cpp for the implementation of this class
//

class CPmcvApp : public CWinApp
{
public:
	CPmcvApp();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPmcvApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// Implementation

	//{{AFX_MSG(CPmcvApp)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_PMCV_H__45F39745_0B0A_4E58_B003_0AE9EC0966C6__INCLUDED_)
