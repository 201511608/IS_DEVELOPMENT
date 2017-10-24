// PmcvDlg.h : header file
//

#if !defined(AFX_PMCVDLG_H__7DAC69B4_74C9_4C6D_99AF_C20F41D34251__INCLUDED_)
#define AFX_PMCVDLG_H__7DAC69B4_74C9_4C6D_99AF_C20F41D34251__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

/////////////////////////////////////////////////////////////////////////////
// CPmcvDlg dialog

class CPmcvDlg : public CDialog
{
// Construction
public:
	CPmcvDlg(CWnd* pParent = NULL);	// standard constructor

// Dialog Data
	//{{AFX_DATA(CPmcvDlg)
	enum { IDD = IDD_PMCV_DIALOG };
	int		m_nShape;
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CPmcvDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	//{{AFX_MSG(CPmcvDlg)
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	virtual void OnOK();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_PMCVDLG_H__7DAC69B4_74C9_4C6D_99AF_C20F41D34251__INCLUDED_)
