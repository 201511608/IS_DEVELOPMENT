// PmcvDlg.cpp : implementation file
//

#include "stdafx.h"
#include "Pmcv.h"
#include "PmcvDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	//{{AFX_DATA(CAboutDlg)
	enum { IDD = IDD_ABOUTBOX };
	//}}AFX_DATA

	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CAboutDlg)
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	//{{AFX_MSG(CAboutDlg)
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
	//{{AFX_DATA_INIT(CAboutDlg)
	//}}AFX_DATA_INIT
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CAboutDlg)
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
	//{{AFX_MSG_MAP(CAboutDlg)
		// No message handlers
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CPmcvDlg dialog

CPmcvDlg::CPmcvDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CPmcvDlg::IDD, pParent)
{
	//{{AFX_DATA_INIT(CPmcvDlg)
	m_nShape = 0;
	//}}AFX_DATA_INIT
	// Note that LoadIcon does not require a subsequent DestroyIcon in Win32
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CPmcvDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CPmcvDlg)
	DDX_Radio(pDX, IDC_RADIO1, m_nShape);
	//}}AFX_DATA_MAP
}

BEGIN_MESSAGE_MAP(CPmcvDlg, CDialog)
	//{{AFX_MSG_MAP(CPmcvDlg)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CPmcvDlg message handlers

BOOL CPmcvDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon
	
	// TODO: Add extra initialization here
	
	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CPmcvDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CPmcvDlg::OnPaint() 
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, (WPARAM) dc.GetSafeHdc(), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

#include "CalcPmcv.h"
#include <math.h>

// The system calls this to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CPmcvDlg::OnQueryDragIcon()
{
	return (HCURSOR) m_hIcon;
}

void CPmcvDlg::OnOK() 
{
  CWaitCursor wt;
	UpdateData(TRUE);

	double dPu  = 1000.0 * 1.0e3;
	double dMuy =  200.0 * 1.0e6;
	double dMuz =  200.0 * 1.0e6;
	double dMu  = sqrt(dMuy*dMuy + dMuz*dMuz);

	CCalcPmcv CalPmcv;
	CalPmcv.SetCommInputData();
	CalPmcv.SetMatlInputData(27.0, 27804.064, 300.0, 200000.0);
	
	if(m_nShape==0)
	{
		double dSize[8] = {0.0};
		dSize[0] = 1000.0; //mm
		dSize[1] = 1000.0; //mm
		
		_UMD_RC_COL_MAINRBAR MainBar;
		MainBar.dDc[0] = 100.0; //mm
		for(int i=0; i<2; ++i)
		{
			MainBar.dBarAreaNa1[i][0] = 387.1;
			MainBar.dBarAreaNa2[i][0] = 387.1;
			MainBar.dBarDiaNa1[i][0]  = 22.2;
			MainBar.dBarDiaNa2[i][0]  = 22.2;
			MainBar.iBarNum[i][0] = 10;		
		}		
		CalPmcv.SetSectInputData(DGN_SECT_SHAPE_INDEX_REG_SB, dSize, MainBar);
	}
	else if(m_nShape==1)
	{
		double dSize[8] = {0.0};	
		dSize[0] = 1000.0; //mm
		
		_UMD_RC_COL_MAINRBAR MainBar;
		MainBar.dDc[0] = 100.0; //mm		
		MainBar.dBarAreaNa1[0][0] = 387.1;
		MainBar.dBarDiaNa1[0][0]  = 22.2;
		MainBar.iBarNum[0][0] = 24;		
	
		CalPmcv.SetSectInputData(DGN_SECT_SHAPE_INDEX_REG_SR, dSize, MainBar);
	}
	else if(m_nShape==2)
	{
		double dSize[8] = {0.0};	
		dSize[0] = 1000.0; //mm
		dSize[1] = 1000.0; //mm

		// Make concrete polygon data.
		_UMD_RC_CON_POLY_DATA ConPoly;
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(0.0000000000000000, 1000.0000000000000));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(0.0000000000000000, 0.0000000000000000));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(1000.0000000000000, 0.0000000000000000));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(1000.0000000000000, 1000.0000000000000));

		double dDia1  =  22.2;
		double dArea1 = 387.1;

		// Make general rebar data.
		_UMD_RC_COL_MAINRBAR MainBar;
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(100.00000000000000, 900.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(188.88888888890000, 900.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(277.77777777780000, 900.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(366.66666666670000, 900.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(455.55555555555998, 900.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(544.44444444444002, 900.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(633.33333333330006, 900.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(722.22222222220000, 900.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(811.11111111109994, 900.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(900.00000000000000, 900.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(100.00000000000000, 100.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(188.88888888890000, 100.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(277.77777777780000, 100.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(366.66666666670000, 100.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(455.55555555555998, 100.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(544.44444444444002, 100.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(633.33333333330006, 100.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(722.22222222220000, 100.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(811.11111111109994, 100.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(900.00000000000000, 100.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(100.00000000000000, 827.27272727269997, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(100.00000000000000, 754.54545454549998, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(100.00000000000000, 681.81818181819995, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(100.00000000000000, 609.09090909090003, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(100.00000000000000, 536.36363636364001, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(100.00000000000000, 463.63636363635999, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(100.00000000000000, 390.90909090909997, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(100.00000000000000, 318.18181818180000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(100.00000000000000, 245.45454545449999, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(100.00000000000000, 172.72727272729998, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(900.00000000000000, 827.27272727269997, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(900.00000000000000, 754.54545454549998, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(900.00000000000000, 681.81818181819995, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(900.00000000000000, 609.09090909090003, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(900.00000000000000, 536.36363636364001, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(900.00000000000000, 463.63636363635999, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(900.00000000000000, 390.90909090909997, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(900.00000000000000, 318.18181818180000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(900.00000000000000, 245.45454545449999, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(900.00000000000000, 172.72727272729998, dDia1, dArea1));
		
		CalPmcv.SetSectInputDataGen(DGN_SECT_SHAPE_INDEX_REG_GEN, dSize, ConPoly, MainBar);		
	}
	else if(m_nShape==3)
	{
		// Make concrete polygon data.
		_UMD_RC_CON_POLY_DATA ConPoly;
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(0.0000000000000000, 500.00000000000006));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(7.5960000000000036, 413.17590000000001));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(30.153999999999996, 328.99000000000001));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(66.987000000000023, 250.00000000000000));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(116.97800000000001, 178.60599999999999));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(178.60599999999999, 116.97800000000001));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(250.00000000000000, 66.987000000000023));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(328.99000000000001, 30.153999999999996));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(413.17590000000001, 7.5960000000000036));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(499.99999999999994, 0.0000000000000000));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(586.82410000000004, 7.5960000000000036));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(671.00999999999999, 30.153999999999996));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(750.00000000000000, 66.987000000000023));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(821.39400000000001, 116.97800000000001));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(883.02199999999993, 178.60599999999999));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(933.01299999999992, 250.00000000000000));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(969.84600000000000, 328.99000000000001));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(992.40400000000000, 413.17590000000001));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(1000.0000000000000, 500.00000000000000));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(992.40400000000000, 586.82410000000004));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(969.84600000000000, 671.00999999999999));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(933.01299999999992, 750.00000000000000));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(883.02199999999993, 821.39400000000001));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(821.39400000000001, 883.02199999999993));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(750.00000000000000, 933.01299999999992));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(671.00999999999999, 969.84600000000000));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(586.82410000000004, 992.40400000000000));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(500.00000000000006, 1000.0000000000000));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(413.17590000000001, 992.40400000000000));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(328.99000000000001, 969.84600000000000));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(250.00000000000000, 933.01299999999992));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(178.60599999999999, 883.02199999999993));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(116.97800000000001, 821.39400000000001));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(66.987000000000023, 750.00000000000000));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(30.153999999999996, 671.00999999999999));
		ConPoly.OutPoly.aVertex.Add(_UMD_RC_GSEC_VERTEX(7.5960000000000036, 586.82410000000004));
		
		double dDia1  =  22.2;
		double dArea1 = 387.1;

		// Make general rebar data.
		_UMD_RC_COL_MAINRBAR MainBar;
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(500.00000000000000, 900.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(603.52761804099998, 886.37033051560002, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(700.00000000000000, 846.41016151379995, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(782.84271247459992, 782.84271247459992, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(846.41016151379995, 700.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(886.37033051560002, 603.52761804099998, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(900.00000000000000, 499.99999999999989, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(886.37033051560002, 396.47238195900002, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(846.41016151379995, 300.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(782.84271247459992, 217.15728752540002, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(700.00000000000000, 153.58983848619999, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(603.52761804099998, 113.62966948439998, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(499.99999999999983, 100.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(396.47238195900002, 113.62966948439998, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(300.00000000000000, 153.58983848619999, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(217.15728752540002, 217.15728752540002, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(153.58983848619999, 300.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(113.62966948439998, 396.47238195900002, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(100.00000000000000, 500.00000000000028, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(113.62966948439998, 603.52761804099998, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(153.58983848619999, 700.00000000000000, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(217.15728752540002, 782.84271247459992, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(300.00000000000000, 846.41016151379995, dDia1, dArea1));
		MainBar.arGenRbar.Add(_UMD_RC_GEN_RBAR_UNIT(396.47238195900002, 886.37033051560002, dDia1, dArea1));
		
		double dSize[8] = {0.0};	
		dSize[0] = 1000.0; //mm
		dSize[1] = 1000.0; //mm

		CalPmcv.SetSectInputDataGen(DGN_SECT_SHAPE_INDEX_REG_GEN, dSize, ConPoly, MainBar);
	}
	else ASSERT(0);

	CalPmcv.SetLoadInputData(dPu, dMuy, dMuz, dMu);
	CalPmcv.ExecuteCheck(TRUE);


	


}
