# Microsoft Developer Studio Project File - Name="Pmcv" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101

CFG=Pmcv - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Pmcv.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Pmcv.mak" CFG="Pmcv - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Pmcv - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "Pmcv - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "Pmcv - Win32 Release"

# PROP BASE Use_MFC 6
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 6
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MD /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_AFXDLL" /Yu"stdafx.h" /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /Yu"stdafx.h" /FD /c
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x412 /d "NDEBUG" /d "_AFXDLL"
# ADD RSC /l 0x412 /d "NDEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 /nologo /subsystem:windows /machine:I386
# ADD LINK32 /nologo /subsystem:windows /machine:I386

!ELSEIF  "$(CFG)" == "Pmcv - Win32 Debug"

# PROP BASE Use_MFC 6
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 6
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /Yu"stdafx.h" /FD /GZ /c
# ADD CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /D "_MBCS" /Yu"stdafx.h" /FD /GZ /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD BASE RSC /l 0x412 /d "_DEBUG" /d "_AFXDLL"
# ADD RSC /l 0x412 /d "_DEBUG" /d "_AFXDLL"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 /nologo /subsystem:windows /debug /machine:I386 /pdbtype:sept
# ADD LINK32 /nologo /subsystem:windows /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "Pmcv - Win32 Release"
# Name "Pmcv - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\Pmcv.cpp
# End Source File
# Begin Source File

SOURCE=.\Pmcv.rc
# End Source File
# Begin Source File

SOURCE=.\PmcvDlg.cpp
# End Source File
# Begin Source File

SOURCE=.\StdAfx.cpp
# ADD CPP /Yc"stdafx.h"
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\Pmcv.h
# End Source File
# Begin Source File

SOURCE=.\PmcvDlg.h
# End Source File
# Begin Source File

SOURCE=.\Resource.h
# End Source File
# Begin Source File

SOURCE=.\StdAfx.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# Begin Source File

SOURCE=.\res\Pmcv.ico
# End Source File
# Begin Source File

SOURCE=.\res\Pmcv.rc2
# End Source File
# End Group
# Begin Group "Calculation"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\CalcBasePmcv.cpp
# End Source File
# Begin Source File

SOURCE=.\CalcBasePmcv.h
# End Source File
# Begin Source File

SOURCE=.\CalcPmcv.cpp
# End Source File
# Begin Source File

SOURCE=.\CalcPmcv.h
# End Source File
# Begin Source File

SOURCE=.\MathFunc.cpp
# End Source File
# Begin Source File

SOURCE=.\MathFunc.h
# End Source File
# Begin Source File

SOURCE=.\MdgnPmTool.cpp
# End Source File
# Begin Source File

SOURCE=.\MdgnPmTool.h
# End Source File
# Begin Source File

SOURCE=.\MdgnSectTool.cpp
# End Source File
# Begin Source File

SOURCE=.\MdgnSectTool.h
# End Source File
# Begin Source File

SOURCE=.\PolyMaker.cpp
# End Source File
# Begin Source File

SOURCE=.\PolyMaker.h
# End Source File
# Begin Source File

SOURCE=.\StructDef.h
# End Source File
# End Group
# Begin Group "GGeom"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\C3DPoint.cpp
# End Source File
# Begin Source File

SOURCE=.\C3DPoint.h
# End Source File
# Begin Source File

SOURCE=.\GClipPolybyVatti.cpp
# End Source File
# Begin Source File

SOURCE=.\GClipPolybyVatti.h
# End Source File
# Begin Source File

SOURCE=.\GGeometryEngine.cpp
# End Source File
# Begin Source File

SOURCE=.\GGeometryEngine.h
# End Source File
# Begin Source File

SOURCE=.\GLine2D.cpp
# End Source File
# Begin Source File

SOURCE=.\GLine2D.h
# End Source File
# Begin Source File

SOURCE=.\GLine3D.cpp
# End Source File
# Begin Source File

SOURCE=.\GLine3D.h
# End Source File
# Begin Source File

SOURCE=.\GMathUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\GMathUtil.h
# End Source File
# Begin Source File

SOURCE=.\GPoint2D.cpp
# End Source File
# Begin Source File

SOURCE=.\GPoint2D.h
# End Source File
# Begin Source File

SOURCE=.\GPoint3D.cpp
# End Source File
# Begin Source File

SOURCE=.\GPoint3D.h
# End Source File
# Begin Source File

SOURCE=.\GRect.cpp
# End Source File
# Begin Source File

SOURCE=.\GRect.h
# End Source File
# Begin Source File

SOURCE=.\GRenderAll.h
# End Source File
# Begin Source File

SOURCE=.\GVector.cpp
# End Source File
# Begin Source File

SOURCE=.\GVector.h
# End Source File
# Begin Source File

SOURCE=.\LGLtype.h
# End Source File
# End Group
# Begin Source File

SOURCE=.\ReadMe.txt
# End Source File
# End Target
# End Project
