This is the original README from the PathSim source code by Moe Wheatley AE4JY.
The source code has been reworked by Vojtech Bubnik OK1IAK.
-------------------------------------------------------------------------------





==========================================================================
++++++++++++++++++   R E A D M E . T X T   +++++++++++++++++++++++++++++++
++++++++++++++    for PathSim program source code    ++++++++++++++++++++++
==========================================================================
Ver. .10 Initial Beta release 8-12-00
Ver. .20 Initial Beta release 8-22-00
	Added Settings save menu.
Ver. 1.0 Full Release  12-2-00
	Changed Spreading filters to be true Gaussian shaped.
	Changed FIR algorithm for speed.
	Added online help.
	Changed data format of .sim files and are not compatable with 	beta release ones.

Moe Wheatley AE4JY <ae4jy@mindspring.com>  www.mindspring.com/~ae4jy
==========================================================================

System Requirements:

133MHz Pentium ( works on 486 but MUST have floating point processor)
Windows 95/98/NT ( is a 32 bit Windows application )
A soundcard that is supported by Windows.

-----------------------------------------------------------------------
For those interested in the source code for this application, it was written
Using MicroSoft Visual C++ 6.0.  MFC was used for the user interface stuff
and a single worker thread is used to generate and process all the data
streams.

Note: when executing in DEBUG mode, because the optimizer is off, execution
speed will slow significantly.
========================================================================
 VC++ 6.0      Project Source File Descriptions for : WinPSK
========================================================================
ErrorCodes.h    // some table include files
FilterTables.h

fft.cpp		//the FFT class
fft.h
           
IOCntrl.cpp	//DSP thread class
IOCntrl.h

MainFrm.cpp	//main frame of program
MainFrm.h

PlotData.cpp	// Signal data plotting class
PlotData.h

Delay.cpp	//Input BP filter and delay lines
Delay.h

GaussianFIR.cpp	//Gaussian FIR class
GaussianFIR.h

InDlg.cpp	//Input dlg
InDlg.h

OutDlg.cpp	//Output dlg
OutDlg.h

NoiseGen.cpp	//Gaussian noise gen and filter
NoiseGen.h

Sound.cpp	// soundcard interface class
Sound.h

Perform.cpp	//non Class global test bench routines
Perform.h
Wave.cpp	// Wave file interface class
Wave.h

Path.cpp	// The Path simulator class
Path.h

PathSim.cpp	// main app
PathSim.h

PathSimDoc.cpp	// document class for inter-object communications
PathSimDoc.h

PathSimView.cpp	// Form view class for most of the controls
PathSimView.h

SetupDlg.cpp // Program setup dialog
SetupDlg.h

/////////////////////////////////////////////////////////////////////////////
Other standard files:

StdAfx.h, StdAfx.cpp
    These files are used to build a precompiled header (PCH) file
    named Pathsim.pch and a precompiled types file named StdAfx.obj.

Resource.h
    This is the standard header file, which defines new resource IDs.
    Microsoft Developer Studio reads and updates this file.

/////////////////////////////////////////////////////////////////////////////
