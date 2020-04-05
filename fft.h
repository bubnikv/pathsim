// fft.h: interface for the Cfft class.
//
//  This is a slightly modified version of Takuya OOURA's
//     original radix 4 FFT package.
//Copyright(C) 1996-1998 Takuya OOURA
//    (email: ooura@mmm.t.u-tokyo.ac.jp).
//////////////////////////////////////////////////////////////////////

#ifndef PATHSIM_FFT_HPP
#define PATHSIM_FFT_HPP

#include <math.h>
#include "cmplx.h"

namespace PathSim {

static constexpr int FFT_SIZE = 2048;

class Cfft
{
public:
	void CalcFFT(double * InBuf);
	void ResetFFT();
	void SetFFTParams(  int ave, double gain, int type);
	bool GetFFTData( int start, int stop, long* OutBuf);
	Cfft();
	virtual ~Cfft();
private:
	bool m_Overload;
	int m_LastAve;
	double m_Gain;
	double m_Clip;
	bool m_LogMode;
	bool m_LastLogMode;
	int m_AveSize;
	double * SinCosTbl;
	double * WindowTbl;
	double* pFFTAveBuf;
	double* pFFTInBuf;
	int * WorkArea;
	void makewt(int nw, int *ip, double *w);
	void makect(int nc, int *ip, double *c);
	void bitrv2(int n, int *ip, double *a);
	void cftfsub(int n, double *a, double *w);
	void rftfsub(int n, double *a, int nc, double *c);
    void cft1st(int n, double *a, double *w);
    void cftmdl(int n, int l, double *a, double *w);
};

} // namespace PathSim

#endif // PATHSIM_FFT_HPP
