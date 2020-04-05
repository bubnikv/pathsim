#ifndef PATHSIM_PATH_HPP
#define PATHSIM_PATH_HPP

#include <math.h>
#include "cmplx.h"
#include "FilterTables.h"
#include "GaussFIR.h"

namespace PathSim {

class Path  
{
public:
	Path();

	void CalcPath( cmplx* pIn, cmplx* pOut);
	void InitPath( double Spread, double Offset, int blocksize,int numpaths, bool active);

private:
	int m_IIRLength;
	cmplx MakeGaussianDelaySample();
	int m_NoiseSampRate;
	bool m_PathActive;
	int m_inc;
	int m_BlockSize;
	int m_Indx;
	double m_Offset;
	double m_Spread;
	double m_LPGain;
	double m_Timeinc;
	cmplx m_pQue0[INTP_QUE_SIZE];
	cmplx m_pQue1[INTP_QUE_SIZE];
	cmplx m_pQue2[INTP_QUE_SIZE];
	cmplx m_pQue3[INTP_QUE_SIZE];
	int m_FirState0;
	int m_FirState1;
	int m_FirState2;
	int m_FirState3;
	GaussFIR m_lpfir;
};

} // namespace PathSim

#endif // PATHSIM_PATH_HPP
