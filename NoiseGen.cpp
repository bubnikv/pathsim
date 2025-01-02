#include "NoiseGen.h"

#include <math.h>

#ifdef WIN32
#include <stdlib.h>
#endif // WIN32

namespace PathSim {

// Obtained experimentally to compensate for BP filter.
static constexpr double K_ENBW = 1.10;

void NoiseGen::InitNoiseGen()
{
	for (int i = 0; i < HILBPFIR_LENGTH; ++ i)
		m_queue[i] = 0.0;		// fill delay buffer with zero
	m_FirState = HILBPFIR_LENGTH - 1;
}

//////////////////////////////////////////////////////////////////
// Adds bufsize gaussian random doubles with 0 mean and
// RMSlevel = RMS = std to the specified buffer, pIn
//////////////////////////////////////////////////////////////////
void NoiseGen::AddBWLimitedNoise(int bufsize, double *pIn, double siggain,double RMSlevel)
{
	int i = 0;
	double rad;
	double r;
	double u1;
	double u2;
	const double* Kptr;
	double* Firptr;

	RMSlevel = RMSlevel*K_ENBW;	//ENBW gain compensation(measured experimentally)
	while( i<bufsize )
	{
// Generate two uniform random numbers between -1 and +1
// that are inside the unit circle
		do {
			u1 = 1.0 - 2.0 * (double)rand()/(double)RAND_MAX ;
			u2 = 1.0 - 2.0 * (double)rand()/(double)RAND_MAX ;
			r = u1*u1 + u2*u2;
		} while(r >= 1.0 || r == 0.0);
		rad = sqrt(-2.0*log(r)/r);

// Performance test stuff
//HistogramSamp(  u1*rad, -6.0, 6.0, 10000000 );
//HistogramSamp(  1.0 - 2.0 * (double)rand()/(double)RAND_MAX, -1.5, 1.5, 10000000 );
//pIn[i] = 0;
//pIn[i+1] = 0;

#if 1	//if want to use 3 Khz BP filtered noise set to 1 else is flat
// 3KHz BP filter the Gaussian noise(use one of the Hilbert 3Khz coefficient tables)
		m_queue[m_FirState] = (RMSlevel*u1*rad);	//place in circular Queue
		Firptr = m_queue;
		u1 = 0.0;
		Kptr = IHilbertBPFirCoef+HILBPFIR_LENGTH-m_FirState;
		for(int j=0; j<	HILBPFIR_LENGTH;j++)
			u1 += ( (*Firptr++)*(*Kptr++) );
		if( --m_FirState < 0)
			m_FirState = HILBPFIR_LENGTH-1;
//  Add BP filtered noise to signal
		pIn[i] = siggain*pIn[i++] + u1;

		m_queue[m_FirState] = (RMSlevel*u2*rad);	//place in circular Queue
		Firptr = m_queue;
		u2 = 0.0;
		Kptr = IHilbertBPFirCoef+HILBPFIR_LENGTH-m_FirState;
		for(int j=0; j<	HILBPFIR_LENGTH;j++)
			u2 += ( (*Firptr++)*(*Kptr++) );
		if( --m_FirState < 0)
			m_FirState = HILBPFIR_LENGTH-1;
//  Add BP filtered noise to signal
		pIn[i] = siggain*pIn[i++] + u2;
#else
		pIn[i] = siggain*pIn[i++] + RMSlevel*u1*rad;
		pIn[i] = siggain*pIn[i++] + RMSlevel*u2*rad;
#endif
	}
//gDebug1 = CalcTotalRMS(pIn, bufsize);
}

} // namespace PathSim
