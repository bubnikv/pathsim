// Delay.cpp: implementation of the CDelay class.
//   ( also performs Hilbert Real to complex I/Q 3KHz filtering )

#include "Delay.h"
#include "FilterTables.h"

namespace PathSim {

//  Set delay times
void Delay::SetDelays( double t1, double t2)	//delays in msecs
{
    for( int i = 0; i<HILBPFIR_LENGTH; ++ i) {
		m_queue[i].x = 0.0;
		m_queue[i].y = 0.0;
	}
    for( int i = 0; i<BUFSIZE; ++ i) {
		m_delay_line[i].x = 0.0;
		m_delay_line[i].y = 0.0;
	}
	m_Inptr = BUFSIZE-1;
	m_FirState = HILBPFIR_LENGTH-1;
	m_T1ptr = BUFSIZE - (int)(8.0*t1) - 1;
	m_T2ptr = BUFSIZE - (int)(8.0*t2) - 1;
}

//  Uses pointers to create variable delays
void Delay::CreateDelays( cmplx* inbuf, cmplx* t1buf, cmplx* t2buf )
{
    for (int i=0; i<BLOCKSIZE; i++)
	{
		m_delay_line[m_Inptr++] = inbuf[i]; // copy new data from inbuf into delay buffer
		if( m_Inptr >= BUFSIZE )
			m_Inptr = 0;
		t1buf[i] = m_delay_line[m_T1ptr++];
		if( m_T1ptr >= BUFSIZE )
			m_T1ptr = 0;
		t2buf[i] = m_delay_line[m_T2ptr++];
		if( m_T2ptr >= BUFSIZE )
			m_T2ptr = 0;
	}
}

// Hilbert 3KHz BP filters.  Real input and complex I/Q output
//   This FIR bandwidth limits the real input as well as creates a
//   complex I/Q output signal for the rest of the processing chain.
void Delay::CalcBPFilter(double* pIn, cmplx* pOut)
{
	cmplx acc;
	const double* IKptr;
	const double* QKptr;
	cmplx* Firptr;
	int j;
	for(int i=0; i<BLOCKSIZE; i++)
	{
		m_queue[m_FirState].x = pIn[i];	//place real values in circular Queue
		m_queue[m_FirState].y = pIn[i];
        Firptr = m_queue;
		IKptr = IHilbertBPFirCoef+HILBPFIR_LENGTH-m_FirState;
		QKptr = QHilbertBPFirCoef+HILBPFIR_LENGTH-m_FirState;	//
		acc.x = 0.0;
		acc.y = 0.0;
		for(j=0; j<	HILBPFIR_LENGTH;j++)
		{
			acc.x += ( (Firptr->x)*(*IKptr++) );
			acc.y += ( (Firptr++->y)*(*QKptr++) );
		}
		pOut[i].x = acc.x;
		pOut[i].y = acc.y;
		if( --m_FirState < 0)
			m_FirState = HILBPFIR_LENGTH-1;
	}
}

} // namespace PathSim
