#ifndef PATHSIM_DELAY_HPP
#define PATHSIM_DELAY_HPP

#include <math.h>
#include <vector>

#include "cmplx.h"
#include "FilterTables.h"

namespace PathSim {

class Delay
{
public:
	void CreateDelays( cmplx* inbuf, cmplx* t1buf, cmplx* t2buf );
	void SetDelays( double t1, double t2);
	void CalcBPFilter(double* pIn, cmplx* pOut);

private:
    // 50mSecs max delay
    static constexpr int MAXDELAY   = 50 * 8;
    static constexpr int BLOCKSIZE  = 2048;
    static constexpr int BUFSIZE    = BLOCKSIZE + MAXDELAY;

    cmplx 	m_delay_line[BUFSIZE];
    cmplx 	m_queue[HILBPFIR_LENGTH];
    int 	m_FirState = 0;
    int 	m_T1ptr = 0;
    int 	m_T2ptr = 0;
    int 	m_Inptr = 0;
};

} // namespace PathSim

#endif // PATHSIM_DELAY_HPP
