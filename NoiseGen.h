#ifndef PATHSIM_NOISEGEN_HPP
#define PATHSIM_NOISEGEN_HPP

#include "FilterTables.h"

namespace PathSim {

class NoiseGen  
{
public:
	void InitNoiseGen();
	void AddBWLimitedNoise(int bufsize,double* pIn, double siggain, double RMSlevel);

private:
	double 	m_queue[HILBPFIR_LENGTH];
	int 	m_FirState;
};

} // namespace PathSim

#endif // PATHSIM_NOISEGEN_HPP
