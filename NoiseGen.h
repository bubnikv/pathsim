#ifndef PATHSIM_NOISEGEN_HPP
#define PATHSIM_NOISEGEN_HPP

#include "FilterTables.h"

namespace PathSim {

class NoiseGen  
{
public:
	void init(bool band_limited);
	void add_band_limited_noise(int bufsize, double* pInOut, double siggain, double RMSlevel);

private:
	double 	m_queue[HILBPFIR_LENGTH];
	int 	m_queue_pos;
	bool    m_band_limited;
};

} // namespace PathSim

#endif // PATHSIM_NOISEGEN_HPP
