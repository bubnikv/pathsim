// Gaussian low pass filter

#ifndef PATHSIM_GAUSS_FIR_HPP
#define PATHSIM_GAUSS_FIR_HPP

#include <math.h>
#include <vector>
#include "cmplx.h"

namespace PathSim {

class GaussFIR
{
public:
	void 	init(double Fs, double F2sig);
	cmplx 	apply(const cmplx in);

private:
	// Gaussian filter coefficients
	std::vector<double> m_coef;
	// Circular queue of filtered samples.
	std::vector<cmplx> 	m_data;
	int 				m_data_ptr { 0 };
};

} // namespace PathSim

#endif // PATHSIM_GAUSS_FIR_HPP
