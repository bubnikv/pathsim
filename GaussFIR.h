#ifndef PATHSIM_GAUSS_FIR_HPP
#define PATHSIM_GAUSS_FIR_HPP

#include <math.h>
#include <vector>
#include "cmplx.h"

namespace PathSim {

class GaussFIR  
{
public:
	void 	Init(double Fs, double F2sig);
	cmplx 	CalcFilter(const cmplx in);

private:
	std::vector<double> m_coef;
	int 				m_FIRlen;
	std::vector<cmplx> 	m_data;
	int 				m_data_ptr { 0 };

	static constexpr double SQRT2PI = 2.506628274631000502415765284811;
	static constexpr double PI2     = 6.283185307179586476925286766559;
	static constexpr double SQRT2 	= 1.4142135623730950488016887242097;
	// constant determines accuracy of Gaussian LPFIR.
	// larger number makes FIR longer but more accurate
	// This value is good to about -50 dB
	static constexpr double K_GAUSSIAN = 1.4;

	// Gaussian (Normal) distribution function.
	static double dnorm(double x, double mu, double sigma) 
		{ return( 1.0/(SQRT2PI*sigma) )*exp( (-1.0/(2.0*sigma*sigma))*(x-mu)*(x-mu) ); }
};

} // namespace PathSim

#endif // PATHSIM_GAUSS_FIR_HPP
