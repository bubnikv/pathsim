// Gaussian low pass filter

#include "GaussFIR.h"

namespace PathSim {

static constexpr double SQRT2PI = 2.506628274631000502415765284811;
static constexpr double PI2     = 6.283185307179586476925286766559;
static constexpr double SQRT2 	= 1.4142135623730950488016887242097;

// constant determines accuracy of Gaussian LPFIR.
// larger number makes FIR longer but more accurate
// This value is good to about -50 dB
static constexpr double K_GAUSSIAN = 1.4;

// Gaussian (Normal) distribution function.
static inline double dnorm(double x, double mu, double sigma) 
{ 
	return (1.0/(SQRT2PI*sigma))*exp((-1.0/(2.0*sigma*sigma))*(x-mu)*(x-mu));
}

// Calculate length and FIR coefficients for a Gaussian shaped low pass filter.
void GaussFIR::init(double Fs, double F2sig)
{
	double sigma   = (Fs * SQRT2) / (PI2 * F2sig);
	m_fir_len = int(0.5 + K_GAUSSIAN * Fs / F2sig);
	if (! (m_fir_len & 1))
		// make FIR length odd
		++ m_fir_len;
	// Allocate buffer and Coefficient memory based on calculated FIR length.
    m_coef.assign(m_fir_len * 2, 0.);
    m_data.assign(m_fir_len, cmplx());
	// generate the scaled Gaussian shaped impulse response	to create a 0 dB
	//   passband LP filter with a 2 Sigma frequency bandwidth.
	int index = - (m_fir_len - 1) / 2;
	double norm = (1.0 / (SQRT2PI * sigma)) / dnorm(0.0, 0.0, sigma);
	for (int i = 0; i < m_fir_len; ++ i, ++ index) {
		m_coef[i] = norm * dnorm(index, 0.0, sigma);
		// make duplicate for flat FIR
		m_coef[i + m_fir_len] = m_coef[i];
	}
	// used for flat FIR implementation
	m_data_ptr = m_fir_len - 1;
}

// Calculate complex Gaussian FIR filter iteration for one sample.
cmplx GaussFIR::apply(const cmplx in)
{
	m_data[m_data_ptr] = in;
	cmplx acc{ 0., 0. };
    const double *coeff = m_coef.data() + m_fir_len - m_data_ptr;
	for (int i = 0; i < m_fir_len; ++ i, ++ coeff)
		// Filter each vector component by a Gaussian kernel.
        acc += m_data[i] * *coeff;
	if (-- m_data_ptr < 0)
		m_data_ptr += m_fir_len;
	return acc;
}

} // namespace PathSim
