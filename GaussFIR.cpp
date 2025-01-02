#include "GaussFIR.h"

namespace PathSim {

// Calculate length and FIR coefficients for a Gaussian shaped low pass filter.
void GaussFIR::Init(double Fs, double F2sig)
{
	double sigma = (Fs * SQRT2) / (PI2 * F2sig);
	m_FIRlen = int(0.5 + K_GAUSSIAN * Fs / F2sig);
	if (! (m_FIRlen & 1))
		// make FIR length ODD
		++ m_FIRlen;
	// Allocate buffer and Coefficient memory based on calculated FIR length.
    m_coef.assign(m_FIRlen * 2, 0.);
    m_data.assign(m_FIRlen, cmplx());
	// generate the scaled Gaussian shaped impulse response	to create a 0 dB
	//   passband LP filter with a 2 Sigma frequency bandwidth.
	int indx = - (m_FIRlen - 1) / 2;
	for (int i = 0; i < m_FIRlen; ++ i, ++ indx) {
		m_coef[i] = (1.0 / (SQRT2PI * sigma)) * dnorm(indx, 0.0, sigma) / dnorm(0.0, 0.0, sigma);
		// make duplicate for flat FIR
		m_coef[i + m_FIRlen] = m_coef[i];
	}
	// used for flat FIR implementation
	m_data_ptr = m_FIRlen - 1;
}

// Calculate complex Gaussian FIR filter iteration for one sample.
cmplx GaussFIR::CalcFilter(const cmplx in)
{
	m_data[m_data_ptr] = in;
	cmplx  acc;
    const double *coeff = m_coef.data() + m_FIRlen - m_data_ptr;
	for (int i = 0; i < m_FIRlen; ++ i, ++ coeff) {
        acc.r += m_data[i].r * *coeff;
        acc.i += m_data[i].i * *coeff;
	}
	if (-- m_data_ptr < 0)
		m_data_ptr += m_FIRlen;
	return acc;
}

} // namespace PathSim
