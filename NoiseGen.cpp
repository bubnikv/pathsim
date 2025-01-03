#include "NoiseGen.h"

#include <math.h>
#include <string.h>

#ifdef WIN32
#include <stdlib.h>
#endif // WIN32

#include "cmplx.h"

namespace PathSim {

// Obtained experimentally to compensate for BP filter.
static constexpr double K_ENBW = 1.10;

void NoiseGen::init(bool band_limited)
{
	m_band_limited = band_limited;
	memset(m_queue, 0, sizeof(m_queue));
	m_queue_pos = HILBPFIR_LENGTH - 1;
}

// Adds bufsize gaussian random doubles with 0 mean and
// RMSlevel = RMS = std to the specified buffer, pInOut
void NoiseGen::add_band_limited_noise(int bufsize, double *pInOut, double siggain, double RMSlevel)
{
	RMSlevel *= K_ENBW;	// ENBW gain compensation(measured experimentally)
	for (int i = 0; i < bufsize; ) {
		// Generate two normally distributed samples as real and imaginatory components
		// of a normal complex distribution.
		double noise[2];
		{
			double r2;
			cmplx v;
			do {
	            auto rsample = [](){ return 1. - 2. * double(rand())/double(RAND_MAX); };
	            v.set(rsample(), rsample());
				r2 = v.l2();
			} while (r2 >= 1. || r2 == 0.);
			v *= RMSlevel * sqrt(-2. * log(r2) / r2);
			noise[0] = v.r;
			noise[1] = v.i;
		}

		for (int j = 0; j < 2; ++ j) {
			double acc = 0.;
			if (m_band_limited) {
				// 3KHz BP filter the Gaussian noise(use one of the Hilbert 3Khz coefficient tables)
				m_queue[m_queue_pos] = noise[j];
				const double* Firptr = m_queue;
				const double* Kptr = IHilbertBPFirCoef+HILBPFIR_LENGTH-m_queue_pos;
				for (int k = 0; k < HILBPFIR_LENGTH; ++ k)
					acc += (*Firptr++) * (*Kptr++);
				if (-- m_queue_pos < 0)
					m_queue_pos = HILBPFIR_LENGTH - 1;
			} else
				acc = noise[j];
			//  Add BP filtered noise to signal
			pInOut[i] = siggain * pInOut[i ++] + acc;
		}
	}
}

} // namespace PathSim
