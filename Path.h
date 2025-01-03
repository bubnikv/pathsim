#ifndef PATHSIM_PATH_HPP
#define PATHSIM_PATH_HPP

#include <math.h>

#include "cmplx.h"
#include "FilterTables.h"
#include "GaussFIR.h"

namespace PathSim {

// Rayleigh distribution
class Rayleigh {
public:
	enum class SampleRate {
		// Indices correspond to Path::m_upsamplers!
	    Rate_12_point_8_Hz = 3, // Used for 0.1 > Spread >= 0.4
	    Rate_64_Hz = 2,			// Used for 0.4 > Spread >= 2.0
	    Rate_320_Hz = 1,       	// Used for 2.0 > Spread >= 10.0
	    Rate_None = 0, 			// Used for super low Spread < 0.1
	};

	void  		init(double spread, double gain_coeff);
	cmplx 		sample();
	SampleRate  sample_rate() const { return m_sample_rate; }

private:
	double 	 	m_spread;
	double 	 	m_gain;
	SampleRate 	m_sample_rate { SampleRate::Rate_12_point_8_Hz };

	// Gaussian FIR low pass filter
	GaussFIR 	m_lpfir;
};

class Path
{
public:
	Path() {}

	void init_path(double spread, double offset, int blocksize, int numpaths);
	void calc_path(const cmplx* pIn, cmplx* pOut);

private:
	// Polyphase upsampling low pass FIR filter.
	class Upsampler
	{
	public:
		void  init(int rate);
		void  insert_sample(cmplx v);
		cmplx upsample();

		int   insert_at(int index) const { return index % m_rate == 0; }

	private:
		// Samples to be upsampled and low pass filtered by a polyphase filter.
		cmplx m_queue[INTP_QUE_SIZE];
		// Pointer has a span of INTP_VALUE x INTP_QUE_SIZE.
		int   m_ptr;
		// Upsampling rate
		int   m_rate;
	};

	int 		m_block_size;
	int 		m_index { 0 };
	double 		m_offset_frequency;
	double 		m_phase_acc;

	Rayleigh 	m_rayleigh;
	Upsampler 	m_upsamplers[4];
};

} // namespace PathSim

#endif // PATHSIM_PATH_HPP
