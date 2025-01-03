#include "Path.h"

#include <assert.h>

#define _USE_MATH_DEFINES
#include <math.h>

static constexpr const double OFFSET_FREQ_CONST = 2. * M_PI / 8000.;
static constexpr const double KGNB = 0.62665707; // equivalent Noise BW of Gaussian shaped filter

namespace PathSim {

void Rayleigh::init(double spread, double gain_coeff)
{
    assert(spread >= 0 && spread <= 30.0);
    if (spread < 0.)
        spread = 0.;
    else if (spread > 30.0)
        spread = 30.;
    m_spread = spread;

    if (spread < 0.1) {
        // here if spread<.1 so will not use any spread just offset
        m_sample_rate = SampleRate::Rate_None;
        m_gain        = gain_coeff;
    } else {
        double rate;
        if (spread > 2.0) {
            m_sample_rate = SampleRate::Rate_320_Hz;
            rate = 320.0;
        } else if (spread > 0.4) {
            m_sample_rate = SampleRate::Rate_64_Hz;
            rate = 64.;
        } else if (spread >= 0.1) {
            m_sample_rate = SampleRate::Rate_12_point_8_Hz;
            rate = 12.8;
        }
        m_lpfir.init(rate, spread);
        m_gain = gain_coeff * sqrt(rate / (4.0 * spread * KGNB));

        // preload m_lpfir
        for (int i = 0; i < 250; ++ i)
            this->sample();
    }
}

// Create the complex Rayleigh distributed samples by
// creating two Gaussian random distributed numbers for the I and Q
// terms and then passing them through a Gaussian shaped low pass FIR filter.
// The 2 Sigma bandwidth of the LP filter determines the amount of spread.
cmplx Rayleigh::sample()
{
    cmplx out;
    if (m_spread >= 0.1) {
        // Generate two uniform random numbers between -1 and +1 that are inside the unit circle.
        cmplx  v;
        double r2;
        do {
            // Generate vector inside a (-1, 1) x (-1, 1) 
            auto rsample = [](){ return 1. - 2. * double(rand())/double(RAND_MAX); };
            v.set(rsample(), rsample());
            r2 = v.l2();
            // Reject samples outside of the circle.
        } while (r2 >= 1. || r2 == 0.);
        // Convert the uniformly sampled unit radius circle distribution into
        // a complex normal distribution, then low pass filter with a Gaussian shaped FIR filter.
        // r2 has a uniform distribution in (0, 1)
        // and v / sqrt(r2) is a unit vector,
        // thus sqrt(- 2. * log(r2)) has a Rayleigh distribution
        // and v * sqrt(- 2. * log(r2) / r2) has a complex normal distribution.
        out = m_lpfir.apply(v * (m_gain * sqrt(- 2. * log(r2) / r2)));
    } else
        // Not using any spread.
        out.set(m_gain, 0);

    //gDebug1 = CalcCpxRMS( out, 288000);
    //CalcCpxSweepRMS( out, 500);
    return out;
}

void Path::Upsampler::init(int rate)
{
    memset(m_queue, 0, sizeof(m_queue));
    m_ptr = INTP_QUE_SIZE - 1;
    m_rate = rate;
}

void Path::Upsampler::insert_sample(cmplx v)
{
    m_queue[m_ptr / INTP_VALUE] = v;
}

cmplx Path::Upsampler::upsample()
{
    const cmplx*    firptr = m_queue;
    const double*   kptr   = X5IntrpFIRCoef + INTP_FIR_SIZE - m_ptr;
    cmplx           acc{ 0., 0. };
    for (int j = 0; j < INTP_QUE_SIZE; ++ j, kptr += INTP_VALUE, ++ firptr)
        acc += (*firptr) * (*kptr);
    if (-- m_ptr < 0)
        m_ptr = INTP_FIR_SIZE - 1;
    return acc;
}

void Path::init_path(double spread, double offset, int blocksize, int numpaths)
{
    m_block_size        = blocksize;
    m_offset_frequency  = offset;
    m_index             = 0;
    m_phase_acc         = 0.;

    int rate = INTP_VALUE;
    for (int i = 0; i < 4; ++ i, rate *= INTP_VALUE)
       m_upsamplers[i].init(rate);

    m_rayleigh.init(spread, 1. / sqrt(double(numpaths)));
}

// Performs a path calculation on pIn and puts it in pOut
//
//  Two Low Pass filtered Gaussian random numbers are created at
//	12.8, 64 Hz, or 320 Hz rate.  These form the input to a complex
//	interpolation filter that bumps the sample rate up to 8000Hz.
//
//	Two, three, or four stages of X5 upsampling/interpolation are used.
//	The complex noise is then multiplied by the input I/Q signal
//	to produce the spreading/fading simulation.
//
//  Finally a complex NCO is multiplied by the signal to produce a
//	Frequency offset.
void Path::calc_path(const cmplx *pIn, cmplx *pOut)
{
    for (int i = 0; i < m_block_size; ++ i) {
        {
            int j = int(m_rayleigh.sample_rate());
            if (m_upsamplers[j].insert_at(m_index))
                m_upsamplers[j].insert_sample(m_rayleigh.sample());
            for (-- j; j >= 0; -- j)
                if (m_upsamplers[j].insert_at(m_index))
                    m_upsamplers[j].insert_sample(m_upsamplers[j + 1].upsample());
        }
        cmplx fading = m_upsamplers[0].upsample();
//CalcCpxSweepRMS( fading, 8000);
        pOut[i] =
            // Doppler
            cmplx{cos(m_phase_acc), sin(m_phase_acc)} *
            // Fading
            (fading * pIn[i]);
        // keep radian counter bounded
        m_phase_acc = fmod(m_phase_acc + OFFSET_FREQ_CONST * m_offset_frequency, 2. * M_PI);
        if (++ m_index >= INTP_VALUE*INTP_VALUE*INTP_VALUE*INTP_VALUE*m_block_size)
            m_index = 0;
    }
}

} // namespace PathSim
