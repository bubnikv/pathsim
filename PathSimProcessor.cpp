#include "PathSimProcessor.h"

namespace PathSim {

// RMS amplitude out of 32768
// pick so that worst case Gaussian noise
// plus signals will not overflow soundcard
static constexpr double RMS_MAXAMPLITUDE = 4000.0;

static constexpr double RMSAVE = 20.;

PathSimProcessor::~PathSimProcessor()
{
#ifdef PATHSIM_TESTMODE
    EndTest();
#endif
}

void PathSimProcessor::init(const PathSimParams &params)
{
    m_params = params;

    m_pDelay0Buf.assign(BUF_SIZE, cmplx());
    m_pDelay1Buf.assign(BUF_SIZE, cmplx());
    m_pDelay2Buf.assign(BUF_SIZE, cmplx());
    int numpaths = int(m_params.has_path0) + int(m_params.has_path1) + int(m_params.has_path2);
    m_path0.InitPath(m_params.spread0, m_params.offset0, BUF_SIZE, numpaths, m_params.has_path0);
    m_path1.InitPath(m_params.spread1, m_params.offset1, BUF_SIZE, numpaths, m_params.has_path1);
    m_path2.InitPath(m_params.spread2, m_params.offset2, BUF_SIZE, numpaths, m_params.has_path2);
    m_direct_path = numpaths == 0;

    m_noise_gen.InitNoiseGen();
    if (! m_direct_path) { //if not direct path
        // delay/Hilbert filter
        m_delay.SetDelays(m_params.delay1, m_params.delay2);
    }

    m_SigRMS = RMS_MAXAMPLITUDE;
#ifdef PATHSIM_TESTMODE
    InitPerformance();
#endif
}

void PathSimProcessor::process_buffer(double *buffer)
{
#ifdef PATHSIM_TESTMODE
    StartPerformance();
#endif

//calculate sum of squares for RMS calculations
    for (int i = 0; i<BUF_SIZE; i++)
    {
//		SweepGen(  &(buffer[i]),  8000, 10, 4000, 5.0 );
        m_SSum = m_SSum + ( buffer[i]*buffer[i] );
    }
//Simple IIR LP filter the rms averages
    m_SigRMS = (1.0/RMSAVE)*sqrt(m_SSum/BUF_SIZE) + (1.0-1.0/RMSAVE)*m_SigRMS;
    m_SSum = 0.0;
    if (! m_params.has_awgn)
    {
        m_SignalGain = 1.0;
        m_NoiseRMS = 0.0;
    }
    m_state.m_SigRMS = m_SigRMS;
//
    if (! m_direct_path)	//if not direct path
    {
//Bandpass filter into I and Q and get delayed versions of the input data
        m_delay.CalcBPFilter(buffer, m_pDelay0Buf.data());

//double tmp;
//	for( i = 0; i<BUF_SIZE; i++)
//	{
//		tmp = m_pDelay0Buf[i].y;
//		CalcSweepRMS(tmp, 10000);
//		CalcCpxSweepPhz(m_pDelay0Buf[i],10000);
//	}

        m_delay.CreateDelays(m_pDelay0Buf.data(), m_pDelay1Buf.data(), m_pDelay2Buf.data());
//gDebug2 = CalcCpxRunAveRMS(m_pDelay2Buf, BUF_SIZE);

//Calculate each path.
        m_path0.CalcPath(m_pDelay0Buf.data(), m_pDelay0Buf.data());
        m_path1.CalcPath(m_pDelay1Buf.data(), m_pDelay1Buf.data());
        m_path2.CalcPath(m_pDelay2Buf.data(), m_pDelay2Buf.data());
// Sum and Copy just the real part back into the real buffer for output
        for( int i = 0; i<BUF_SIZE; i++)
            buffer[i] = m_pDelay0Buf[i].r + m_pDelay1Buf[i].r + m_pDelay2Buf[i].r;
    }
    if(m_params.has_awgn)	//if AWGN is used, figure out gains for SNR
    {
//		const std::lock_guard<std::mutex> lock(m_mutex);
        if(m_SNR>=1.0)
        {
            m_SignalGain = RMS_MAXAMPLITUDE/m_SigRMS;
            m_NoiseRMS = (m_SignalGain*m_SigRMS)/m_SNR;
            m_noise_gen.AddBWLimitedNoise(BUF_SIZE, buffer, m_SignalGain, m_NoiseRMS);
        }
        else
        {
            m_SignalGain = (RMS_MAXAMPLITUDE*m_SNR)/m_SigRMS;
            m_NoiseRMS = RMS_MAXAMPLITUDE;
            m_noise_gen.AddBWLimitedNoise(BUF_SIZE, buffer, m_SignalGain, RMS_MAXAMPLITUDE);
        }
    }
    m_state.m_NoiseRMS = m_NoiseRMS;

#ifdef PATHSIM_TESTMODE
    StopPerformance();
#endif
}

} // namespace PathSim
