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

void PathSimProcessor::init(const PathSimParams& params)
{
    m_params = params;

    int numpaths = int(m_params.paths.size());
    m_direct_path = numpaths == 0;

    m_paths.assign(numpaths, { {}, {BUF_SIZE, cmplx{}} });
    m_noise_gen.init(true);

    if (! m_direct_path) {
        m_hilbert.init();
        m_delay.init();
        for (const PathParams& p : params.paths) {
            int i = int(&p - params.paths.data());
            m_paths[i].path.init_path(p.spread, p.offset, BUF_SIZE, numpaths);
            if (i > 0)
                m_delay.add_delay(p.delay);
        }
    }

    m_SigRMS = RMS_MAXAMPLITUDE;
}

void PathSimProcessor::process_buffer(double *buffer)
{
    {
        // Calculate sum of squares for RMS calculations
        double acc = 0.;
        for (int i = 0; i < BUF_SIZE; ++ i)
            acc += buffer[i] * buffer[i];
        // Simple IIR LP filter the rms averages
        m_SigRMS = (1.0 / RMSAVE) * sqrt(acc / BUF_SIZE) + (1.0 - 1.0 / RMSAVE) * m_SigRMS;
    }
    if (! m_params.noise.has_awgn) {
        m_SignalGain = 1.0;
        m_NoiseRMS   = 0.0;
    }
    m_state.m_SigRMS = m_SigRMS;

    if (! m_direct_path) {
        // Bandpass filter into I and Q and get delayed versions of the input data
        static_assert(m_hilbert.BLOCKSIZE == this->BUF_SIZE, "Buffer length has to be satisfied");
        m_hilbert.filter_block(buffer, m_paths.front().buffer.data());
        static_assert(m_delay.BLOCKSIZE == this->BUF_SIZE, "Buffer length has to be satisfied");
        std::vector<std::vector<cmplx>*> buffers(m_paths.size(), nullptr);
        for (size_t i = 1; i < m_paths.size(); ++ i)
            buffers.emplace_back(&m_paths[i].buffer);
        m_delay.delay_block(m_paths.front().buffer, buffers);
        // Calculate each path.
        for (PathWithBuffer& path : m_paths)
            path.path.calc_path(path.buffer.data(), path.buffer.data());
        // Sum and Copy just the real part back into the real buffer for output
        for (int i = 0; i < BUF_SIZE; ++ i) {
            double acc = 0;
            for (PathWithBuffer& path : m_paths)
                acc += path.buffer[i].r;
            buffer[i] = acc;
        }
    }
    if (m_params.noise.has_awgn) {
        // if AWGN is used, figure out gains for SNR
        if (m_SNR >= 1.0) {
            m_SignalGain = RMS_MAXAMPLITUDE / m_SigRMS;
            m_NoiseRMS   = m_SignalGain * m_SigRMS / m_SNR;
        } else {
            m_SignalGain = RMS_MAXAMPLITUDE * m_SNR / m_SigRMS;
            m_NoiseRMS   = RMS_MAXAMPLITUDE;
        }
        m_noise_gen.add_band_limited_noise(BUF_SIZE, buffer, m_SignalGain, m_NoiseRMS);
    }
    m_state.m_NoiseRMS = m_NoiseRMS;
}

} // namespace PathSim
