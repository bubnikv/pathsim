#ifndef PATHSIM_PROCESSOR_HPP
#define PATHSIM_PROCESSOR_HPP

#include "Path.h"
#include "PathSimParams.h"
#include "NoiseGen.h"
#include "Delay.h"

#include <math.h>
#include <utility>

namespace PathSim {

class PathSimProcessor
{
public:
    ~PathSimProcessor();

    void init(const PathSimParams &params);

	// Size of data chunks to process one at a time.
	static constexpr int BUF_SIZE = 2048;
    void process_buffer(double *buffer);

private:
    double                  m_SigRMS 		{ 0. };
    double                  m_SignalGain 	{ 0. };
    double                  m_NoiseRMS  	{ 0. };
    double                  m_SNR 			{ 0. };

    struct PathWithBuffer {
        Path                path;
        std::vector<cmplx>  buffer;
    };
    std::vector<PathWithBuffer> m_paths;
    Hilbert                 m_hilbert;
    Delay                   m_delay;
    NoiseGen                m_noise_gen;

    PathSimParams           m_params;

    // Direct path is active if m_params.has_path0/1/2 are all false.
    bool                    m_direct_path  { false };

    // Signal / Noise RMS, for diagnostics.
    State                   m_state;
};

} // namespace PathSim

#endif // PATHSIM_PROCESSOR_HPP
