#ifndef PATHSIM_PROCESSOR_HPP
#define PATHSIM_PROCESSOR_HPP

#include "Path.h"
#include "PathSimParams.h"
#include "NoiseGen.h"
#include "Delay.h"

#include <math.h>

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
    double CalcRMS( double* buf);
    double CalcRMS2( double* buf);

    double                  m_SigRMS 		{ 0. };
    double                  m_SignalGain 	{ 0. };
    double                  m_NoiseRMS  	{ 0. };
    double                  m_SNR 			{ 0. };
    double                  m_SSum			{ 0. };

    std::vector<cmplx>      m_pDelay0Buf;
    std::vector<cmplx>      m_pDelay1Buf;
    std::vector<cmplx>      m_pDelay2Buf;
    Delay                   m_delay;
    Path                    m_path0;
    Path                    m_path1;
    Path                    m_path2;
    NoiseGen                m_noise_gen;

    PathSimParams           m_params;
    // Direct path is active if m_params.has_path0/1/2 are all false.
    bool                    m_direct_path  { false };

    // For visualization?
    State                   m_state;
};

} // namespace PathSim

#endif // PATHSIM_PROCESSOR_HPP
