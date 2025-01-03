#ifndef PATHSIM_PARAMS_HPP
#define PATHSIM_PARAMS_HPP

#include <string>
#include <vector>

namespace PathSim {

struct PathParams
{
	// Delay of the path in ms, 0 to 30ms.
	double  delay       { 0. };
	// Raileigh spread
	double	spread	 	{ 0. };
	// Frequency offset in Hz
	double	offset 		{ 0. };
};

struct NoiseParams
{
	// Is white gaussian noise generator enabled?
	bool 	has_awgn	{ false };
	// Signal to noise (SNR) value, dB
    double	snr		 	{ 10. };
};

struct PathSimParams
{
	// Description of the profile
	std::string 			title;
	// Name of the profile on command line
	std::string 			cmdline_param;

    std::vector<PathParams> paths;
    NoiseParams             noise;
};

extern const std::vector<PathSimParams>& default_params();

struct State {
	double  m_SigRMS 			{ 0. };
	double  m_NoiseRMS 			{ 0. };
};

} // namespace PathSim

#endif // PATHSIM_PARAMS_HPP
