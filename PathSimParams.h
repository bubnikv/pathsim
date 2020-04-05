#ifndef PATHSIM_PARAMS_HPP
#define PATHSIM_PARAMS_HPP

#include <string>

namespace PathSim {

struct PathSimParams
{
	std::string title;
	std::string cmdline_param;

	// Is white gaussian noise generator enabled?
	bool 	has_awgn	{ false };
	// Signal to noise (SNR) value
    double	snr		 	{ 10. };

	bool 	has_path0 	{ false };
	double	spread0	 	{ 1. };
	double	offset0 	{ 0. };

	bool 	has_path1 	{ false };
	double	delay1  	{ 0. };
	double	spread1 	{ 1. };
	double	offset1 	{ 0. };

	bool 	has_path2 	{ false };
	double	delay2  	{ 0. };
	double	spread2 	{ 1. };
	double	offset2 	{ 0. };
};

static constexpr int default_params_cnt = 18;
extern PathSimParams defaut_params[default_params_cnt];

struct State {
	double  m_SigRMS 			{ 0. };
	double  m_NoiseRMS 			{ 0. };
};

} // namespace PathSim

#endif // PATHSIM_PARAMS_HPP
