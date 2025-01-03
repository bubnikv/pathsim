#include "PathSimParams.h"

namespace PathSim {

std::vector<PathSimParams> s_default_params {
	// Description,                          cmdline,
	{ "Direct Path", 						"direct" },
																// noise enabled, 10dB SNR
	{ "AWGN S/N",							"awgn10", 			{}, {true,  10.} },
																// { { DELAY, SPREAD, OFFSET } ... }
	{ "CCIR 520-2 (Doppler Fading)",		"ccir-doppler", 	{ {0.,   0.2,  10.0}, {0.5,   0.2, 5.} } },
	{ "CCIR 520-2 (Flutter Fading)",		"ccir-flutter", 	{ {0.,  10.0,   0.0}, {0.5,  10. , 0.} } },
	{ "CCIR 520-2 (Good Conditions)",		"ccir-good", 		{ {0.,   0.1,   0.0}, {0.5,   0.1, 0.} } },
	{ "CCIR 520-2 (Moderate Conditions)", 	"ccir-moderate", 	{ {0.,   0.5,   0.0}, {1. ,   0.5, 0.} } },
	{ "CCIR 520-2 (Poor Conditions)",		"ccir-poor", 		{ {0.,   0.5,   0.0}, {2. ,   0.5, 0.} } },

    { "Low-Latitude Disturbed", 			"lowlat-disturbed", { {0.,  10.0,   0.0}, {6. ,  10. , 0.} } },
    { "Low-Latitude Moderate",				"lowlat-moderate", 	{ {0.,   1.5,   0.0}, {2. ,   1.5, 0.} } },
    { "Low-Latitude Quiet",					"lowlat-queit", 	{ {0.,   0.5,   0.0}, {0.5,   0.5, 0.} } },

	{ "High-Latitude Disturbed",			"highlat-disturbed",{ {0.,  30.0,   0.0}, {7. ,  30. , 0.} } },
	{ "High-Latitude Moderate", 			"highlat-moderate", { {0.,  10.0,   0.0}, {3. ,  10. , 0.} } },
	{ "High-Latitude Quiet", 				"highlat-quiet", 	{ {0.,   0.5,   0.0}, {1. ,   0.5, 0.} } },

	{ "Mid-Latitude Disturbed", 			"midlat-disturbed", { {0.,   1.0,   0.0}, {2. ,   1. , 0.} } },
	{ "Mid-Latitude Disturbed NVIS", 		"midlat-dist-nvis", { {0.,   1.0,   0.0}, {7. ,   1. , 0.} } },
	{ "Mid-Latitude Moderate", 				"midlat-moderate", 	{ {0.,   1.0,   0.5} } },
	{ "Mid-Latitude Quiet", 				"midlat-quiet", 	{ {0.,   0.5,   0.0}, {0.5,   0.5, 0.} } },

    { "Frequency Shifter", 					"freq-shifter", 	{ {0.,   0.0, 500.0} } },
};

const std::vector<PathSimParams>& default_params() 
{
	return s_default_params;
}

} // namespace PathSim 
