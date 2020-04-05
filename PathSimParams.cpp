#include "PathSimParams.h"

namespace PathSim {

PathSimParams defaut_params[default_params_cnt] = {
	// Title,                               Cmdline,            AWGN,  S/N, PATH0,SPREAD,OFFSET,PATH1,DELAY,SPREAD,OFFSET
	{ "Direct Path", 						"direct" },
	{ "AWGN S/N",							"awgn10", 			true,  10. },

	{ "CCIR 520-2 (Doppler Fading)",		"ccir-doppler", 	false,  0.,  true,   0.2,  10.0, true,  0.5,   0.2, 5.  },
	{ "CCIR 520-2 (Flutter Fading)",		"ccir-flutter", 	false,  0.,  true,  10.0,   0.0, true,  0.5,  10. , 0.  },
	{ "CCIR 520-2 (Good Conditions)",		"ccir-good", 		false,  0.,  true,   0.1,   0.0, true,  0.5,   0.1, 0.  },
	{ "CCIR 520-2 (Moderate Conditions)", 	"ccir-moderate", 	false,  0.,  true,   0.5,   0.0, true,  1. ,   0.5, 0.  },
	{ "CCIR 520-2 (Poor Conditions)",		"ccir-poor", 		false,  0.,  true,   0.5,   0.0, true,  2. ,   0.5, 0.  },

    { "Low-Latitude Disturbed", 			"lowlat-disturbed", false,  0.,  true,  10.0,   0.0, true,  6. ,  10. , 0.  },
    { "Low-Latitude Moderate",				"lowlat-moderate", 	false,  0.,  true,   1.5,   0.0, true,  2. ,   1.5, 0.  },
    { "Low-Latitude Quiet",					"lowlat-queit", 	false,  0.,  true,   0.5,   0.0, true,  0.5,   0.5, 0.  },

	{ "High-Latitude Disturbed",			"highlat-disturbed",false,  0.,  true,  30.0,   0.0, true,  7. ,  30. , 0.  },
	{ "High-Latitude Moderate", 			"highlat-moderate", false,  0.,  true,  10.0,   0.0, true,  3. ,  10. , 0.  },
	{ "High-Latitude Quiet", 				"highlat-quiet", 	false,  0.,  true,   0.5,   0.0, true,  1. ,   0.5, 0.  },

	{ "Mid-Latitude Disturbed", 			"midlat-disturbed", false,  0.,  true,   1.0,   0.0, true,  2. ,   1. , 0.  },
	{ "Mid-Latitude Disturbed NVIS", 		"midlat-dist-nvis", false,  0.,  true,   1.0,   0.0, true,  7. ,   1. , 0.  },
	{ "Mid-Latitude Moderate", 				"midlat-moderate", 	false,  0.,  true,   1.0,   0.5, false, 1. ,   1. , 0.5 },
	{ "Mid-Latitude Quiet", 				"midlat-quiet", 	false,  0.,  true,   0.5,   0.0, true,  0.5,   0.5, 0.  },

    { "Frequency Shifter", 					"freq-shifter", 	false,  0.,  true,   0.0, 500.0 },
};

} // namespace PathSim
