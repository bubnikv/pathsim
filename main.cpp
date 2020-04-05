#include "PathSimProcessor.h"
#include "PathSimParams.h"

#include <iostream>
#include "cxxopts.h"
#include "AudioFile.h"

using namespace PathSim;

int main(int argc, char **argv)
{
	try {
	    cxxopts::Options options(argv[0], "\
pathsim - HF ionospheric propagation simulator implementing a Watterson channel model\n\
Ported by Vojtech Bubnik OK1AK to a command line tool\n\
from a Windows GUI application by Moe Wheatley, AE4JY, 2000\n\
");
		options
            .positional_help("input_file output_file")
			.show_positional_help();

		options
            .allow_unrecognised_options()
			.add_options()
            ("help", "Print help")
            ("i,input_file", "Input audio file", cxxopts::value<std::string>())
            ("o,output_file", "Output audio file", cxxopts::value<std::string>())
//                ->default_value("a.out")->implicit_value("b.def"), "BIN")
            ("positional",
                "Positional arguments: these are the arguments that are entered "
                "without an option", cxxopts::value<std::vector<std::string>>());

		{
		    auto group = options.add_options("Propagation condition");
		    for (int i = 0; i < default_params_cnt; ++ i)
		        group(defaut_params[i].cmdline_param, defaut_params[i].title);
		}

	    options.add_options("Propagation")
	        ("snr", "Signal to Noise Ratio (SNR)", cxxopts::value<double>())
	        ("spread", "Frequency spread of the 1st path [Hz]", cxxopts::value<double>())
	        ("offset", "Frequency offset of the 1st path [Hz]", cxxopts::value<double>())
	        ("delay2", "Delay of the 2nd path [Hz]", cxxopts::value<double>())
	        ("spread2", "Frequency spread of the 2nd path [Hz]", cxxopts::value<double>())
	        ("offset2", "Frequency offset of the 2nd path [Hz]", cxxopts::value<double>())
	        ("delay3", "Delay of the 2nd path [Hz]", cxxopts::value<double>())
	        ("spread3", "Frequency spread of the 3rd path [Hz]", cxxopts::value<double>())
	        ("offset3", "Frequency offset of the 3rd path [Hz]", cxxopts::value<double>());

        options.parse_positional({"input_file", "output_file", "positional"});

	    auto result = options.parse(argc, argv);

        if (result.count("help")) {
			std::cout << options.help({"", "Propagation condition", "Propagation"}) << std::endl;
			exit(0);
		}

        if (argc > 1 || 
            ((result.count("input_file") != 1 || result.count("output_file") != 1) &&
             result.arguments().size() != 2))  {
            std::cerr << "pathsim: input / output not specified" << std::endl;
            exit(-1);
        }

        std::string input_file  = result.count("input_file")  > 0 ? result["input_file"] .as<std::string>() : result.arguments().front().value();
        std::string output_file = result.count("output_file") > 0 ? result["output_file"].as<std::string>() : result.arguments()[1].value();

        PathSimParams params;
        // Parse the propagation condition parameter sets.
        for (int i = 0; i < default_params_cnt; ++ i)
            if (result.count(defaut_params[i].cmdline_param) > 0) {
                if (params.cmdline_param.empty())
                    params = defaut_params[i];
                else {
                    std::cerr << "pathsim: parameter " << defaut_params[i].cmdline_param << " overrides " << params.cmdline_param << std::endl;
                    std::cerr << "Use just one propagation condition parameter.";
                    return -1;
                }
            }

        if (result.count("snr")) {
            params.has_awgn = true;
            params.snr = result["snr"].as<double>();
        }

        if (result.count("spread")) {
            params.has_path0 = true;
            params.spread0   = result["spread"].as<double>();
        }
        if (result.count("offset")) {
            params.has_path0 = true;
            params.offset0   = result["offset"].as<double>();
        }

        if (result.count("delay2")) {
            params.has_path1 = true;
            params.delay2   = result["delay2"].as<double>();
        }
        if (result.count("spread2")) {
            params.has_path1 = true;
            params.spread1   = result["spread2"].as<double>();
        }
        if (result.count("offset2")) {
            params.has_path1 = true;
            params.offset1   = result["offset2"].as<double>();
        }

        if (result.count("delay2")) {
            params.has_path2 = true;
            params.delay2   = result["delay2"].as<double>();
        }
        if (result.count("spread2")) {
            params.has_path2 = true;
            params.spread2   = result["spread2"].as<double>();
        }
        if (result.count("offset2")) {
            params.has_path2 = true;
            params.offset2   = result["offset2"].as<double>();
        }

        AudioFile<double> audio_file;
        bool loaded  = audio_file.load(input_file);
        int  len     = audio_file.getNumSamplesPerChannel();
        int  nblocks = (len + PathSimProcessor::BUF_SIZE - 1) / PathSimProcessor::BUF_SIZE;

        PathSimProcessor processor;
        processor.init(params);
        audio_file.samples.front().resize(PathSimProcessor::BUF_SIZE * nblocks, 0.);
        for (int i = 0; i < nblocks; ++ i)
            processor.process_buffer(audio_file.samples.front().data() + i * PathSimProcessor::BUF_SIZE);
        audio_file.samples.front().resize(len);
        for (int k = 1; k < audio_file.getNumChannels(); ++ k)
            audio_file.samples[k] = audio_file.samples.front();
        audio_file.save(result["output_file"].as<std::string>(), AudioFileFormat::Wave);
	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}

	return 0;
}

#if 0
//////////////////////////////////////////////////////////////////////
// Called when the process has new FFT data to display every .256 seconds.
//////////////////////////////////////////////////////////////////////
afx_msg LRESULT CPathSimView::OnDataRdy(UINT x, long y)
{
	CPathSimDoc* pDoc;
	std::string str;
	double Ngain, Sgain;
	pDoc = GetDocument();
	ASSERT_VALID( pDoc );
	m_pCIOCntrl->ReadGains(Ngain, Sgain);
	str.Format("%g",Ngain );
	m_AWGNGainCtrl.SetWindowText( str);
	str.Format("%g", Sgain );
	m_SignalGainCtrl.SetWindowText( str);
	str.Format("%g", pDoc->m_SigRMS );
	m_SigRMSCtrl.SetWindowText( str);
	str.Format("%g", pDoc->m_NoiseRMS );
	m_NoiseRMSCtrl.SetWindowText( str);

#ifdef TESTMODE
	
	str.Format("%g", gDebug1);
	m_TestTextCtrl1.SetWindowText( str);
	str.Format("%g", gDebug2);
	m_TestTextCtrl2.SetWindowText( str);
#endif
	//go draw the data plots
	if(m_pCPlotData && ( pDoc->m_SimState==SIMSTATE_ON ))
		m_pCPlotData->DrawPlot();
	return 0;
}
#endif
