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
		    for (const PathSimParams &params : default_params())
                group(params.cmdline_param, params.title);
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
        for (const PathSimParams &p : default_params())
            if (result.count(p.cmdline_param) > 0) {
                if (params.cmdline_param.empty())
                    params = p;
                else {
                    std::cerr << "pathsim: parameter " << p.cmdline_param << " overrides " << params.cmdline_param << std::endl;
                    std::cerr << "Use just one propagation condition parameter.";
                    return -1;
                }
            }

        if (result.count("snr"))
            params.noise = { true, result["snr"].as<double>() };

        for (int i = 0; i < 3; ++ i) {
            auto init_path = [i, &params](){
                if (int(params.paths.size()) <= i)
                    params.paths.push_back({});
            };
            std::string sidx;
            if (i > 0)
                sidx = std::to_string(i);
            std::string delay = std::string("delay") + sidx;
            if (result.count(delay)) {
                init_path();
                params.paths.back().delay = result[delay].as<double>();
            }
            std::string spread = std::string("spread") + sidx;
            if (result.count(spread)) {
                init_path();
                params.paths.back().spread = result[spread].as<double>();
            }
            std::string offset = std::string("offset") + sidx;
            if (result.count(offset)) {
                init_path();
                params.paths.back().offset = result[offset].as<double>();
            }
        }

        AudioFile<double> audio_file;
        bool loaded  = audio_file.load(input_file);
        if (! loaded) {
            std::cout << "failed loading input file: " << input_file << std::endl;
            exit(1);
        }
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
