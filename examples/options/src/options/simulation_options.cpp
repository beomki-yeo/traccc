/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// options
#include "traccc/options/simulation_options.hpp"

traccc::simulation_options::simulation_options(po::options_description& desc) {

    desc.add_options()("run_digitization",
                       po::value<bool>()->default_value(false),
                       "run digitization to generate cell file");
    desc.add_options()("digitization_config_file", po::value<std::string>(),
                       "specify digitization config file");
}

void traccc::simulation_options::read(const po::variables_map& vm) {

    run_digitization = vm["run_digitization"].as<bool>();
    if (vm.count("digitization_config_file")) {
        digitization_config_file =
            vm["digitization_config_file"].as<std::string>();
    }
}
