/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// options
#include "traccc/options/common_options.hpp"

traccc::common_options::common_options(
    po::options_description& desc){

    desc.add_options()("input-csv", po::value<bool>()->default_value(false),
                       "Use csv input file")(
        "input-binary", po::value<bool>()->default_value(false),
        "Use binary input file");
    desc.add_options()("events", po::value<unsigned int>()->required(),
                       "number of events");
    desc.add_options()("skip", po::value<int>()->default_value(0),
                       "number of events to skip");
}

void traccc::common_options::read(const po::variables_map& vm) {    

    if (vm["input-csv"].as<bool>() == true) {
        data_format == traccc::data_format::csv;
    } else if (vm["input-binary"].as<bool>() == true) {
        data_format == traccc::data_format::binary;
    }
    events = vm["events"].as<unsigned int>();
    skip = vm["skip"].as<int>();
}