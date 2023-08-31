/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Boost
#include <boost/program_options.hpp>

namespace traccc {

namespace po = boost::program_options;

struct simulation_options {
    bool run_digitization;
    std::string digitization_config_file;

    simulation_options(po::options_description& desc);
    void read(const po::variables_map& vm);
};

}  // namespace traccc