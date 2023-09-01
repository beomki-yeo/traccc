/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/experimental/read_cells.hpp"

#include "../csv/experimental/read_cells.hpp"
#include "traccc/io/utils.hpp"

namespace traccc::io::experimental {

void read_cells(cell_reader_output &out, std::size_t event,
                std::string_view directory, data_format format,
                const digitization_map *dconfig) {

    switch (format) {
        case data_format::csv: {
            read_cells(out,
                       data_directory() + directory.data() +
                           get_event_filename(event, "-cells.csv"),
                       format, dconfig);
            break;
        }

        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

void read_cells(cell_reader_output &out, std::string_view filename,
                data_format format, const digitization_map *dconfig) {

    switch (format) {
        case data_format::csv:
            return csv::experimental::read_cells(out, filename, dconfig);

        default:
            throw std::invalid_argument("Unsupported data format");
    }
}

}  // namespace traccc::io::experimental
