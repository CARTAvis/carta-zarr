/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "value_reader.h"

#include "store_context.h"

#include <tensorstore/array.h>
#include <tensorstore/cast.h>
#include <tensorstore/context.h>
#include <tensorstore/data_type.h>
#include <tensorstore/tensorstore.h>
#include <tensorstore/util/result.h>

#include <exception>
#include <string>

namespace carta::zarr::internal::zarr {
namespace {

Error MakeError(ErrorCode code, std::string message, std::string node_path = {}) {
    return Error{code, std::move(message), std::move(node_path)};
}

}  // namespace

Result<std::vector<double>> ReadNumericValues(const std::filesystem::path& array_directory,
                                              const StoreContextPtr& context,
                                              std::string_view node) {
    try {
        auto opened = OpenZarrArray(array_directory, context, node);
        if (!opened) {
            return opened.error();
        }

        // Every real Zarr type a coordinate may be stored in is read as double. Converting rather
        // than requiring float64 is what the probe already promises: it accepts any real type, and
        // an image whose coordinates are float32 or integer has to open rather than fail here.
        // Conversion rides the read's own copy, as it does for pixels.
        auto converted = tensorstore::Cast(opened.value(), tensorstore::dtype_v<double>);
        if (!converted.ok()) {
            return MakeError(ErrorCode::unsupported_data_type,
                             "Array is not readable as double: " + converted.status().ToString(), std::string(node));
        }

        auto read_result = tensorstore::Read(converted.value()).result();
        if (!read_result.ok()) {
            return MakeError(ErrorCode::io_error, "TensorStore read failed: " + read_result.status().ToString(),
                             std::string(node));
        }

        const auto& array = read_result.value();
        std::vector<double> result;
        result.reserve(array.num_elements());
        tensorstore::IterateOverArrays([&result](const double* val) { result.push_back(*val); }, tensorstore::c_order,
                                       array);
        return result;
    } catch (const std::exception& e) {
        return MakeError(ErrorCode::io_error, e.what(), std::string(node));
    }
}

}  // namespace carta::zarr::internal::zarr
