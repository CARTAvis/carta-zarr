/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_ERROR_H_
#define CARTA_ZARR_ERROR_H_

#include <string>

namespace carta::zarr {

enum class ErrorCode {
    not_found,
    not_zarr,
    unsupported_transport,
    unsupported_zarr_version,
    unsupported_schema,
    unsupported_schema_version,
    ambiguous_schema,
    invalid_argument,
    invalid_metadata,
    unsupported_data_type,
    unsupported_codec,
    invalid_slice,
    buffer_too_small,
    io_error,
    decode_error,
    cancelled,
    not_implemented,
};

// The code's name, as it would be spelled in a message or a log. DiagnosticCodeName says the same of
// a diagnostic's, and the codes the two share are spelled the same way.
inline const char* ErrorCodeName(ErrorCode code) noexcept {
    switch (code) {
        case ErrorCode::not_found:
            return "not_found";
        case ErrorCode::not_zarr:
            return "not_zarr";
        case ErrorCode::unsupported_transport:
            return "unsupported_transport";
        case ErrorCode::unsupported_zarr_version:
            return "unsupported_zarr_version";
        case ErrorCode::unsupported_schema:
            return "unsupported_schema";
        case ErrorCode::unsupported_schema_version:
            return "unsupported_schema_version";
        case ErrorCode::ambiguous_schema:
            return "ambiguous_schema";
        case ErrorCode::invalid_argument:
            return "invalid_argument";
        case ErrorCode::invalid_metadata:
            return "invalid_metadata";
        case ErrorCode::unsupported_data_type:
            return "unsupported_data_type";
        case ErrorCode::unsupported_codec:
            return "unsupported_codec";
        case ErrorCode::invalid_slice:
            return "invalid_slice";
        case ErrorCode::buffer_too_small:
            return "buffer_too_small";
        case ErrorCode::io_error:
            return "io_error";
        case ErrorCode::decode_error:
            return "decode_error";
        case ErrorCode::cancelled:
            return "cancelled";
        case ErrorCode::not_implemented:
            return "not_implemented";
    }
    return "unknown";
}

struct Error {
    // Always set where the library makes one. Initialized all the same, so an Error a caller declares
    // and fills in later is never read with a code nobody gave it; it is the first code, as `Error{}`
    // already made it.
    ErrorCode code = ErrorCode::not_found;
    std::string message;
    std::string node_path;
};

}  // namespace carta::zarr

#endif  // CARTA_ZARR_ERROR_H_
