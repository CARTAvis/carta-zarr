/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

// Result's accessors, as a consumer reaches them. They are named after std::expected's so that a
// reader already knows what each does; this pins that each does it.

#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <variant>

#include <carta-zarr/descriptor.h>
#include <carta-zarr/result.h>

#include "support/check.h"

namespace {

using carta::zarr::Error;
using carta::zarr::ErrorCode;
using carta::zarr::Result;
using carta::zarr::testing::Require;

Error AnError() {
    return Error{ErrorCode::io_error, "the read failed", "SKY"};
}

void TestAValue() {
    Result<std::string> result{std::string("pixels")};
    Require(result.has_value() && static_cast<bool>(result), "a value did not report itself");
    Require(*result == "pixels" && result->size() == 6 && result.value() == "pixels",
            "*, -> and value() did not reach the value");
    Require(result.value_or("fallback") == "pixels", "value_or replaced a value that was there");

    *result = "changed";
    Require(result.value() == "changed", "* did not reach the value itself");
}

void TestAnError() {
    Result<std::string> result{AnError()};
    Require(!result.has_value() && !static_cast<bool>(result), "an error reported a value");
    Require(result.error().code == ErrorCode::io_error && result.error().node_path == "SKY",
            "error() did not reach the error");
    Require(result.value_or("fallback") == "fallback", "value_or did not fall back on an error");

    bool threw = false;
    try {
        static_cast<void>(result.value());
    } catch (const std::bad_variant_access&) {
        threw = true;
    }
    Require(threw, "value() on an error did not throw, as std::expected's does");
}

// The rvalue overloads move rather than copy, which a move-only type shows: copying one would not
// compile, and a copy that happened to compile would leave the source's pointer where it was.
void TestTheRvalueOverloadsMove() {
    Result<std::unique_ptr<int>> held{std::make_unique<int>(7)};
    const auto taken = std::move(held).value_or(nullptr);
    Require(taken && *taken == 7, "value_or on an rvalue did not move the value out");

    Result<std::unique_ptr<int>> absent{AnError()};
    Require(std::move(absent).value_or(nullptr) == nullptr, "value_or on an rvalue did not fall back");

    Result<int> failed{AnError()};
    const Error moved = std::move(failed).error();
    Require(moved.message == "the read failed", "error() on an rvalue did not hand the error over");
}

void TestVoid() {
    const Result<void> done;
    Require(done.has_value() && static_cast<bool>(done), "a default Result<void> did not report success");

    Result<void> failed{AnError()};
    Require(!failed.has_value() && failed.error().code == ErrorCode::io_error,
            "a Result<void> did not carry its error");
    const Error moved = std::move(failed).error();
    Require(moved.node_path == "SKY", "error() on an rvalue Result<void> did not hand the error over");
}

// Every code has a name of its own, and a code a diagnostic shares is spelled as the diagnostic's is:
// the two are what a consumer writes into one log.
void TestErrorCodeNames() {
    std::set<std::string> names;
    for (int code = 0; code <= static_cast<int>(ErrorCode::not_implemented); ++code) {
        const std::string name = carta::zarr::ErrorCodeName(static_cast<ErrorCode>(code));
        Require(name != "unknown", "error code " + std::to_string(code) + " has no name");
        Require(names.insert(name).second, "two error codes are both named " + name);
    }
    using carta::zarr::DiagnosticCode;
    for (const auto& [error, diagnostic] :
         {std::pair{ErrorCode::invalid_metadata, DiagnosticCode::invalid_metadata},
          std::pair{ErrorCode::unsupported_data_type, DiagnosticCode::unsupported_data_type},
          std::pair{ErrorCode::ambiguous_schema, DiagnosticCode::ambiguous_schema}}) {
        Require(std::string(carta::zarr::ErrorCodeName(error)) == carta::zarr::DiagnosticCodeName(diagnostic),
                "an error and a diagnostic that share a code spell it differently");
    }
}

}  // namespace

int main() {
    try {
        TestAValue();
        TestAnError();
        TestTheRvalueOverloadsMove();
        TestVoid();
        TestErrorCodeNames();
        std::cout << "carta-zarr result tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "carta-zarr result tests failed: " << error.what() << '\n';
        return 1;
    }
}
