/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef CARTA_ZARR_RESULT_H_
#define CARTA_ZARR_RESULT_H_

#include "carta-zarr/error.h"

#include <cassert>
#include <optional>
#include <type_traits>
#include <utility>
#include <variant>

namespace carta::zarr {

// A value or the Error that stands in for it.
//
// Named and shaped after C++23's std::expected so that a consumer already knows how to read one, and
// so that replacing it with that type later is a rename: has_value, value, error, value_or, operator*
// and operator-> mean what they mean there. value() on an error throws std::bad_variant_access, as
// std::expected's throws bad_expected_access; operator*, operator-> and error() on the wrong
// alternative are preconditions, checked by assert.
//
// [[nodiscard]] because every entry point of this library reports failure here and nowhere else. A
// reduction returns Result<void>, and a call whose result is dropped is a failure nobody hears
// about: the sink simply stops being called, which reads as a reduction that had nothing more to say.
template <typename T>
class [[nodiscard]] Result {
public:
    Result(T value) : _value(std::move(value)) {}
    Result(Error error) : _value(std::move(error)) {}

    bool has_value() const noexcept { return std::holds_alternative<T>(_value); }
    explicit operator bool() const noexcept { return has_value(); }

    T& value() & { return std::get<T>(_value); }
    const T& value() const& { return std::get<T>(_value); }
    T&& value() && { return std::get<T>(std::move(_value)); }

    T& operator*() & noexcept {
        assert(has_value());
        return *std::get_if<T>(&_value);
    }
    const T& operator*() const& noexcept {
        assert(has_value());
        return *std::get_if<T>(&_value);
    }
    T&& operator*() && noexcept {
        assert(has_value());
        return std::move(*std::get_if<T>(&_value));
    }
    T* operator->() noexcept {
        assert(has_value());
        return std::get_if<T>(&_value);
    }
    const T* operator->() const noexcept {
        assert(has_value());
        return std::get_if<T>(&_value);
    }

    template <typename U>
    T value_or(U&& fallback) const& {
        return has_value() ? **this : static_cast<T>(std::forward<U>(fallback));
    }
    template <typename U>
    T value_or(U&& fallback) && {
        return has_value() ? std::move(**this) : static_cast<T>(std::forward<U>(fallback));
    }

    Error& error() & noexcept {
        assert(!has_value());
        return *std::get_if<Error>(&_value);
    }
    const Error& error() const& noexcept {
        assert(!has_value());
        return *std::get_if<Error>(&_value);
    }
    Error&& error() && noexcept {
        assert(!has_value());
        return std::move(*std::get_if<Error>(&_value));
    }

private:
    static_assert(!std::is_same_v<std::decay_t<T>, Error>, "a Result of an Error cannot tell its value from its error");

    std::variant<T, Error> _value;
};

template <>
class [[nodiscard]] Result<void> {
public:
    Result() noexcept = default;
    Result(Error error) : _error(std::move(error)) {}

    bool has_value() const noexcept { return !_error.has_value(); }
    explicit operator bool() const noexcept { return has_value(); }

    Error& error() & noexcept {
        assert(!has_value());
        return *_error;
    }
    const Error& error() const& noexcept {
        assert(!has_value());
        return *_error;
    }
    Error&& error() && noexcept {
        assert(!has_value());
        return std::move(*_error);
    }

private:
    std::optional<Error> _error;
};

}  // namespace carta::zarr

#endif  // CARTA_ZARR_RESULT_H_
