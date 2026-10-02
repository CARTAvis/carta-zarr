/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "observation.h"

#include "attributes.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace carta::zarr::internal::xradio {

namespace {

constexpr double kSecondsPerDay = 86400.0;
// MJD of 1970-01-01, the day days_from_civil counts from.
constexpr double kUnixEpochMjd = 40587.0;

// Days from 1970-01-01 to a proleptic Gregorian date (Howard Hinnant's days_from_civil).
std::int64_t DaysFromCivil(std::int64_t year, unsigned month, unsigned day) {
    year -= month <= 2 ? 1 : 0;
    const std::int64_t era = (year >= 0 ? year : year - 399) / 400;
    const auto year_of_era = static_cast<unsigned>(year - (era * 400));
    const unsigned day_of_year = (((153 * (month > 2 ? month - 3 : month + 9)) + 2) / 5) + day - 1;
    const unsigned day_of_era = (year_of_era * 365) + (year_of_era / 4) - (year_of_era / 100) + day_of_year;
    return (era * 146097) + static_cast<std::int64_t>(day_of_era) - 719468;
}

bool IsLeapYear(std::int64_t year) {
    return (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
}

// A date read left to right, each step failing on anything that is not where it should be.
class DateReader {
public:
    explicit DateReader(std::string_view text) : _text(text) {}

    bool done() const { return _at == _text.size(); }

    // Takes `wanted` when it is next.
    bool Take(char wanted) {
        if (_at < _text.size() && _text.at(_at) == wanted) {
            ++_at;
            return true;
        }
        return false;
    }

    // Exactly `width` decimal digits.
    std::optional<unsigned> Digits(std::size_t width) {
        if (_at + width > _text.size()) {
            return std::nullopt;
        }
        unsigned value = 0;
        for (std::size_t i = 0; i < width; ++i) {
            const char digit = _text.at(_at + i);
            if (digit < '0' || digit > '9') {
                return std::nullopt;
            }
            value = (value * 10) + static_cast<unsigned>(digit - '0');
        }
        _at += width;
        return value;
    }

    // One or more decimal digits read as the fraction they spell after a point.
    std::optional<double> Fraction() {
        double value = 0.0;
        double scale = 0.1;
        const std::size_t first = _at;
        for (; _at < _text.size() && _text.at(_at) >= '0' && _text.at(_at) <= '9'; ++_at) {
            value += (_text.at(_at) - '0') * scale;
            scale /= 10.0;
        }
        return _at == first ? std::nullopt : std::optional<double>(value);
    }

private:
    std::string_view _text;
    std::size_t _at = 0;
};

// The MJD an ISO 8601 date names, read in whatever scale it was written in: YYYY-MM-DD, then
// optionally a T or a space and hh:mm, :ss and a fraction, and a trailing Z. A zone offset is not
// applied but refused, as is anything else, because an epoch guessed at is worse than none.
std::optional<double> MjdOfIsoDate(std::string_view text) {
    DateReader date(text);
    const auto year = date.Digits(4);
    const auto month = date.Take('-') ? date.Digits(2) : std::nullopt;
    const auto day = date.Take('-') ? date.Digits(2) : std::nullopt;
    if (!year || !month || !day || *month < 1 || *month > 12 || *day < 1) {
        return std::nullopt;
    }
    constexpr std::array<unsigned, 12> kDaysInMonth{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    const unsigned days_in_month = kDaysInMonth.at(*month - 1) + (*month == 2 && IsLeapYear(*year) ? 1 : 0);
    if (*day > days_in_month) {
        return std::nullopt;
    }

    double seconds = 0.0;
    if (date.Take('T') || date.Take(' ')) {
        const auto hour = date.Digits(2);
        const auto minute = date.Take(':') ? date.Digits(2) : std::nullopt;
        if (!hour || !minute || *hour > 23 || *minute > 59) {
            return std::nullopt;
        }
        seconds = (*hour * 3600.0) + (*minute * 60.0);
        if (date.Take(':')) {
            const auto second = date.Digits(2);
            // 60 is a leap second, which UTC has.
            if (!second || *second > 60) {
                return std::nullopt;
            }
            seconds += *second;
            if (date.Take('.')) {
                const auto fraction = date.Fraction();
                if (!fraction) {
                    return std::nullopt;
                }
                seconds += *fraction;
            }
        }
        date.Take('Z');
    }
    if (!date.done()) {
        return std::nullopt;
    }
    return static_cast<double>(DaysFromCivil(*year, *month, *day)) + kUnixEpochMjd + (seconds / kSecondsPerDay);
}

// The MJD a number names, by the format written beside it: XRADIO writes "mjd", and a number with no
// format has always been read as one. Seconds since 1970 are what XRADIO's time coordinate holds.
std::optional<double> MjdOfNumber(double value, const std::string& format) {
    if (format.empty() || format == "MJD") {
        return value;
    }
    if (format == "UNIX") {
        return kUnixEpochMjd + (value / kSecondsPerDay);
    }
    return std::nullopt;
}

}  // namespace

ObservationInfo DescribeObservation(const zarr::ArrayMetadata& image) {
    ObservationInfo observation;
    observation.object_name = AttributeString(image.attributes, "object_name");
    observation.observer = AttributeString(image.attributes, "observer");
    if (const auto* const telescope = MemberObject(image.attributes, "telescope")) {
        observation.telescope_name = AttributeString(*telescope, "name");

        const auto* const direction = MemberObject(*telescope, "direction");
        const auto* const distance = MemberObject(*telescope, "distance");
        const auto* const direction_data = direction == nullptr ? nullptr : MemberArray(*direction, "data");
        const auto* const distance_data = distance == nullptr ? nullptr : MemberArray(*distance, "data");
        // Every element is checked before it is converted, not just the array around it: a value of
        // the wrong type throws out of nlohmann, and this is optional metadata reached while
        // describing an image, where the caller is holding a Result and expecting a diagnostic at
        // worst rather than an exception.
        if (direction_data != nullptr && direction_data->size() >= 2 && direction_data->at(0).is_number() &&
            direction_data->at(1).is_number() && distance_data != nullptr && !distance_data->empty() &&
            distance_data->at(0).is_number()) {
            const double lon = direction_data->at(0).get<double>();
            const double lat = direction_data->at(1).get<double>();
            const double radius = distance_data->at(0).get<double>();
            observation.observatory_position = std::array<double, 3>{
                radius * std::cos(lat) * std::cos(lon), radius * std::cos(lat) * std::sin(lon), radius * std::sin(lat)};
        }
    }
    if (const auto* const obsdate = MemberObject(image.attributes, "obsdate")) {
        std::string format;
        if (const auto* const attributes = MemberObject(*obsdate, "attrs")) {
            observation.timesys = Upper(AttributeString(*attributes, "scale"));
            format = Upper(AttributeString(*attributes, "format"));
        }
        // A date is kept as written and also reported as the MJD it names, in the same scale, so
        // that a consumer building an epoch has one field to read however the writer spelled it.
        if (const auto* const data = Member(*obsdate, "data"); data != nullptr) {
            if (data->is_string()) {
                observation.date_obs = data->get<std::string>();
                observation.mjd_obs = MjdOfIsoDate(observation.date_obs);
            } else if (data->is_number()) {
                observation.mjd_obs = MjdOfNumber(data->get<double>(), format);
            }
        }
    }
    return observation;
}

}  // namespace carta::zarr::internal::xradio
