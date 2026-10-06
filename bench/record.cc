/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "record.h"

#include "mode.h"

#include <nlohmann/json.hpp>

#include <array>
#include <cinttypes>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <set>
#include <stdexcept>
#include <system_error>
#include <utility>

#include <unistd.h>

#ifndef CARTA_ZARR_BENCH_COMMIT
#define CARTA_ZARR_BENCH_COMMIT "unknown"
#endif
#ifndef CARTA_ZARR_BENCH_BUILD_TYPE
#define CARTA_ZARR_BENCH_BUILD_TYPE "unknown"
#endif
#ifndef CARTA_ZARR_BENCH_TUNING
#define CARTA_ZARR_BENCH_TUNING "unknown"
#endif

namespace carta::zarr::bench {

namespace {

std::string Text(double value) {
    std::array<char, 32> buffer{};
    std::snprintf(buffer.data(), buffer.size(), "%.6f", value);
    return buffer.data();
}

template <typename T>
std::string Text(const std::optional<T>& value) {
    if (!value) {
        return "";
    }
    if constexpr (std::is_floating_point_v<T>) {
        return Text(*value);
    } else {
        return std::to_string(*value);
    }
}

std::string Hex(std::uint64_t value) {
    std::array<char, 17> buffer{};
    std::snprintf(buffer.data(), buffer.size(), "%016" PRIx64, value);
    return buffer.data();
}

std::string Bool(bool value) {
    return value ? "true" : "false";
}

using Column = std::pair<const char*, std::function<std::string(const Row&)>>;

// Every column, in order. The header is the first of each pair and a row the second.
const std::vector<Column>& Columns() {
    static const std::vector<Column> columns{
        {"csv_version", [](const Row&) { return std::to_string(kCsvVersion); }},
        {"run_id", [](const Row& r) { return r.run_id; }},
        {"run_key", [](const Row& r) { return r.run_key; }},
        {"label", [](const Row& r) { return r.label; }},
        {"timestamp_utc", [](const Row& r) { return r.timestamp_utc; }},
        {"host", [](const Row& r) { return r.host; }},
        {"bench_commit", [](const Row& r) { return r.bench_commit; }},
        {"build_type", [](const Row& r) { return r.build_type; }},
        {"tuning", [](const Row& r) { return r.tuning; }},
        {"dataset", [](const Row& r) { return r.dataset; }},
        {"dataset_identity_hash", [](const Row& r) { return r.dataset_identity_hash; }},
        {"image_id", [](const Row& r) { return r.image_id; }},
        {"shape", [](const Row& r) { return r.shape; }},
        {"chunk_shape", [](const Row& r) { return r.chunk_shape; }},
        {"shard_shape", [](const Row& r) { return r.shard_shape; }},
        {"chunk_decoded_bytes", [](const Row& r) { return std::to_string(r.chunk_decoded_bytes); }},
        {"codec", [](const Row& r) { return r.codec; }},
        {"stripe_effective", [](const Row& r) { return r.stripe_effective; }},
        {"compression_ratio", [](const Row& r) { return r.compression_ratio; }},
        {"mode", [](const Row& r) { return r.mode; }},
        {"processes", [](const Row& r) { return std::to_string(r.processes); }},
        {"io_threads", [](const Row& r) { return std::to_string(r.io_threads); }},
        {"decode_threads", [](const Row& r) { return std::to_string(r.decode_threads); }},
        {"cache_bytes", [](const Row& r) { return r.cache_bytes; }},
        {"read_budget_bytes", [](const Row& r) { return std::to_string(r.read_budget_bytes); }},
        {"seed", [](const Row& r) { return std::to_string(r.seed); }},
        {"ops", [](const Row& r) { return std::to_string(r.ops); }},
        {"region_fraction", [](const Row& r) { return r.region_fraction; }},
        {"histogram_method", [](const Row& r) { return r.histogram_method; }},
        {"animation_frames", [](const Row& r) { return r.animation_frames; }},
        {"animation_fps", [](const Row& r) { return r.animation_fps; }},
        {"animation_prefetch", [](const Row& r) { return r.animation_prefetch; }},
        {"trial_timeout_s", [](const Row& r) { return std::to_string(r.trial_timeout_s); }},
        {"cold_method", [](const Row& r) { return r.cold_method; }},
        {"cold_ok", [](const Row& r) { return Bool(r.cold_ok); }},
        {"trial", [](const Row& r) { return std::to_string(r.trial); }},
        {"process_index", [](const Row& r) { return std::to_string(r.process_index); }},
        {"op_index", [](const Row& r) { return Text(r.op_index); }},
        {"position", [](const Row& r) { return r.position; }},
        {"overlap", [](const Row& r) { return Bool(r.overlap); }},
        {"shares_chunks", [](const Row& r) { return Bool(r.shares_chunks); }},
        {"status", [](const Row& r) { return r.status; }},
        {"error", [](const Row& r) { return r.error; }},
        {"seconds", [](const Row& r) { return Text(r.seconds); }},
        {"t_start_s", [](const Row& r) { return Text(r.t_start_s); }},
        {"t_end_s", [](const Row& r) { return Text(r.t_end_s); }},
        {"elements", [](const Row& r) { return std::to_string(r.elements); }},
        {"logical_bytes", [](const Row& r) { return std::to_string(r.logical_bytes); }},
        {"storage_read_bytes", [](const Row& r) { return Text(r.storage_read_bytes); }},
        {"user_cpu_s", [](const Row& r) { return Text(r.user_cpu_s); }},
        {"sys_cpu_s", [](const Row& r) { return Text(r.sys_cpu_s); }},
        {"peak_rss_bytes", [](const Row& r) { return Text(r.peak_rss_bytes); }},
        {"checksum", [](const Row& r) { return r.checksum ? Hex(*r.checksum) : std::string(); }},
        {"frames_played", [](const Row& r) { return Text(r.frames_played); }},
        {"frame_first_s", [](const Row& r) { return Text(r.frame_first_s); }},
        {"frame_median_s", [](const Row& r) { return Text(r.frame_median_s); }},
        {"frame_max_s", [](const Row& r) { return Text(r.frame_max_s); }},
        {"late_frames", [](const Row& r) { return Text(r.late_frames); }},
        {"late_max_s", [](const Row& r) { return Text(r.late_max_s); }},
        {"prefetches", [](const Row& r) { return Text(r.prefetches); }},
        {"late_prefetches", [](const Row& r) { return Text(r.late_prefetches); }},
    };
    return columns;
}

std::size_t ColumnIndex(const char* name) {
    const auto& columns = Columns();
    for (std::size_t index = 0; index < columns.size(); ++index) {
        if (std::string_view(columns[index].first) == name) {
            return index;
        }
    }
    return columns.size();
}

// Quoted when it must be, as RFC 4180 has it. A line break becomes a space instead, so that a row
// stays a line: an error message is the only field that might hold one, and reads as well without.
std::string Field(std::string value) {
    for (auto& character : value) {
        if (character == '\n' || character == '\r') {
            character = ' ';
        }
    }
    if (value.find_first_of(",\"") == std::string::npos) {
        return value;
    }
    std::string quoted = "\"";
    for (const char character : value) {
        quoted += character;
        if (character == '"') {
            quoted += '"';
        }
    }
    return quoted + "\"";
}

std::uint64_t Fnv(const std::string& text) {
    std::uint64_t hash = 0xCBF29CE484222325ull;
    for (const unsigned char character : text) {
        hash ^= character;
        hash *= 0x100000001B3ull;
    }
    return hash;
}

std::string Host() {
    std::array<char, 256> name{};
    if (gethostname(name.data(), name.size() - 1) != 0) {
        return "";
    }
    return name.data();
}

std::optional<double> Number(const std::string& text) {
    if (text.empty()) {
        return std::nullopt;
    }
    try {
        return std::stod(text);
    } catch (...) {
        return std::nullopt;
    }
}

}  // namespace

DatasetFacts ReadManifest(const std::string& dataset) {
    DatasetFacts facts;
    std::ifstream file(std::filesystem::path(dataset) / "bench-manifest.json");
    if (!file) {
        return facts;
    }
    const auto manifest = nlohmann::json::parse(file, nullptr, false);
    if (manifest.is_discarded() || !manifest.is_object()) {
        return facts;
    }
    const auto text = [](const nlohmann::json& value) -> std::string {
        if (value.is_null()) {
            return "";
        }
        return value.is_string() ? value.get<std::string>() : value.dump();
    };
    facts.identity_hash = text(manifest.value("identity_hash", nlohmann::json()));
    if (const auto layout = manifest.find("layout"); layout != manifest.end() && layout->is_object()) {
        facts.codec = text(layout->value("codec", nlohmann::json()));
    }
    if (const auto storage = manifest.find("storage"); storage != manifest.end() && storage->is_object()) {
        facts.stripe_effective = text(storage->value("stripe_effective", nlohmann::json()));
    }
    if (const auto result = manifest.find("result"); result != manifest.end() && result->is_object()) {
        const auto ratio = result->value("compression_ratio", nlohmann::json());
        if (ratio.is_number()) {
            std::array<char, 32> buffer{};
            std::snprintf(buffer.data(), buffer.size(), "%.3f", ratio.get<double>());
            facts.compression_ratio = buffer.data();
        }
    }
    return facts;
}

Row RowTemplate(const RunOptions& options, Mode mode, ColdMethod cold, const DatasetFacts& facts,
                const std::string& run_id, const std::string& build) {
    Row row;
    row.run_id = run_id;
    row.label = options.label;
    row.host = Host();
    row.bench_commit = CARTA_ZARR_BENCH_COMMIT;
    row.build_type = CARTA_ZARR_BENCH_BUILD_TYPE;
    row.tuning = CARTA_ZARR_BENCH_TUNING;
    row.dataset = options.dataset;
    row.dataset_identity_hash = facts.identity_hash;
    row.image_id = options.image_id;
    row.codec = facts.codec;
    row.stripe_effective = facts.stripe_effective;
    row.compression_ratio = facts.compression_ratio;
    row.mode = ModeName(mode);
    row.processes = options.processes;
    row.io_threads = options.context.io_threads;
    row.decode_threads = options.context.decode_threads;
    row.cache_bytes = options.context.cache_bytes ? std::to_string(*options.context.cache_bytes) : "default";
    row.read_budget_bytes = options.read_budget_bytes;
    row.seed = options.seed;
    row.ops = options.OpsFor(mode);
    Workload::For(mode, options)->WriteSettings(row);
    row.trial_timeout_s = options.trial_timeout.count();
    row.cold_method = ColdMethodName(cold);

    // Everything that decides what a trial measures, and nothing that does not: a dataset is known by
    // what the generator says decides its bytes when there is a manifest, since the same layout
    // written again elsewhere is the same measurement, and by its path when there is not; and the
    // build measuring it by the bytes it runs, since a reader changed since is another measurement.
    std::error_code ignored;
    const auto identity = facts.identity_hash.empty()
                            ? std::filesystem::absolute(options.dataset, ignored).lexically_normal().string()
                            : facts.identity_hash;
    std::string key = std::to_string(kCsvVersion);
    for (const auto& part : {identity,
                             build,
                             row.tuning,
                             options.image_id,
                             row.mode,
                             std::to_string(row.processes),
                             std::to_string(row.io_threads),
                             std::to_string(row.decode_threads),
                             row.cache_bytes,
                             std::to_string(row.read_budget_bytes),
                             std::to_string(row.seed),
                             std::to_string(row.ops),
                             row.region_fraction,
                             row.histogram_method,
                             row.animation_frames,
                             row.animation_fps,
                             row.animation_prefetch,
                             std::to_string(row.trial_timeout_s),
                             row.cold_method,
                             row.label}) {
        key += '|';
        key += part;
    }
    row.run_key = Hex(Fnv(key));
    return row;
}

std::string CsvHeader() {
    std::string header;
    for (const auto& [name, field] : Columns()) {
        if (!header.empty()) {
            header += ',';
        }
        header += name;
    }
    return header;
}

std::string FormatRow(const Row& row) {
    std::string line;
    bool first = true;
    for (const auto& [name, field] : Columns()) {
        if (!first) {
            line += ',';
        }
        first = false;
        line += Field(field(row));
    }
    return line;
}

std::vector<std::string> ParseCsvLine(const std::string& line) {
    std::vector<std::string> fields(1);
    bool quoted = false;
    for (std::size_t index = 0; index < line.size(); ++index) {
        const char character = line[index];
        if (quoted) {
            if (character == '"' && index + 1 < line.size() && line[index + 1] == '"') {
                fields.back() += '"';
                ++index;
            } else if (character == '"') {
                quoted = false;
            } else {
                fields.back() += character;
            }
        } else if (character == '"') {
            quoted = true;
        } else if (character == ',') {
            fields.emplace_back();
        } else if (character != '\r') {
            fields.back() += character;
        }
    }
    return fields;
}

std::optional<RowSummary> SummariseRow(const std::string& line) {
    static const auto run_key = ColumnIndex("run_key");
    static const auto trial = ColumnIndex("trial");
    static const auto status = ColumnIndex("status");
    static const auto seconds = ColumnIndex("seconds");
    static const auto t_end_s = ColumnIndex("t_end_s");
    static const auto processes = ColumnIndex("processes");
    static const auto ops = ColumnIndex("ops");
    static const auto process_index = ColumnIndex("process_index");
    static const auto op_index = ColumnIndex("op_index");
    const auto fields = ParseCsvLine(line);
    if (fields.size() != Columns().size()) {
        return std::nullopt;
    }
    RowSummary summary;
    summary.run_key = fields[run_key];
    try {
        summary.trial = static_cast<unsigned>(std::stoul(fields[trial]));
        summary.processes = static_cast<unsigned>(std::stoul(fields[processes]));
        summary.ops = static_cast<unsigned>(std::stoul(fields[ops]));
        summary.process_index = static_cast<unsigned>(std::stoul(fields[process_index]));
        if (!fields[op_index].empty()) {
            summary.op_index = static_cast<unsigned>(std::stoul(fields[op_index]));
        }
    } catch (...) {
        return std::nullopt;
    }
    summary.status = fields[status];
    summary.seconds = Number(fields[seconds]);
    summary.t_end_s = Number(fields[t_end_s]);
    return summary;
}

namespace {

// Cut a file back to its last newline. The bench ends every line it writes with one, so a last line
// without it is a write that stopped part-way -- a row, or the header, of which only some bytes
// reached the disk. It is no row at all, and a row appended after it would join it into one line
// that is neither. False, with `error` set, only when the file cannot be read or cut.
bool DropLineCutShort(const std::string& path, std::string& error) {
    std::error_code failed;
    if (!std::filesystem::exists(path, failed)) {
        return true;
    }
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        error = "cannot read " + path;
        return false;
    }
    const std::string text{std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
    if (text.empty() || text.back() == '\n') {
        return true;
    }
    const auto last_newline = text.find_last_of('\n');
    const auto keep = last_newline == std::string::npos ? std::uintmax_t{0} : std::uintmax_t{last_newline + 1};
    std::filesystem::resize_file(path, keep, failed);
    if (failed) {
        error = "cannot drop the line " + path + " was cut short in: " + failed.message();
        return false;
    }
    std::fprintf(stderr, "note: dropped the last line of %s, which a write stopped part-way through\n", path.c_str());
    return true;
}

}  // namespace

std::optional<CsvOutput> CsvOutput::Open(const std::string& path, std::string& error) {
    if (path.empty()) {
        CsvOutput output(stdout, {});
        std::fprintf(stdout, "%s\n", CsvHeader().c_str());
        std::fflush(stdout);
        return output;
    }

    // What each trial's rows say about it: whether any was an error, how many operations it was
    // meant to have, and which of them have a row.
    struct Coverage {
        bool error = false;
        std::optional<std::pair<unsigned, unsigned>> grid;
        std::set<std::pair<unsigned, unsigned>> reported;
    };
    std::map<std::string, std::map<unsigned, Coverage>> trials;
    bool has_header = false;
    if (!DropLineCutShort(path, error)) {
        return std::nullopt;
    }
    if (std::ifstream existing(path); existing) {
        std::string line;
        if (std::getline(existing, line) && !line.empty()) {
            if (line != CsvHeader()) {
                error = path + " was written by another version of carta-zarr-bench: its header differs";
                return std::nullopt;
            }
            has_header = true;
        }
        for (std::size_t number = 2; std::getline(existing, line); ++number) {
            const auto row = SummariseRow(line);
            // Every line left ends in a newline, so none was cut short by a write stopping; one that
            // is still not a row is something else, and appending to the file would bury it.
            if (!row) {
                error = "line " + std::to_string(number) + " of " + path +
                        " is not a row carta-zarr-bench writes; repair the file or move it aside";
                return std::nullopt;
            }
            auto& coverage = trials[row->run_key][row->trial];
            const std::pair<unsigned, unsigned> grid{row->processes, row->ops};
            // Rows of one trial disagreeing about its size cannot all be its rows.
            coverage.error =
                coverage.error || row->status == "error" || !row->op_index || (coverage.grid && *coverage.grid != grid);
            coverage.grid = grid;
            if (row->op_index && row->process_index < row->processes && *row->op_index < row->ops) {
                coverage.reported.emplace(row->process_index, *row->op_index);
            }
        }
    }

    std::FILE* file = std::fopen(path.c_str(), "a");
    if (file == nullptr) {
        error = "cannot open " + path + " for appending";
        return std::nullopt;
    }
    CsvOutput output(file, path);
    if (!has_header) {
        std::fprintf(file, "%s\n", CsvHeader().c_str());
        std::fflush(file);
    }
    for (const auto& [key, outcomes] : trials) {
        for (const auto& [trial, coverage] : outcomes) {
            const auto expected =
                coverage.grid ? std::uint64_t{coverage.grid->first} * coverage.grid->second : std::uint64_t{0};
            if (!coverage.error && expected > 0 && coverage.reported.size() == expected) {
                output._completed[key].insert(trial);
            } else {
                output._partial.emplace(key, trial);
            }
        }
    }
    return output;
}

CsvOutput::CsvOutput(CsvOutput&& other) noexcept
    : _file(std::exchange(other._file, nullptr)),
      _path(std::move(other._path)),
      _completed(std::move(other._completed)),
      _partial(std::move(other._partial)) {}

CsvOutput::~CsvOutput() {
    if (_file != nullptr && _file != stdout) {
        std::fclose(_file);
    }
}

void CsvOutput::Remove(const std::string& run_key, unsigned trial) {
    std::fclose(_file);
    _file = nullptr;
    // Written beside the file and moved over it, so that a run stopped part-way leaves the old file
    // whole rather than half of a new one.
    const auto replacement = _path + ".partial";
    {
        std::ifstream existing(_path);
        std::ofstream kept(replacement, std::ios::trunc);
        for (std::string line; std::getline(existing, line);) {
            const auto row = SummariseRow(line);
            if (!row || row->run_key != run_key || row->trial != trial) {
                kept << line << '\n';
            }
        }
        kept.flush();
        if (!kept) {
            throw std::runtime_error("cannot write " + replacement);
        }
    }
    std::filesystem::rename(replacement, _path);
    _file = std::fopen(_path.c_str(), "a");
    if (_file == nullptr) {
        throw std::runtime_error("cannot reopen " + _path + " for appending");
    }
}

void CsvOutput::Write(const std::string& run_key, unsigned trial, const std::vector<std::string>& rows) {
    if (_partial.erase({run_key, trial}) > 0) {
        Remove(run_key, trial);
    }
    for (const auto& row : rows) {
        std::fprintf(_file, "%s\n", row.c_str());
    }
    std::fflush(_file);
}

std::size_t ItemSize(DataType type) {
    switch (type) {
        case DataType::boolean:
        case DataType::int8:
        case DataType::uint8:
            return 1;
        case DataType::int16:
        case DataType::uint16:
        case DataType::float16:
            return 2;
        case DataType::int32:
        case DataType::uint32:
        case DataType::float32:
            return 4;
        case DataType::int64:
        case DataType::uint64:
        case DataType::float64:
        case DataType::complex64:
            return 8;
        case DataType::complex128:
            return 16;
        case DataType::unknown:
            break;
    }
    return 0;
}

namespace {

std::string Shape(const std::vector<AxisDescriptor>& axes, const std::vector<std::uint64_t>& lengths) {
    std::string text;
    for (std::size_t index = 0; index < axes.size() && index < lengths.size(); ++index) {
        if (!text.empty()) {
            text += ';';
        }
        text += axes[index].name + "=" + std::to_string(lengths[index]);
    }
    return text;
}

}  // namespace

void DescribeImage(Row& row, const Image& image) {
    const auto& descriptor = image.descriptor();
    const auto& geometry = image.chunk_geometry();
    std::vector<std::uint64_t> lengths;
    for (const auto& axis : descriptor.axes) {
        lengths.push_back(axis.length);
    }
    row.image_id = descriptor.id;
    row.shape = Shape(descriptor.axes, lengths);
    row.chunk_shape = Shape(descriptor.axes, geometry.chunk_shape);
    row.shard_shape = geometry.sharded ? Shape(descriptor.axes, geometry.shard_shape) : "";
    // The bench reads with ReadOptions' defaults but for its deadline and budget, neither of which
    // changes what a chunk decodes to.
    row.chunk_decoded_bytes = image.DecodedChunkBytes();
    if (row.codec.empty()) {
        row.codec = geometry.compressor;
    }
}

std::string Describe(const Error& error) {
    std::string text = std::string(ErrorCodeName(error.code)) + ": " + error.message;
    if (!error.node_path.empty()) {
        text += " (" + error.node_path + ")";
    }
    return text;
}

}  // namespace carta::zarr::bench
