/*
 * This file is part of the CARTA Image Viewer: https://github.com/CARTAvis
 * Copyright 2026 Academia Sinica Institute of Astronomy and Astrophysics (ASIAA)
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "image.h"

#include "../../zarr/array_metadata.h"
#include "attributes.h"
#include "beam_table.h"
#include "coordinates.h"
#include "direction.h"
#include "flag.h"
#include "linear_axis.h"
#include "observation.h"
#include "probe_report.h"
#include "qualification.h"

#include <carta-zarr/descriptor.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <optional>
#include <utility>

namespace carta::zarr::internal::xradio {
namespace {

namespace zarr_metadata = ::carta::zarr::internal::zarr;

constexpr std::string_view kVersion = "1.2";

constexpr int kKnownImages = 6;

int KnownImageRank(std::string_view image_id) {
    static constexpr std::array<std::string_view, kKnownImages> known{
        "SKY", "MODEL", "RESIDUAL", "POINT_SPREAD_FUNCTION", "PRIMARY_BEAM", "MASK_DECONVOLVE"};
    const auto* const found = std::find(known.begin(), known.end(), image_id);
    return found == known.end() ? static_cast<int>(known.size()) : static_cast<int>(found - known.begin());
}

// The dataset carries a well-formed array for every sky-plane coordinate. Both XRADIO readers write
// all five unconditionally (xds_from_casacore.py and xds_from_fits.py both assign coords["time"]), so
// a missing coordinate means the store is malformed.
void RequirePresentCoordinates(ProbeReport& report) {
    for (const auto& coordinate : kCoordinates) {
        if (OnPlane(coordinate, Plane::sky)) {
            report.RequireCoordinate(coordinate);
        }
    }
}

void AppendDiagnostics(ImageDescriptor& descriptor, std::vector<Diagnostic> diagnostics) {
    descriptor.diagnostics.insert(descriptor.diagnostics.end(), std::make_move_iterator(diagnostics.begin()),
                                  std::make_move_iterator(diagnostics.end()));
}

std::vector<std::string> FindDataGroups(const nlohmann::json& root_attributes, std::string_view image_id) {
    std::vector<std::string> data_groups;
    const auto* const groups = MemberObject(root_attributes, "data_groups");
    if (groups == nullptr) {
        return data_groups;
    }

    for (const auto& [group_name, group] : groups->items()) {
        if (!group.is_object()) {
            continue;
        }
        const auto references_image = std::any_of(group.begin(), group.end(), [&](const nlohmann::json& value) {
            return value.is_string() && value.get<std::string>() == image_id;
        });
        if (references_image) {
            data_groups.push_back(group_name);
        }
    }
    return data_groups;
}

// XRADIO writes the role on the variable's own "type" attribute, lowercased ("sky", "model",
// "residual"); the same attribute spells "flag" for pixel masks, which are never images. The v2
// schema proposes a separate "image_type" attribute that no released XRADIO writes yet, so it is
// only consulted as a forward-compatible fallback.
//
// One rule for the listing and for describing, so that an entry and the image it opens as cannot
// disagree about which part the image plays.
std::string ImageRoleOf(const zarr_metadata::ArrayMetadata& image) {
    const std::string declared_role = AttributeString(image.attributes, "type");
    std::string role = declared_role == "flag" ? std::string{} : declared_role;
    if (role.empty()) {
        role = AttributeString(image.attributes, "image_type");
    }
    return role;
}

std::optional<SpectralCoordinate> DescribeSpectralCoordinate(const Store& store,
                                                             const std::vector<double>& frequency_values,
                                                             ImageDescriptor& descriptor) {
    if (frequency_values.empty()) {
        return std::nullopt;
    }

    SpectralCoordinate spectral;
    spectral.channel_frequencies = frequency_values;
    nlohmann::json frequency_attributes = nlohmann::json::object();
    if (const auto& frequency_metadata = store.ReadNodeMetadata("frequency"); frequency_metadata) {
        if (const auto* const attributes = MemberObject(frequency_metadata.value(), "attributes")) {
            frequency_attributes = *attributes;
        }
    }

    // The axis reports the same unit, from the same rule.
    spectral.unit = CoordinateUnit(store, SkyCoordinate(AxisRole::spectral));
    // The channel a linear description is measured from, which is the first one unless the
    // reference frequency names another. What that measure carries -- its frame under `attrs`, its
    // value under `data` -- is read from the one lookup.
    double reference_value = spectral.channel_frequencies.front();
    if (const auto* const reference_frequency = MemberObject(frequency_attributes, "reference_frequency")) {
        if (const auto* const attributes = MemberObject(*reference_frequency, "attrs")) {
            spectral.system = Upper(AttributeString(*attributes, "observer"));
        }
        if (const auto value = MemberNumber(*reference_frequency, "data")) {
            reference_value = *value;
        }
    }
    if (const auto* const rest_frequency = MemberObject(frequency_attributes, "rest_frequency")) {
        spectral.rest_frequency = MemberNumber(*rest_frequency, "data");
    }

    auto fit = FitSpectralAxis(spectral.channel_frequencies, reference_value, "frequency");
    spectral.reference_pixel = fit.reference_pixel;
    spectral.reference_value = fit.reference_value;
    spectral.increment = fit.increment;
    AppendDiagnostics(descriptor, std::move(fit.diagnostics));
    return spectral;
}

std::optional<TemporalCoordinate> DescribeTemporalCoordinate(const Store& store, std::vector<double> values) {
    if (values.empty()) {
        return std::nullopt;
    }

    TemporalCoordinate temporal;
    temporal.values = std::move(values);
    if (const auto& metadata = store.ReadNodeMetadata("time"); metadata) {
        if (const auto* const attributes = MemberObject(metadata.value(), "attributes")) {
            temporal.unit = AttributeString(*attributes, "units");
            temporal.scale = Upper(AttributeString(*attributes, "scale"));
            temporal.format = Upper(AttributeString(*attributes, "format"));
        }
    }
    return temporal;
}

// What discovery found, and the two things beyond the listing that deciding the match needs: whether
// anything in this store was malformed rather than merely unopenable, and which diagnostic says so;
// and whether anything in it is evidence of an image dataset at all.
//
// Here rather than on ImageDiscovery because it is this profile reasoning about its own store, and
// ImageDiscovery is what every profile answers with. A second profile would decide its own match
// from its own reasons; shaping that for it now would be guessing at a caller that does not exist.
struct Discovered {
    ImageDiscovery discovery;
    std::optional<std::size_t> first_malformation;
    // A malformed node says the store is broken, not what it is. An image variable this profile
    // recognised, openable or not, or a node bearing one of XRADIO's image names, says it is an
    // image dataset; so does a root that says so, which InspectImages asks.
    bool image_evidence = false;
};

// Which variables of this store are images this profile will open. Only InspectImages calls it:
// deciding whether the store matches means enumerating it, and the answer to both questions comes
// back together so that nothing has to enumerate twice or cache the result.
Result<Discovered> DiscoverImages(const Store& store) {
    const auto& inventory = store.Inventory();
    if (!inventory) {
        return inventory.error();
    }

    Discovered found;
    ImageDiscovery& result = found.discovery;
    for (const auto& entry : inventory.value()) {
        auto qualified = QualifyNode(store, entry);
        if (qualified.listed || KnownImageRank(entry.name) < kKnownImages) {
            found.image_evidence = true;
        }
        if (qualified.diagnostic) {
            if (qualified.malformed && !found.first_malformation) {
                found.first_malformation = result.diagnostics.size();
            }
            result.diagnostics.push_back(*qualified.diagnostic);
        }
        if (!qualified.listed) {
            continue;
        }
        std::vector<Diagnostic> said;
        if (qualified.diagnostic) {
            said.push_back(std::move(*qualified.diagnostic));
        }
        // Metadata the inventory has already parsed, and a listed node is always an array whose
        // metadata did parse: the two kinds that are not never qualify as listed.
        const auto& metadata = entry.array->value();
        // Only for an image that opens: the axes promise what OpenImage would report, and a
        // variable it refuses reports nothing -- an aperture-plane one would have only some of them.
        auto axes = qualified.openable ? DescribeAxes(store, metadata) : std::vector<AxisDescriptor>{};
        result.images.push_back(
            ImageEntry{entry.name, qualified.openable, std::move(said), ImageRoleOf(metadata), std::move(axes)});
    }

    const auto sort_images = [](std::vector<ImageEntry>& images) {
        std::sort(images.begin(), images.end(), [](const auto& left, const auto& right) {
            const int left_rank = KnownImageRank(left.id);
            const int right_rank = KnownImageRank(right.id);
            return left_rank == right_rank ? left.id < right.id : left_rank < right_rank;
        });
    };
    sort_images(result.images);
    const auto openable = std::find_if(result.images.begin(), result.images.end(),
                                       [](const ImageEntry& image) { return image.openable; });
    if (openable != result.images.end()) {
        result.default_image_id = openable->id;
    }
    return found;
}

}  // namespace

Result<SchemaInspection> InspectImages(const Store& store) {
    ProbeReport report(store, "image dataset");

    const auto& root_attributes = store.RootAttributes();
    auto found = DiscoverImages(store);
    if (!found) {
        return found.error();
    }
    auto& discovery = found.value().discovery;

    // What a refusal is reported as is the first diagnostic, and node enumeration order decides
    // which that is otherwise -- which is alphabetical, and has nothing to do with why the store was
    // refused. A store holding an aperture-plane variable and a malformed one would name the
    // aperture plane.
    if (!discovery.default_image_id && found.value().first_malformation) {
        auto& said = discovery.diagnostics;
        const auto reason = said.begin() + static_cast<std::ptrdiff_t>(*found.value().first_malformation);
        std::rotate(said.begin(), reason, reason + 1);
    }
    report.SetDiagnostics(discovery.diagnostics);

    // What was found goes back with what was decided from it, whichever way the decision went.
    const auto finish = [&](SchemaMatchKind kind) -> Result<SchemaInspection> {
        auto probe = report.Finish(kind, std::string(kVersion));
        if (!probe) {
            return probe.error();
        }
        return SchemaInspection{std::move(probe.value()), std::move(discovery)};
    };

    if (!discovery.default_image_id) {
        // Nothing here this profile will open, and which answer that is depends on why. A store of
        // complex or aperture-plane variables is well formed and simply not for this library. One
        // whose images disagree with their coordinates, or whose metadata will not parse, is a
        // store this profile recognised and found broken -- and reporting that as a non-match tells
        // whoever picked the file that nothing knew what it was, which is untrue and no use.
        //
        // But only a store with something of an image about it. A malformed node alone is not that:
        // a visibility dataset with a weights array that will not parse is Zarr of something else,
        // and calling it a broken image sent CARTA's file browser to list it as one.
        const bool image_dataset = found.value().image_evidence ||
                                   AttributeString(root_attributes, "type") == "image_dataset";
        return finish(found.value().first_malformation && image_dataset ? SchemaMatchKind::invalid
                                                                          : SchemaMatchKind::no_match);
    }

    // Once discovery found an openable image, validate the metadata needed by the image reader.
    const auto& first_image = *discovery.default_image_id;
    const auto& array_result = store.ReadArrayMetadata(first_image);
    if (report.RequireArrayMetadata(array_result, first_image)) {
        RequirePresentCoordinates(report);
        // Unconditionally: every XRADIO writes it, so a store without it is malformed rather than a
        // store with no direction to report.
        report.RequireCoordinateSystem(root_attributes);
    }
    return finish(report.ok() ? SchemaMatchKind::match : SchemaMatchKind::invalid);
}

Result<DescribedImage> DescribeImage(const Store& store, std::string_view image_id) {
    // Refused on its metadata before any value is read, so that an image this profile will not open
    // is reported as that rather than as whatever reading its coordinates ran into.
    if (auto qualified = RequireQualified(store, image_id); !qualified) {
        return qualified.error();
    }
    // And held to its own document, which is what its pixels are read from. Nothing below reads
    // them, so without this an image whose own document disagreed with the root's copy opened and
    // then failed every read. Its coordinates are held to theirs as they are read.
    if (const auto& verified = store.VerifyArray(image_id); !verified) {
        return verified.error();
    }
    auto values = ReadCoordinateValues(store);
    if (!values) {
        return values.error();
    }
    return DescribeImageFrom(store, image_id, values.value());
}

Result<DescribedImage> DescribeImageFrom(const Store& store, std::string_view image_id,
                                         const CoordinateValues& values) {
    // Asked rather than decided again. This used to classify the variable itself, on a weaker rule
    // than the one the listing was built with -- l and m rather than the whole axis set -- so the
    // two could disagree about what an image is.
    if (auto qualified = RequireQualified(store, image_id); !qualified) {
        return qualified.error();
    }
    // Parsed already, and the qualification just accepted it: the store hands back what it holds.
    const auto& image = store.ReadArrayMetadata(image_id).value();

    ImageDescriptor descriptor;
    descriptor.id = std::string(image_id);
    descriptor.stored_type = zarr_metadata::ParseDataType(image.data_type);
    descriptor.unit = AttributeString(image.attributes, "units");
    descriptor.image_role = ImageRoleOf(image);

    const auto& root_attrs = store.RootAttributes();
    descriptor.data_groups = FindDataGroups(root_attrs, image_id);
    descriptor.axes = DescribeAxes(store, image);

    // A direction axis is linear by construction, so its increment is reported even when the samples
    // are not evenly spaced; the fit says so in a diagnostic rather than withholding the value.
    auto direction = DescribeDirection(root_attrs, values.l, values.m, descriptor.diagnostics);
    if (!direction) {
        return direction.error();
    }
    descriptor.direction = std::move(direction.value());

    // Absent only when no samples were handed in, which reading them never does: every coordinate is
    // required, and the probe refused a dataset missing one. Same for the two below.
    if (auto spectral = DescribeSpectralCoordinate(store, values.frequency, descriptor); spectral) {
        descriptor.spectral = std::move(spectral);
    }

    if (!values.polarization.empty()) {
        PolarizationCoordinate pol;
        pol.labels = values.polarization;
        descriptor.polarization = std::move(pol);
    }

    if (auto temporal = DescribeTemporalCoordinate(store, values.time); temporal) {
        descriptor.temporal = std::move(temporal);
    }

    // Parse Observation & Telescope Metadata from the selected image's attributes.
    descriptor.observation = DescribeObservation(image);

    // Read from the metadata already in hand rather than asked of the store again, and assigned
    // unconditionally: every array ParseArrayMetadata accepted has a layout, so there is no failure
    // here to handle and no absence to represent.
    const auto layout = zarr_metadata::ParseStorageLayout(image);

    auto pixel_mask = DetermineFlag(store, image, image_id, descriptor.diagnostics);
    if (!pixel_mask) {
        return pixel_mask.error();
    }
    descriptor.pixel_mask_id = pixel_mask.value();
    descriptor.has_pixel_mask = !descriptor.pixel_mask_id.empty();
    // The flag is read only by a masked read, so an image whose flag disagreed with the root's copy
    // opened, read unmasked, and failed every masked read.
    ChunkGeometry flag_geometry;
    if (descriptor.has_pixel_mask) {
        const auto& verified = store.VerifyArray(descriptor.pixel_mask_id);
        if (!verified) {
            return verified.error();
        }
        flag_geometry = BuildChunkGeometry(descriptor, zarr_metadata::ParseStorageLayout(verified.value()));
    }

    auto geometry = BuildChunkGeometry(descriptor, layout);
    return DescribedImage{std::move(descriptor), std::move(geometry), std::move(flag_geometry)};
}

Result<std::vector<Beam>> ReadBeams(const Store& store, std::string_view image_id) {
    const auto& sky_meta = store.ReadNodeMetadata(image_id);
    if (!sky_meta) {
        return sky_meta.error();
    }
    std::string beam_array_name;
    if (const auto* const attributes = MemberObject(sky_meta.value(), "attributes")) {
        beam_array_name = AttributeString(*attributes, "beam_fit_params");
    }
    if (beam_array_name.empty()) {
        // No beam array associated
        return std::vector<Beam>{};
    }

    // Bound to the table's own document, so its layout and its unit are read from what its values
    // were decoded with rather than from the root's copy of it.
    const auto& table = store.ReadNumericArray(beam_array_name);
    if (!table) {
        return table.error();
    }

    // The image named a beam table, so what the table needs to be read is required rather than
    // optional. A label array that cannot be read -- absent, malformed, or written with a codec
    // this build does not carry -- used to be discarded here, and the parameter indices it did not
    // yield then returned an empty beam list: the answer for an image that has no beam at all,
    // which this one is not.
    auto parameter_labels = store.ReadStringArray1D("beam_params_label");
    if (!parameter_labels) {
        return parameter_labels.error();
    }
    return DescribeBeams(table.value(), parameter_labels.value());
}

}  // namespace carta::zarr::internal::xradio
