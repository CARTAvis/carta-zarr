/* This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
   Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
   Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
   SPDX-License-Identifier: GPL-3.0-or-later
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
#include "schema_version.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <tuple>
#include <utility>

#include <carta-zarr/descriptor.h>

namespace carta::zarr::internal::xradio {
namespace {

namespace zarr_metadata = ::carta::zarr::internal::zarr;

// Where an image no data group has as its own is listed: XRADIO's own names in their order, then every
// other name after them.
std::size_t KnownImageRank(std::string_view image_id) {
    return static_cast<std::size_t>(std::find(kKnownImages.begin(), kKnownImages.end(), image_id) -
                                    kKnownImages.begin());
}

bool IsKnownImage(std::string_view image_id) {
    return KnownImageRank(image_id) < kKnownImages.size();
}

// Where each data group's own image is listed: the base group's first, which is the store's main
// image, then every other group's by group name. A group's own image is its `sky`, or its
// `aperture` when it has no sky. The group order is by name because the JSON object holding the
// groups keeps no order this library reads; XRADIO writes `base` first.
std::map<std::string, std::size_t> GroupImageRanks(const nlohmann::json& root_attributes) {
    std::map<std::string, std::size_t> ranks;
    const auto* const groups = MemberObject(root_attributes, "data_groups");
    if (groups == nullptr) {
        return ranks;
    }
    const auto rank_of = [&](const nlohmann::json& group, std::size_t rank) {
        std::string image = AttributeString(group, "sky");
        if (image.empty()) {
            image = AttributeString(group, "aperture");
        }
        auto named = NormalizeNodeName(image);
        if (!image.empty() && named) {
            ranks.emplace(std::move(named).value(), rank);
        }
    };
    if (const auto* const base = MemberObject(*groups, "base")) {
        rank_of(*base, 0);
    }
    std::size_t rank = 1;
    for (const auto& [name, group] : groups->items()) {
        if (name != "base" && group.is_object()) {
            rank_of(group, rank++);
        }
    }
    return ranks;
}

// The dataset carries a well-formed array for every coordinate of each plane it has an image of. Both
// XRADIO readers write all of a plane's unconditionally (xds_from_casacore.py and xds_from_fits.py
// both assign coords["time"]), so a missing coordinate means the store is malformed.
void RequirePresentCoordinates(ProbeReport& report, bool sky, bool aperture) {
    for (const auto& coordinate : kCoordinates) {
        if ((sky && OnPlane(coordinate, Plane::sky)) || (aperture && OnPlane(coordinate, Plane::aperture))) {
            report.RequireCoordinate(coordinate);
        }
    }
}

// The aperture plane's linear coordinate, from u's and v's samples and the `crval` XRADIO writes
// beside them. A degenerate axis leaves its reference pixel and increment at zero, as a direction
// axis does, and says so.
ApertureCoordinate DescribeAperture(const Store& store, const std::vector<double>& u_values,
                                    const std::vector<double>& v_values, std::vector<Diagnostic>& diagnostics) {
    ApertureCoordinate aperture;
    aperture.unit = CoordinateUnit(store, PlaneCoordinate(AxisRole::spatial_x, Plane::aperture));
    const auto set_axis = [&](const std::vector<double>& values, std::size_t axis, std::string_view name) {
        std::optional<double> reference;
        if (const auto& metadata = store.ReadArrayMetadata(name); metadata) {
            reference = MemberNumber(metadata.value().attributes, "crval");
        }
        if (!reference) {
            diagnostics.push_back(Diagnostic{DiagnosticCode::invalid_metadata,
                                             "Aperture coordinate has no crval; its reference value is taken as 0",
                                             std::string(name)});
        }
        auto fit = FitApertureAxis(values, reference, name);
        aperture.reference_value.at(axis) = reference.value_or(0.0);
        if (fit.increment) {
            aperture.increment.at(axis) = *fit.increment;
        }
        if (fit.reference_pixel) {
            aperture.reference_pixel.at(axis) = *fit.reference_pixel;
        }
        diagnostics.insert(diagnostics.end(), std::make_move_iterator(fit.diagnostics.begin()),
                           std::make_move_iterator(fit.diagnostics.end()));
    };
    set_axis(u_values, 0, "u");
    set_axis(v_values, 1, "v");
    return aperture;
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

std::optional<SpectralCoordinate> DescribeSpectralCoordinate(const Store& store,
                                                             const std::vector<double>& frequency_values,
                                                             ImageDescriptor& descriptor) {
    if (frequency_values.empty()) {
        return std::nullopt;
    }

    SpectralCoordinate spectral;
    spectral.values = frequency_values;
    nlohmann::json frequency_attributes = nlohmann::json::object();
    if (const auto& frequency_metadata = store.ReadNodeMetadata("frequency"); frequency_metadata) {
        if (const auto* const attributes = MemberObject(frequency_metadata.value(), "attributes")) {
            frequency_attributes = *attributes;
        }
    }

    spectral.unit = CoordinateUnit(store, SkyCoordinate(AxisRole::spectral));
    spectral.reference_frame = AttributeString(frequency_attributes, "frame");
    // The channel a linear description is measured from, which is the first one unless the
    // reference frequency names another.
    double reference_value = spectral.values.front();
    if (const auto* const reference_frequency = MemberObject(frequency_attributes, "reference_frequency")) {
        if (const auto value = MemberNumber(*reference_frequency, "data")) {
            reference_value = *value;
        }
    }
    if (const auto* const rest_frequency = MemberObject(frequency_attributes, "rest_frequency")) {
        spectral.rest_frequency = MemberNumber(*rest_frequency, "data");
    }

    auto fit = FitSpectralAxis(spectral.values, reference_value, "frequency");
    spectral.reference_pixel = fit.reference_pixel;
    spectral.reference_value = fit.reference_value;
    spectral.increment = fit.increment;
    AppendDiagnostics(descriptor, std::move(fit.diagnostics));

    // An axis of one channel has no spacing to fit an increment from. Schema 0.0.2 states the width
    // of every channel, and for the one channel that is the increment, measured from it. An axis of
    // more channels keeps what its channels say, whatever the width says.
    if (!spectral.increment && spectral.values.size() == 1) {
        if (const auto* const width = MemberObject(frequency_attributes, "channel_width")) {
            const auto* const attributes = MemberObject(*width, "attrs");
            const std::string unit = attributes == nullptr ? std::string{} : AttributeString(*attributes, "units");
            const auto value = MemberNumber(*width, "data");
            if (value && std::isfinite(*value) && *value != 0.0 && (unit.empty() || unit == spectral.unit)) {
                spectral.reference_pixel = 1.0;
                spectral.reference_value = spectral.values.front();
                spectral.increment = *value;
            }
        }
    }
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
    // recognized, openable or not, or a node bearing one of XRADIO's image names, says it is an
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
        if (qualified.listed || IsKnownImage(entry.name)) {
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
        auto axes = qualified.openable ? DescribeAxes(metadata) : std::vector<AxisDescriptor>{};
        result.images.push_back(ImageEntry{entry.name, qualified.openable, std::move(said), std::move(axes)});
    }

    // The data groups' own images first, as the store ranks them; then the rest, XRADIO's names in
    // their order and every other name after them. The first that opens is the default.
    const auto group_ranks = GroupImageRanks(store.RootAttributes());
    const auto rank_of = [&](const std::string& id) {
        const auto found = group_ranks.find(id);
        return std::tuple<std::size_t, std::size_t, const std::string&>(
            found == group_ranks.end() ? std::numeric_limits<std::size_t>::max() : found->second, KnownImageRank(id),
            id);
    };
    std::sort(result.images.begin(), result.images.end(),
              [&](const ImageEntry& left, const ImageEntry& right) { return rank_of(left.id) < rank_of(right.id); });
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
    // A store that is not read is still enumerated, to tell an image dataset -- refused for what it
    // names -- from Zarr of something else that names a version too, as XRADIO's measurement sets do,
    // or none.
    auto checked = CheckSchemaVersion(store);
    auto found = DiscoverImages(store);
    if (!found) {
        return found.error();
    }
    auto& discovery = found.value().discovery;
    // A malformed node alone is not evidence of an image dataset: a visibility dataset with a weights
    // array that will not parse is Zarr of something else, and calling it a broken image sent CARTA's
    // file browser to list it as one.
    const bool image_dataset =
        found.value().image_evidence || AttributeString(root_attributes, "type") == "image_dataset";

    // What a refusal is reported as is the first diagnostic, and node enumeration order decides
    // which that is otherwise -- which is alphabetical, and has nothing to do with why the store was
    // refused. A store holding an aperture-plane variable and a malformed one would name the
    // aperture plane.
    if (!discovery.default_image_id && found.value().first_malformation) {
        auto& said = discovery.diagnostics;
        const auto reason = said.begin() + static_cast<std::ptrdiff_t>(*found.value().first_malformation);
        std::rotate(said.begin(), reason, reason + 1);
    }
    // What was found goes back with what was decided from it, whichever way the decision went.
    const auto finish = [&](SchemaMatchKind kind) -> Result<SchemaInspection> {
        auto probe = report.Finish(kind, checked.version);
        if (!probe) {
            return probe.error();
        }
        return SchemaInspection{std::move(probe.value()), std::move(discovery)};
    };

    if (checked.refusal) {
        // Refused for its version ahead of anything else wrong with it, which is reported by rules
        // the store may not follow. Something that is not an image dataset is not refused at all.
        if (!image_dataset) {
            report.SetDiagnostics(discovery.diagnostics);
            return finish(SchemaMatchKind::no_match);
        }
        std::vector<Diagnostic> said{std::move(*checked.refusal)};
        said.insert(said.end(), discovery.diagnostics.begin(), discovery.diagnostics.end());
        report.SetDiagnostics(std::move(said));
        return finish(SchemaMatchKind::invalid);
    }

    if (!discovery.default_image_id) {
        report.SetDiagnostics(discovery.diagnostics);
        // Nothing here this profile will open, and which answer that is depends on why. A store of
        // complex images, on either plane, is well formed and simply not for this library. One
        // whose images disagree with their coordinates, or whose metadata will not parse, is a
        // store this profile recognized and found broken -- and reporting that as a non-match tells
        // whoever picked the file that nothing knew what it was, which is untrue and no use. But only
        // a store with something of an image about it.
        return finish(found.value().first_malformation && image_dataset ? SchemaMatchKind::invalid
                                                                        : SchemaMatchKind::no_match);
    }

    // Once discovery found an openable image, validate the metadata needed by the image reader.
    report.SetDiagnostics(discovery.diagnostics);
    const auto& first_image = *discovery.default_image_id;
    const auto& array_result = store.ReadArrayMetadata(first_image);
    if (report.RequireArrayMetadata(array_result, first_image)) {
        // The planes of every image that opens, whichever opens first: each needs its own coordinates.
        bool sky = false;
        bool aperture = false;
        for (const auto& entry : discovery.images) {
            if (!entry.openable) {
                continue;
            }
            const auto& metadata = store.ReadArrayMetadata(entry.id);
            const auto plane = metadata ? PlaneOf(metadata.value()) : std::nullopt;
            sky = sky || plane == Plane::sky;
            aperture = aperture || plane == Plane::aperture;
        }
        RequirePresentCoordinates(report, sky, aperture);
        // Whenever an image of the sky opens: every XRADIO writes it then, so a store without it is
        // malformed rather than a store with no direction to report. An aperture-plane image has no
        // direction, and XRADIO writes none for a dataset of them alone.
        if (sky) {
            report.RequireCoordinateSystem(root_attributes);
        }
    }
    if (!report.ok()) {
        return finish(SchemaMatchKind::invalid);
    }
    return finish(SchemaMatchKind::match);
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
    // Qualified, so on a plane; the store's listing holds its metadata.
    const auto plane = PlaneOf(store.ReadArrayMetadata(image_id).value()).value_or(Plane::sky);
    auto values = ReadCoordinateValues(store, plane);
    if (!values) {
        return values.error();
    }
    return DescribeImageFrom(store, image_id, values.value());
}

Result<DescribedImage> DescribeImageFrom(const Store& store, std::string_view image_id,
                                         const CoordinateValues& values) {
    // Asked rather than decided again, so that describing an image and the listing cannot disagree
    // about what an image is.
    if (auto qualified = RequireQualified(store, image_id); !qualified) {
        return qualified.error();
    }
    if (const auto checked = CheckSchemaVersion(store); checked.refusal) {
        const auto& refusal = *checked.refusal;
        const auto code = refusal.code == DiagnosticCode::unsupported_schema_version
                            ? ErrorCode::unsupported_schema_version
                            : ErrorCode::invalid_schema;
        return Error{code, refusal.message, refusal.node_path};
    }
    // Parsed already, and the qualification just accepted it: the store hands back what it holds.
    const auto& image = store.ReadArrayMetadata(image_id).value();

    ImageDescriptor descriptor;
    descriptor.id = std::string(image_id);
    descriptor.stored_type = zarr_metadata::ParseDataType(image.data_type);
    descriptor.unit = AttributeString(image.attributes, "units");
    const auto& root_attrs = store.RootAttributes();
    descriptor.sub_type = AttributeString(image.attributes, "sub_type");
    descriptor.data_groups = FindDataGroups(root_attrs, image_id);
    descriptor.axes = DescribeAxes(image);

    // A direction axis is linear by construction, so its increment is reported even when the samples
    // are not evenly spaced; the fit says so in a diagnostic rather than withholding the value. So is
    // an aperture-plane axis, and an image is on one plane or the other (ADR 0024).
    if (PlaneOf(image).value_or(Plane::sky) == Plane::aperture) {
        descriptor.aperture = DescribeAperture(store, values.u, values.v, descriptor.diagnostics);
    } else {
        auto direction = DescribeDirection(root_attrs, values.l, values.m, descriptor.diagnostics);
        if (!direction) {
            return direction.error();
        }
        descriptor.direction = std::move(direction.value());
    }

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
    descriptor.user_keywords = DescribeUserKeywords(root_attrs, image, image_id, descriptor.diagnostics);

    // Read from the metadata already in hand rather than asked of the store again, and assigned
    // unconditionally: every array ParseArrayMetadata accepted has a layout, so there is no failure
    // here to handle and no absence to represent.
    const auto layout = zarr_metadata::ParseStorageLayout(image);

    auto pixel_mask = DetermineFlag(store, image, image_id);
    if (!pixel_mask) {
        return pixel_mask.error();
    }
    descriptor.pixel_mask_id = pixel_mask.value();
    // The flag is read only by a masked read, so an image whose flag disagreed with the root's copy
    // opened, read unmasked, and failed every masked read.
    ChunkGeometry flag_geometry;
    if (!descriptor.pixel_mask_id.empty()) {
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
    // A beam is the sky's: the aperture plane has none. XRADIO 1.2.5 and 1.3.0 name a
    // BEAM_FIT_PARAMS_APERTURE on an aperture-plane image all the same, and write no such array.
    if (const auto& image = store.ReadArrayMetadata(image_id); image && PlaneOf(image.value()) == Plane::aperture) {
        return std::vector<Beam>{};
    }
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
    // this build does not carry -- fails the read: an empty beam list is the answer for an image
    // that has no beam at all, which this one is not.
    auto parameter_labels = store.ReadStringArray1D("beam_params_label");
    if (!parameter_labels) {
        return parameter_labels.error();
    }
    return DescribeBeams(table.value(), parameter_labels.value());
}

}  // namespace carta::zarr::internal::xradio
