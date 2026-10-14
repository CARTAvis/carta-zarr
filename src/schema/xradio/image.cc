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
#include "schema_row.h"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <optional>
#include <utility>

#include <carta-zarr/descriptor.h>

namespace carta::zarr::internal::xradio {
namespace {

namespace zarr_metadata = ::carta::zarr::internal::zarr;

// Where an image is listed: XRADIO's own names in the row's order, then every other name after them.
std::size_t KnownImageRank(const SchemaRow& row, std::string_view image_id) {
    const auto& known = row.known_images;
    return static_cast<std::size_t>(std::find(known.begin(), known.end(), image_id) - known.begin());
}

bool IsKnownImage(const SchemaRow& row, std::string_view image_id) {
    return KnownImageRank(row, image_id) < row.known_images.size();
}

// XRADIO names a spectral frame that astropy also has by astropy's name, lowercased, and every other
// one by casacore's (CASACORE_TO_OBSERVER in XRADIO's image conventions). Uppercasing gives back
// casacore's name for all of them but the geocentric frame, which astropy calls GCRS and casacore GEO.
std::string SpectralSystemOf(const std::string& observer) {
    auto system = Upper(observer);
    return system == "GCRS" ? std::string{"GEO"} : system;
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

// XRADIO writes the role on the variable's own "type" attribute, lowercased; the same attribute
// spells "flag" for pixel masks, which are never images. A sky image's part is written two ways. From
// a list of tclean's images XRADIO 1.2 typed it there ("model", "residual"). From a dict of them, and
// in every schema 0.0.2 store, each is typed "sky" and its part is the name of the data group whose
// sky it is, so a sky image another group than "base" names as its sky, and "base" does not, plays
// that group's part.
//
// One rule for the listing and for describing, so that an entry and the image it opens as cannot
// disagree about which part the image plays.
std::string ImageRoleOf(const zarr_metadata::ArrayMetadata& image, std::string_view image_id,
                        const nlohmann::json& root_attributes) {
    const std::string role = AttributeString(image.attributes, "type");
    if (role == "flag") {
        return {};
    }
    const auto* const groups = role == "sky" ? MemberObject(root_attributes, "data_groups") : nullptr;
    if (groups == nullptr) {
        return role;
    }
    // The base group's sky is the sky, whatever other group shares it.
    if (const auto* const base = MemberObject(*groups, "base");
        base != nullptr && AttributeString(*base, "sky") == image_id) {
        return role;
    }
    for (const auto& [group_name, group] : groups->items()) {
        if (AttributeString(group, "sky") == image_id) {
            return group_name;
        }
    }
    return role;
}

std::optional<SpectralCoordinate> DescribeSpectralCoordinate(const Store& store,
                                                             const std::vector<double>& frequency_values,
                                                             const SchemaRow& row, ImageDescriptor& descriptor) {
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
    spectral.unit = CoordinateUnit(store, SkyCoordinate(AxisRole::spectral), row);
    if (row.spectral_frame == SpectralFrameSource::frame_attribute) {
        spectral.system = AttributeString(frequency_attributes, "frame");
    }
    // The channel a linear description is measured from, which is the first one unless the
    // reference frequency names another. What that measure carries -- for an XRADIO 1.2 store its
    // frame under `attrs`, and its value under `data` -- is read from the one lookup.
    double reference_value = spectral.channel_frequencies.front();
    if (const auto* const reference_frequency = MemberObject(frequency_attributes, "reference_frequency")) {
        const auto* const attributes = MemberObject(*reference_frequency, "attrs");
        if (attributes != nullptr && row.spectral_frame == SpectralFrameSource::reference_observer) {
            spectral.system = SpectralSystemOf(AttributeString(*attributes, "observer"));
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
Result<Discovered> DiscoverImages(const Store& store, const SchemaRow& row) {
    const auto& inventory = store.Inventory();
    if (!inventory) {
        return inventory.error();
    }

    Discovered found;
    ImageDiscovery& result = found.discovery;
    const auto& root_attributes = store.RootAttributes();
    for (const auto& entry : inventory.value()) {
        auto qualified = QualifyNode(store, entry);
        if (qualified.listed || IsKnownImage(row, entry.name)) {
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
        auto axes = qualified.openable ? DescribeAxes(store, metadata, row) : std::vector<AxisDescriptor>{};
        result.images.push_back(ImageEntry{entry.name, qualified.openable, std::move(said),
                                           ImageRoleOf(metadata, entry.name, root_attributes), std::move(axes)});
    }

    const auto sort_images = [&row](std::vector<ImageEntry>& images) {
        std::sort(images.begin(), images.end(), [&row](const auto& left, const auto& right) {
            const auto left_rank = KnownImageRank(row, left.id);
            const auto right_rank = KnownImageRank(row, right.id);
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
    // A store naming a version no row reads is still enumerated, by the newest row, to tell an image
    // dataset of that version -- refused for its version -- from Zarr of something else that carries
    // the attribute too, as XRADIO's measurement sets do.
    auto chosen = ChooseSchemaRow(store);
    const auto& row = chosen.row != nullptr ? *chosen.row : SchemaRows().front();
    auto found = DiscoverImages(store, row);
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
        auto probe = report.Finish(kind, chosen.version);
        if (!probe) {
            return probe.error();
        }
        return SchemaInspection{std::move(probe.value()), std::move(discovery)};
    };

    if (chosen.row == nullptr) {
        // Refused for its version ahead of anything else wrong with it, which is reported by rules
        // the store may not follow. Something that is not an image dataset is not refused at all.
        if (!image_dataset) {
            report.SetDiagnostics(discovery.diagnostics);
            return finish(SchemaMatchKind::no_match);
        }
        auto said = std::move(chosen.diagnostics);
        said.insert(said.end(), discovery.diagnostics.begin(), discovery.diagnostics.end());
        report.SetDiagnostics(std::move(said));
        return finish(SchemaMatchKind::invalid);
    }

    if (!discovery.default_image_id) {
        report.SetDiagnostics(discovery.diagnostics);
        // Nothing here this profile will open, and which answer that is depends on why. A store of
        // complex or aperture-plane variables is well formed and simply not for this library. One
        // whose images disagree with their coordinates, or whose metadata will not parse, is a
        // store this profile recognised and found broken -- and reporting that as a non-match tells
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
        RequirePresentCoordinates(report);
        // Unconditionally: every XRADIO writes it, so a store without it is malformed rather than a
        // store with no direction to report.
        report.RequireCoordinateSystem(root_attributes);
    }
    if (!report.ok()) {
        return finish(SchemaMatchKind::invalid);
    }
    // Which rules read a dataset that matched, after what discovery said: a refusal is reported as
    // its first diagnostic, and which row would have read the store is not the reason for one.
    auto said = report.diagnostics();
    said.insert(said.end(), chosen.diagnostics.begin(), chosen.diagnostics.end());
    report.SetDiagnostics(std::move(said));
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
    auto values = ReadCoordinateValues(store);
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
    const auto chosen = ChooseSchemaRow(store);
    if (chosen.row == nullptr) {
        const auto& refusal = chosen.diagnostics.front();
        return Error{ErrorCode::unsupported_schema_version, refusal.message, refusal.node_path};
    }
    const auto& row = *chosen.row;
    // Parsed already, and the qualification just accepted it: the store hands back what it holds.
    const auto& image = store.ReadArrayMetadata(image_id).value();

    ImageDescriptor descriptor;
    descriptor.id = std::string(image_id);
    descriptor.stored_type = zarr_metadata::ParseDataType(image.data_type);
    descriptor.unit = AttributeString(image.attributes, "units");
    const auto& root_attrs = store.RootAttributes();
    descriptor.image_role = ImageRoleOf(image, image_id, root_attrs);
    descriptor.data_groups = FindDataGroups(root_attrs, image_id);
    descriptor.axes = DescribeAxes(store, image, row);

    // A direction axis is linear by construction, so its increment is reported even when the samples
    // are not evenly spaced; the fit says so in a diagnostic rather than withholding the value.
    auto direction = DescribeDirection(root_attrs, values.l, values.m, descriptor.diagnostics);
    if (!direction) {
        return direction.error();
    }
    descriptor.direction = std::move(direction.value());

    // Absent only when no samples were handed in, which reading them never does: every coordinate is
    // required, and the probe refused a dataset missing one. Same for the two below.
    if (auto spectral = DescribeSpectralCoordinate(store, values.frequency, row, descriptor); spectral) {
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
    descriptor.observation = DescribeObservation(image, row);

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
    // this build does not carry -- fails the read: an empty beam list is the answer for an image
    // that has no beam at all, which this one is not.
    auto parameter_labels = store.ReadStringArray1D("beam_params_label");
    if (!parameter_labels) {
        return parameter_labels.error();
    }
    return DescribeBeams(table.value(), parameter_labels.value());
}

}  // namespace carta::zarr::internal::xradio
