# Troubleshooting

Check every `Result` before accessing its value, as the
[core concepts](concepts.md#rules-that-apply-to-every-operation) require. Report
the error code name, message and node path together:

```cpp
const auto& error = result.error();  // Only after checking !result.
std::cerr << ErrorCodeName(error.code) << ": " << error.message;
if (!error.node_path.empty()) {
    std::cerr << " [" << error.node_path << ']';
}
std::cerr << '\n';
```

## Error codes

| Symptom or code | What to check |
|---|---|
| `not_found`, `not_zarr` | Check that the path exists and names a Zarr group with root metadata |
| `unsupported_transport` | Use a local filesystem path; remote stores are not implemented |
| `unsupported_zarr_version` | Check that the dataset is Zarr v3, not v2 |
| `unsupported_schema` | Check that the data is an XRADIO image dataset |
| `unsupported_schema_version` | Check the dataset's root `schema_version`: only 0.0.2, or none, can be read |
| An entry has `openable == false` | Inspect its diagnostics; complex, aperture-plane and unsupported boolean candidates cannot be opened |
| `invalid_metadata` during `OpenImage` | Check that the image, its coordinates and its declared flag agree, and that consolidated metadata matches each array's own document |
| `invalid_slice` | Check that the ranges fit the descriptor axes and use valid counts and strides |
| `buffer_too_small` during a read | Check that the destination's size, in elements, covers the product of the selected counts |
| `buffer_too_small` creating read-ahead | Give each cache room for two runs of decoded chunks |
| `io_error`, `decode_error`, `unsupported_codec` | Inspect the node and message for permissions, inaccessible chunks or an encoding problem |
| `cancelled` | Expected when a cancellation callback, progress callback or deadline stopped the operation |

## Diagnostics

A diagnostic is different from an error: the dataset or descriptor is still
returned. For example, an optional observation field may be unavailable while
the pixels remain readable. `DiagnosticCodeName` provides the diagnostic's name.
The probe's diagnostics can also explain `kind == invalid`, even though the
probe call itself succeeded. An `unversioned_schema` diagnostic means the
dataset was read by the older XRADIO rules: it declares no version, or declares
0.0.2 without the fields 0.0.2 requires. Consult
[schema versions](reference/schema-versions.md) before relying on
version-specific fields such as `sub_type` or `channel_width`.

## Surprising pixel or statistic values

- Flag application is on by default. A true flag becomes NaN; disable
  `apply_pixel_mask` only when you intend to inspect the stored image values.
- Statistics count finite pixels separately from flagged or non-finite ones.
  When `num_pixels` is zero there is no mean; check it before dividing.
- Unfinished reduction blocks are estimates over pixels visited so far. Their
  final replacement has `complete == true`.
- Cube histogram bins are approximate after re-aggregation; their exact
  distribution can depend on worker count.

## Stale data after changing files

Treat an open dataset as a read-only snapshot. Metadata and decoded chunks are
cached and are not revalidated on every read. Do not mutate the dataset while
using its handles. After replacing data, reopen with fresh resources (a new
Context and Dataset) rather than expecting an existing session to refresh.

## Reporting a problem

Include the library version (the `CartaZarr_VERSION` your CMake project found,
or the commit you built), the operation, the selected ranges, the error or
diagnostics, and a small dataset that reproduces it if you can share one.
