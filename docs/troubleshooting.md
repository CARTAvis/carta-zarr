# Troubleshooting

Check every `Result` before accessing its value. Report the error code name,
message and node path together:

```cpp
const auto& error = result.error();  // Only after checking !result.
std::cerr << ErrorCodeName(error.code) << ": " << error.message
          << " [" << error.node_path << "]\n";
```

| Symptom or code | What to check |
|---|---|
| `not_found`, `not_zarr` | The path exists and names a Zarr group with root metadata |
| `unsupported_transport` | Use a local filesystem path; remote stores are not implemented |
| `unsupported_zarr_version` | The dataset uses Zarr v3 |
| `unsupported_schema` | The data follows the XRADIO image profile |
| `unsupported_schema_version` | The image dataset declares schema 0.0.2 or has no declared version; other declared versions are refused |
| `unversioned_schema` diagnostic | The reader used older XRADIO rules: no version was declared, or declared 0.0.2 lacks its required version-specific fields |
| An entry has `openable == false` | Inspect its diagnostics; complex, aperture-plane and unsupported boolean candidates cannot be opened |
| `invalid_metadata` during `OpenImage` | The image, coordinates and declared flag agree; consolidated metadata matches each array's own document |
| `invalid_slice` | Ranges fit the descriptor axes and use valid counts and strides |
| `buffer_too_small` during a read | Destination size is an element count and covers the product of selected counts |
| `buffer_too_small` creating read-ahead | Cache capacity can retain two runs of decoded chunks |
| `io_error`, `decode_error`, `unsupported_codec` | Inspect the exact node and message for permissions, inaccessible chunks or an encoding problem |
| `cancelled` | A cancellation callback, progress sink or deadline stopped the operation |

A diagnostic is different from an error: the dataset or descriptor is still
returned. For example, an optional observation field may be unavailable while
the pixels remain readable. `DiagnosticCodeName` provides the diagnostic's name.
The probe's diagnostics can also explain `kind == invalid`, even though the
probe call itself succeeded. A 0.0.2 version string with `unversioned_schema`
means a fallback occurred; consult [schema compatibility](guides/compatibility.md)
before relying on version-specific fields such as `sub_type` or `channel_width`.

## Surprising pixel or statistic values

- Flag application is on by default. A true flag becomes NaN; disable
  `apply_pixel_mask` only when you intend to inspect the stored image values.
- Statistics count finite pixels separately from excluded or non-finite ones.
  Do not derive a mean with a zero denominator.
- Unfinished reduction blocks are estimates over pixels visited so far. Their
  final replacement has `complete == true`.
- Cube histogram bins are approximate after re-aggregation; their exact
  distribution can depend on worker count.

## Stale data after changing files

Treat an open dataset as a read-only snapshot. Metadata and decoded chunks are
cached and are not revalidated on every read. Do not mutate the dataset while
using its handles. After replacing data, reopen with fresh resources (a new
Context and Dataset) rather than expecting an existing session to refresh.

For bug reports, include the library version, operation, selected ranges,
error or diagnostics and a small reproducible dataset if available.
