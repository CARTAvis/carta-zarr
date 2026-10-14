# Read XRADIO images from C++

carta-zarr is a read-only C++17 library for XRADIO images stored as Zarr v3.
Open a dataset, choose an image, and read pixels or compute statistics without
exposing TensorStore to your application.

```{warning}
This is **development documentation** for an unreleased library. Before 1.0,
minor versions may change both the source API and ABI. The intended integration
with carta-backend is also under development.
```

**Start here:** [Install the library](installation.md), then follow the
[complete first-read example](quickstart.md). If you already have it installed,
learn how [datasets and images](concepts.md) relate.

## What can I read?

| Supported | Scope |
|---|---|
| Schema | XRADIO image schema 0.0.2 and older unversioned XRADIO 1.2 stores |
| Storage | Local filesystem Zarr v3 groups, including consolidated metadata |
| Images | Real-valued sky-plane variables with time, frequency, polarization, l and m axes |
| Output | Pixels converted to float32; reductions accumulated in double |
| Flags | Applied by default; flagged pixels become NaN |

Writing, remote object stores, complex-valued images and aperture-plane images
are outside the current supported set. Discovery can list unsupported image
candidates with diagnostics explaining why they cannot be opened.
See [schema compatibility](guides/compatibility.md) for version selection,
fallback diagnostics and upgrading an application to carta-zarr 0.2.

## Choose your next task

- [List images and read a channel or ROI](guides/reading.md).
- [Reduce regions and compute histograms](guides/statistics.md).
- [Control memory, cancellation and animation read-ahead](guides/resources.md).
- [Build and run the C++ examples](examples.md).
- [Diagnose a failed open or read](troubleshooting.md).
- [Prepare an XRADIO dataset](guides/datasets.md).
- [Tune storage layouts](storage-tuning.md).
- [Understand the architecture](architecture.md).

```{toctree}
:hidden:
:caption: Getting started

installation
quickstart
concepts
examples
```

```{toctree}
:hidden:
:caption: How-to guides

guides/reading
guides/statistics
guides/resources
guides/datasets
guides/compatibility
troubleshooting
storage-tuning
```

```{toctree}
:hidden:
:caption: Reference

api/index
architecture
glossary
changelog
contributing-docs
```
