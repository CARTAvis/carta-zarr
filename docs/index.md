# Read XRADIO images from C++

carta-zarr is a read-only C++17 library for radio-astronomy images stored as
[Zarr](https://zarr.dev/) v3 in the [XRADIO](https://github.com/casangi/xradio)
layout, the image schema defined by the XRADIO Python package. Open a dataset, choose an
image, and read pixels or compute statistics without exposing TensorStore, the
storage engine underneath, to your application. Terms are defined in the
[glossary](glossary.md).

```{warning}
This documentation follows the `dev` branch of a library before its 1.0 release.
Before 1.0, minor versions may change both the source API and ABI. The intended integration
with carta-backend is also under development.
```

**Start here:** [Install the library](installation.md), then follow the
[complete first-read example](quickstart.md). If you already have it installed,
learn how [datasets and images](concepts.md) relate.

## What can I read?

| Area | Current support |
|---|---|
| Schema | XRADIO image schema 0.0.2, as XRADIO 1.2.5 and later write it |
| Storage | Local filesystem Zarr v3 groups, including consolidated metadata |
| Images | Real-valued variables with time, frequency, polarization, and either direction-cosine (`l`, `m`) or aperture-plane (`u`, `v`) axes |
| Output | Pixels converted to float32; reductions accumulated in double |
| Flags | Applied by default; flagged pixels become NaN |

Writing, remote object stores and complex-valued images are outside the current
supported set. A real-valued image of the aperture plane (`u`, `v`) opens, with an
aperture coordinate in place of a direction. Listing the images in a dataset can show
unsupported candidates, with diagnostics explaining why each cannot be opened.
See [schema versions and reading rules](reference/schema-versions.md) for version
selection and fallback diagnostics. Changes since 0.1 are listed in the [changelog](changelog.md).

## Choose your next task

**Read and compute**

- [List images and read a channel or region](guides/reading.md).
- [Reduce regions and compute histograms](guides/statistics.md).
- [Control memory, cancellation and read-ahead](guides/resources.md).
- [Build and run the C++ examples](examples.md).

**Data and storage**

- [Prepare an XRADIO dataset](guides/datasets.md).
- [Tune storage layouts](storage-tuning.md).
- [Measure your own storage](guides/measure-storage.md) with the benchmark tools.

**When something goes wrong**

- [Diagnose a failed open or read](troubleshooting.md).

**Background**

- [Understand the architecture](architecture.md): what the library does and
  what it leaves to the application.

```{toctree}
:hidden:
:caption: Get started

installation
quickstart
concepts
examples
```

```{toctree}
:hidden:
:caption: Read and compute

guides/reading
guides/statistics
guides/resources
```

```{toctree}
:hidden:
:caption: Data and storage

guides/datasets
storage-tuning
guides/measure-storage
measurements/chunk-size
```

```{toctree}
:hidden:
:caption: Reference

api/index
reference/schema-versions
reference/zarr-bench
architecture
glossary
changelog
```

```{toctree}
:hidden:
:caption: Help

troubleshooting
```

```{toctree}
:hidden:
:caption: Project

contributing-docs
```
