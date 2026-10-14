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
See [schema versions and reading rules](reference/schema-versions.md) for version
selection and fallback diagnostics, or [upgrade to 0.2](upgrading/0.2.md) from 0.1.

## Choose your next task

**Read and compute**

- [List images and read a channel or region](guides/reading.md).
- [Reduce regions and compute histograms](guides/statistics.md).
- [Control memory, cancellation and animation read-ahead](guides/resources.md).
- [Build and run the C++ examples](examples.md).

**Data and storage**

- [Prepare an XRADIO dataset](guides/datasets.md).
- [Tune storage layouts](storage-tuning.md).

**When something goes wrong, or changes**

- [Diagnose a failed open or read](troubleshooting.md).
- [Upgrade from 0.1 to 0.2](upgrading/0.2.md).
- [Understand the architecture](architecture.md).

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
```

```{toctree}
:hidden:
:caption: Troubleshooting

troubleshooting
```

```{toctree}
:hidden:
:caption: Reference

api/index
reference/schema-versions
glossary
changelog
```

```{toctree}
:hidden:
:caption: Background

architecture
```

```{toctree}
:hidden:
:caption: Upgrading

upgrading/0.2
```

```{toctree}
:hidden:
:caption: Contributing

contributing-docs
```
