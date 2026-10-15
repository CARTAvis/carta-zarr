# Contributing to the documentation

## Build locally

### Prerequisites

You need Git, Doxygen 1.18 and, for the recommended build, uv. CI pins Doxygen
1.18 because other releases may change the XML that `docs/conf.py`
post-processes.

- **macOS:** `brew install git uv doxygen`.
- **Ubuntu:** `sudo apt-get install git`, then install uv using its
  [official instructions](https://docs.astral.sh/uv/getting-started/installation/).
  The packaged Doxygen is older: download `doxygen-1.18.0.linux.bin.tar.gz` from
  [Doxygen's releases](https://github.com/doxygen/doxygen/releases/tag/Release_1_18_0)
  and put its `bin` directory on `PATH`.

Check that `doxygen --version` prints 1.18.0. If Homebrew installed another
version, use the release download as on Ubuntu.

### With uv (recommended)

uv downloads the Python it needs, so the documentation packages stay out of
your system Python.

Run the commands below from the root of the checkout you want to preview.

Build and serve:

```sh
uv run --no-project --python 3.12 \
  --with-requirements docs/requirements.txt \
  python -m sphinx -b html -W --keep-going docs build-docs/html

uv run --no-project --python 3.12 \
  python -m http.server 8000 --bind 127.0.0.1 --directory build-docs/html
```

Open `http://localhost:8000` in a browser **on that computer**. Keep the terminal
running; press Ctrl+C to stop serving. After updating the checkout with `git pull`,
rerun the build command and refresh the browser. uv caches the tools between
invocations; no activation step or `uv init` is needed.

### With venv and pip

With the [prerequisites](#prerequisites) other than uv, run from the
repository root with Python 3.11 or newer, which Sphinx requires (CI uses 3.12):

```sh
python3 -m venv build-docs/venv
build-docs/venv/bin/python -m pip install -r docs/requirements.txt
build-docs/venv/bin/python -m sphinx -b html -W --keep-going docs build-docs/html
python3 -m http.server 8000 --bind 127.0.0.1 --directory build-docs/html
```

Open `http://localhost:8000`. The docs build runs Doxygen automatically and does
not configure or compile TensorStore. Build output and virtual environments
stay under the ignored `build-docs/` directory.

## Sources and checks

- Edit Markdown under `docs/` for tutorials and guides.
- Edit public header comments for API contracts. A small Doxygen input filter
  recognises the existing `//` prose without changing the installed headers.
- Glossary and changelog pages read the root files rather than maintaining
  copies. The storage-tuning page is edited in place, like any guide. During
  the site build, repository-only links in the changelog and storage guide
  point to GitHub.
- The quickstart and advanced guides include the [runnable examples](examples.md)
  directly. The package-consumer test builds all four against the installed
  package, runs them on committed unversioned, 0.0.2 and continuum fixtures,
  and checks missing-argument failures; its assertions are in
  `tests/run_package_consumer.cmake`. Run it after building the library with
  `ctest --test-dir build -R '^carta-zarr-package-consumer$' --output-on-failure`.
  Keep example changes in step with the public API.
- Keep diagrams as Mermaid text next to their explanation. CI renders them to SVG
  while building, so a syntax error fails the build. Locally, `npm install
  @mermaid-js/mermaid-cli@11.4.2`, then set `CARTA_ZARR_MMDC` to its `mmdc` and
  `PUPPETEER_EXECUTABLE_PATH` to a Chrome to do the same; without them diagrams
  render in the browser from the pinned CDN module and are not checked.
- Run the strict Sphinx HTML build to catch documentation warnings and broken
  document references. External links are audited separately by the
  workflow's weekly `linkcheck` job (also runnable from the Actions tab) or
  locally with `-b linkcheck`; it never blocks a pull request, and transient
  external failures need manual interpretation.

## Publishing and versions

The documentation workflow validates pull requests and uploads an HTML artifact.
Pushes to `dev` and manual runs on `dev` can deploy it to GitHub Pages. In the
repository's **Settings → Pages**, select **GitHub Actions** as the publishing
source. Deployment uses the `github-pages` environment and scoped permissions;
pull requests do not receive deployment permissions.

The expected project URL after enabling and successfully deploying Pages is
`https://cartavis.github.io/carta-zarr/`. Site settings have to be enabled by a
repository maintainer; adding the workflow alone does not establish a live site.

Until the first release, the site explicitly describes development code. At the
first release, introduce versioned publication of immutable tagged documentation
alongside the development guide. The current workflow does not yet publish
multiple versions.
