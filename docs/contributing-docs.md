# Maintain the documentation

## Build locally

Install Doxygen (`brew install doxygen` on macOS or `sudo apt-get install doxygen`
on Ubuntu), then run from the repository root with Python 3.11 or newer:

```sh
python3 -m venv build-docs/venv
build-docs/venv/bin/python -m pip install -r docs/requirements.txt
build-docs/venv/bin/python -m sphinx -b html -W --keep-going docs build-docs/html
python3 -m http.server 8000 --directory build-docs/html
```

Open `http://localhost:8000`. The docs build runs Doxygen automatically and does
not configure or compile TensorStore. Build output and virtual environments
stay under the ignored `build-docs/` directory.

## Sources and checks

- Edit Markdown under `docs/` for tutorials and guides.
- Edit public header comments for API contracts. A small Doxygen input filter
  recognises the existing `//` prose without changing the installed headers.
- Glossary and changelog pages read the root files rather than maintaining
  copies. Storage tuning retains its existing page. During the site build,
  repository-only links in the changelog and storage guide point to GitHub.
- The quickstart includes `examples/first_read.cc` directly. The package-consumer
  test builds and runs it against a committed fixture; keep example changes in
  step with the public API.
- Keep diagrams as Mermaid text next to their explanation. Rendering uses the
  pinned Mermaid browser module and requires access to its CDN.
- Run the strict Sphinx HTML build to catch documentation warnings and broken
  document references. To audit external links separately, run the same command
  with `-b linkcheck`; transient external failures need manual interpretation.

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
