"""The user guide and public API build independently of the C++ library."""

import os
import re
import shlex
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from sphinx.errors import ExtensionError


DOCS = Path(__file__).resolve().parent
ROOT = DOCS.parent
XML = ROOT / "build-docs" / "doxygen"
REPO_URL = "https://github.com/CARTAvis/carta-zarr"

project = "carta-zarr"
copyright = "2026, ASIAA, AUI and IDIA"
author = "CARTA contributors"
version = re.search(r"project\(carta_zarr VERSION ([\d.]+)", (ROOT / "CMakeLists.txt").read_text())[1]
release = f"{version} (development)"

extensions = ["myst_parser", "breathe", "sphinxcontrib.mermaid"]
source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
# read-api.md and sidecar.md are local drafts that .gitignore keeps out of the
# repository; Sphinx does not read .gitignore, so a checkout holding them must skip them.
exclude_patterns = ["adr/**", "archive/**", "README.md", "read-api.md", "sidecar.md"]
myst_enable_extensions = ["colon_fence"]
myst_heading_anchors = 3
primary_domain = "cpp"
highlight_language = "cpp"
breathe_projects = {"carta-zarr": str(XML / "xml")}
breathe_default_project = "carta-zarr"
html_theme = "furo"
html_title = "carta-zarr documentation"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "source_repository": REPO_URL + "/",
    "source_branch": "dev",
    "source_directory": "docs/",
    "light_css_variables": {"color-brand-primary": "#176477", "color-brand-content": "#176477"},
    "dark_css_variables": {"color-brand-primary": "#6bd0de", "color-brand-content": "#6bd0de"},
}
# GitHub adds heading anchors with JavaScript, so a plain fetch never finds them.
linkcheck_anchors_ignore_for_url = [r"https://github\.com/"]
linkcheck_retries = 2
# The in-browser Mermaid module; mermaid-cli (11.4.2 in CI) is a separate package with its own numbering.
mermaid_version = "11.4.1"
# With CARTA_ZARR_MMDC naming a mermaid-cli executable (and PUPPETEER_EXECUTABLE_PATH a
# Chrome), diagrams render to SVG at build time: syntax errors fail the build and the
# site needs no CDN. Without it they render in the browser from the pinned CDN module.
if os.environ.get("CARTA_ZARR_MMDC"):
    mermaid_cmd = os.environ["CARTA_ZARR_MMDC"]
    mermaid_output_format = "svg"
    mermaid_params = ["--puppeteerConfigFile", str(DOCS / "mermaid-puppeteer.json"), "--backgroundColor", "white"]


def generate_api(app):
    config = (DOCS / "Doxyfile").read_text()
    filter_command = shlex.join([sys.executable, str(DOCS / "header_filter.py")])
    config += f'\nOUTPUT_DIRECTORY = "{XML}"\nINPUT_FILTER = "{filter_command}"\n'
    subprocess.run(["doxygen", "-"], input=config, text=True, cwd=ROOT, check=True)
    # Recent Doxygen versions append primary-template members to the explicit
    # Result<void> specialization, even though it does not inherit them. Keep
    # only its own declarations, avoiding both bogus methods and duplicate IDs.
    index = ET.parse(XML / "xml" / "index.xml")
    for compound in index.findall("compound"):
        if compound.findtext("name") != "carta::zarr::Result< void >":
            continue
        refid = compound.attrib["refid"]
        path = XML / "xml" / f"{refid}.xml"
        tree = ET.parse(path)
        for section in tree.findall("compounddef/sectiondef"):
            for member in list(section):
                if not member.attrib["id"].startswith(refid + "_"):
                    section.remove(member)
        tree.write(path, encoding="utf-8", xml_declaration=True)


def repository_links(app, docname, source):
    # Reuse the root changelog and glossary while resolving repo-relative links
    # for the website. Register the dependency so incremental builds notice
    # edits, and fail loudly when a source no longer has the shape this expects.
    if docname == "changelog":
        path = ROOT / "CHANGELOG.md"
        app.env.note_dependency(str(path))
        text = path.read_text().replace("](docs/adr/", f"]({REPO_URL}/blob/dev/docs/adr/")
        source[0] = text
        leftover = re.findall(r"\]\((?:docs|tools|\.\.)/[^)]*\)", text)
    elif docname == "storage-tuning":
        text = source[0].replace("](../tools/", f"]({REPO_URL}/blob/dev/tools/")
        source[0] = text.replace("](adr/", f"]({REPO_URL}/blob/dev/docs/adr/")
        leftover = re.findall(r"\]\((?:tools|\.\./)[^)]*\)", source[0])
    elif docname == "glossary":
        # glossary.md includes GLOSSARY.md from line 2, skipping its title.
        path = ROOT / "GLOSSARY.md"
        app.env.note_dependency(str(path))
        if not path.read_text().startswith("# "):
            raise ExtensionError("GLOSSARY.md must start with a '# ' title: glossary.md skips its first line")
        return
    else:
        return
    if leftover:
        raise ExtensionError(f"{docname}: repository-relative links not rewritten for the website: {leftover}")


def setup(app):
    app.connect("builder-inited", generate_api)
    app.connect("source-read", repository_links)
