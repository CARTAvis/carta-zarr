"""The user guide and public API build independently of the C++ library."""

import re
import shlex
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


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
mermaid_version = "11.4.1"


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
    # Reuse the root changelog while resolving its repo-relative ADR links for
    # the website. Register the dependency so incremental builds notice edits.
    if docname == "changelog":
        path = ROOT / "CHANGELOG.md"
        app.env.note_dependency(str(path))
        source[0] = path.read_text().replace("](docs/adr/", f"]({REPO_URL}/blob/dev/docs/adr/")
    elif docname == "storage-tuning":
        source[0] = source[0].replace("](../tools/", f"]({REPO_URL}/blob/dev/tools/")
        source[0] = source[0].replace("](adr/", f"]({REPO_URL}/blob/dev/docs/adr/")


def setup(app):
    app.connect("builder-inited", generate_api)
    app.connect("source-read", repository_links)
