"""Expose the public headers' existing prose to Doxygen without editing the API."""

import re
import sys
from pathlib import Path


for line in Path(sys.argv[1]).read_text().splitlines(keepends=True):
    # Doxygen needs /// before a declaration and ///< after a field. Keep existing
    # Doxygen comments and ordinary trailing comments on braces untouched.
    line = re.sub(r"^(\s*)//(?![/!])", r"\1///", line)
    line = re.sub(r";(\s*)//(?![/!])", r";\1///<", line)
    # Keep punctuation as Unicode: Breathe renders Doxygen's <ndash/> as a
    # literal HTML entity rather than a dash in ordinary paragraph text.
    if line.lstrip().startswith("///"):
        line = re.sub(r"(?<=\s)--(?=\s|$)", "—", line)
    sys.stdout.write(line)
