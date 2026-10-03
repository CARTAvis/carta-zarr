#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "numpy==2.3.1",
#   "zarr==3.2.1",
# ]
#
# [tool.uv]
# # As generate.py: never the system interpreter, and the packages it imports at the versions it pins.
# python-preference = "only-managed"
# ///

"""What tools/zarr-bench/generate.py may do to the place it writes a dataset, asked of Output directly.

--force deletes whatever is at the output, and the dataset being read must survive that whatever the
output is spelt as. That used to hold because main() checked for an overlap before it deleted: two
statements whose order was the guarantee. Output.claim is now the only way to delete or create the
output, and it refuses an overlap first; this asks it, on directories made here, what it does with
each kind of output, and what is left when writing fails part-way."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools" / "zarr-bench"))

import generate  # noqa: E402

DIGEST = "0123456789abcdef"


def dataset(path: Path, digest: str = DIGEST, complete: bool = True) -> Path:
    """A directory that looks like a dataset generate.py wrote, with a manifest saying which."""
    path.mkdir(parents=True)
    (path / "zarr.json").write_text("{}")
    (path / generate.MANIFEST_NAME).write_text(json.dumps({"identity_hash": digest, "complete": complete}))
    return path


class ClaimTest(unittest.TestCase):
    def setUp(self) -> None:
        self._directory = tempfile.TemporaryDirectory()
        self.root = Path(self._directory.name).resolve()
        self.source = dataset(self.root / "source.zarr", digest="the source's own")

    def tearDown(self) -> None:
        self._directory.cleanup()

    def claim(self, path: Path, digest: str = DIGEST, force: bool = True) -> generate.Output:
        return generate.Output.claim(path, [("--source", self.source)], digest, force)

    def test_an_output_over_what_is_read_is_refused_before_anything_is_touched(self) -> None:
        alias = self.root / "alias.zarr"
        alias.symlink_to(self.source)
        for output in (self.source, alias, self.root, self.source / "inside"):
            with self.subTest(output=str(output)):
                with self.assertRaises(SystemExit) as refused:
                    self.claim(output)
                self.assertIn("--source", str(refused.exception))
                self.assertTrue((self.source / "zarr.json").is_file(), f"{output} took the source with it")
                self.assertFalse((self.source / "inside").exists(), "an output inside the source was made")

    def test_the_same_dataset_is_reused_rather_than_written_again(self) -> None:
        output = dataset(self.root / "out.zarr")
        claimed = self.claim(output, force=False)
        self.assertTrue(claimed.reused)
        self.assertEqual(claimed.path, output)
        self.assertTrue((output / "zarr.json").is_file(), "a dataset that was already this one was replaced")

    def test_an_unfinished_dataset_is_not_reused(self) -> None:
        output = dataset(self.root / "out.zarr", complete=False)
        claimed = self.claim(output)
        self.assertFalse(claimed.reused)
        self.assertEqual(list(output.iterdir()), [], "an unfinished dataset was left to be written over")

    def test_another_dataset_is_replaced_only_when_forced(self) -> None:
        output = dataset(self.root / "out.zarr", digest="another")
        with self.assertRaises(SystemExit) as refused:
            self.claim(output, force=False)
        self.assertIn("--force", str(refused.exception))
        self.assertTrue((output / "zarr.json").is_file(), "another dataset was deleted without --force")

        claimed = self.claim(output, force=True)
        self.assertFalse(claimed.reused)
        self.assertTrue(output.is_dir() and list(output.iterdir()) == [], "--force did not leave an empty output")

    def test_a_new_output_is_made_with_its_parents(self) -> None:
        output = self.root / "a" / "b" / "out.zarr"
        claimed = self.claim(output, force=False)
        self.assertFalse(claimed.reused)
        self.assertTrue(output.is_dir())

    def test_writing_that_fails_leaves_nothing_behind(self) -> None:
        for failure in (RuntimeError("a write failed"), SystemExit("refused part-way"), KeyboardInterrupt()):
            with self.subTest(failure=type(failure).__name__):
                claimed = self.claim(self.root / "out.zarr")
                with self.assertRaises(type(failure)):
                    with claimed.writing() as path:
                        (path / "zarr.json").write_text("{}")
                        raise failure
                self.assertFalse((self.root / "out.zarr").exists(), "an unfinished dataset was left behind")
                self.assertTrue((self.source / "zarr.json").is_file())

    def test_writing_that_succeeds_keeps_what_it_wrote(self) -> None:
        claimed = self.claim(self.root / "out.zarr")
        with claimed.writing() as path:
            (path / "zarr.json").write_text("{}")
        self.assertTrue((self.root / "out.zarr" / "zarr.json").is_file())

    def test_a_synthetic_cube_reads_its_template_and_a_rewrite_its_source(self) -> None:
        rewrite = generate.parse_arguments(["--source", str(self.source), "--output", "out", "--chunk", "l=1"])
        self.assertEqual(generate.inputs_read(rewrite), [("--source", self.source)])
        template = dataset(self.root / "template.zarr", digest="a template")
        synthetic = generate.parse_arguments(["--synthetic", "--shape", "l=4,m=4", "--template", str(template),
                                              "--output", "out", "--chunk", "l=1"])
        self.assertEqual(generate.inputs_read(synthetic), [("--template", template)])


if __name__ == "__main__":
    unittest.main()
