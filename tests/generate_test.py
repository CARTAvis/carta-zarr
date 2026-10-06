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

# This file is part of the CARTA Image Viewer: https://github.com/CARTAvis/carta-zarr
# Copyright 2026- Academia Sinica Institute of Astronomy and Astrophysics (ASIAA),
# Associated Universities, Inc. (AUI) and the Inter-University Institute for Data Intensive Astronomy (IDIA)
# SPDX-License-Identifier: GPL-3.0-or-later

"""What tools/zarr-bench/generate.py may do to the place it writes a dataset, asked of Output directly.

--force deletes whatever is at the output, and the dataset being read must survive that whatever the
output is spelt as. That used to hold because main() checked for an overlap before it deleted: two
statements whose order was the guarantee. Output.claim is now the only way to delete or create the
output, and it refuses an overlap first; this asks it, on directories made here, what it does with
each kind of output, and what is left when writing fails part-way."""

from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any

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

    def test_an_output_holding_what_the_source_links_to_is_refused(self) -> None:
        # A dataset's arrays are as often as not links to where the bytes are. --force deletes the
        # output, and an output that holds what a link of the source points at takes those bytes
        # with it however far from the source's own directory they are.
        elsewhere = self.root / "elsewhere"
        (elsewhere / "SKY").mkdir(parents=True)
        (elsewhere / "SKY" / "zarr.json").write_text("{}")
        (self.source / "SKY").symlink_to(elsewhere / "SKY")
        for output in (elsewhere, elsewhere / "SKY", elsewhere / "SKY" / "inside"):
            with self.subTest(output=str(output)):
                with self.assertRaises(SystemExit) as refused:
                    self.claim(output)
                self.assertIn("--source", str(refused.exception))
                self.assertTrue((elsewhere / "SKY" / "zarr.json").is_file(), f"{output} took the linked pixels")
                self.assertFalse((elsewhere / "SKY" / "inside").exists(), "an output inside the link was made")
        self.claim(self.root / "beside")

    def test_a_link_back_up_the_source_does_not_loop(self) -> None:
        (self.source / "loop").symlink_to(self.source)
        self.claim(self.root / "beside")

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

    def test_the_same_dataset_is_written_again_when_forced(self) -> None:
        output = dataset(self.root / "out.zarr")
        claimed = self.claim(output, force=True)
        self.assertFalse(claimed.reused, "--force reused what it was asked to replace")
        self.assertEqual(list(output.iterdir()), [], "--force did not leave an empty output")

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



class FindFlagTest(unittest.TestCase):
    """The flag a rewrite rewrites beside its image is the one carta-zarr masks the image with, by the
    rule of src/schema/xradio/flag.cc: the image's own `flag` attribute, else the one the root's
    data_groups give the image as its sky, else the only flag-typed boolean array of its shape. A
    boolean is a Zarr bool or xarray's int8 marked dtype "bool", which is what XRADIO writes. A flag
    left behind keeps the source's layout, and the layouts compared would each be a mix."""

    def setUp(self) -> None:
        self._directory = tempfile.TemporaryDirectory()
        self.root = Path(self._directory.name)

    def tearDown(self) -> None:
        self._directory.cleanup()

    def array(self, name: str, data_type: str = "float32", **attributes: Any) -> dict[str, Any]:
        metadata = {"zarr_format": 3, "node_type": "array", "shape": [1, 2, 3, 4, 5], "data_type": data_type,
                    "dimension_names": ["time", "frequency", "polarization", "l", "m"], "attributes": attributes}
        (self.root / name).mkdir()
        (self.root / name / "zarr.json").write_text(json.dumps(metadata))
        return metadata

    def find(self, image: dict[str, Any], groups: dict[str, Any] | None = None) -> str | None:
        root = {"attributes": {"data_groups": groups}} if groups is not None else {"attributes": {}}
        return generate.find_flag(self.root, "SKY", image, generate.member_names(self.root), root)

    def test_xarrays_encoding_of_a_bool_is_a_flag(self) -> None:
        image = self.array("SKY")
        self.array("FLAG_SKY", "int8", type="flag", dtype="bool")
        self.assertEqual(self.find(image), "FLAG_SKY")

    def test_an_int8_that_is_not_marked_bool_is_not(self) -> None:
        image = self.array("SKY")
        self.array("FLAG_SKY", "int8", type="flag")
        self.assertIsNone(self.find(image))

    def test_a_data_group_declares_its_skys_flag(self) -> None:
        image = self.array("SKY")
        self.array("FLAG_SKY", "int8", type="flag", dtype="bool")
        self.array("OTHER_FLAG", "bool", type="flag")
        self.assertEqual(self.find(image, {"base": {"sky": "SKY", "flag": "FLAG_SKY"}}), "FLAG_SKY")

    def test_the_images_own_attribute_outranks_a_data_group(self) -> None:
        image = self.array("SKY", flag="MINE")
        self.array("MINE", "bool", type="flag")
        self.array("FLAG_SKY", "bool", type="flag")
        self.assertEqual(self.find(image, {"base": {"sky": "SKY", "flag": "FLAG_SKY"}}), "MINE")

    def test_two_groups_naming_different_flags_is_refused(self) -> None:
        image = self.array("SKY")
        self.array("A", "bool", type="flag")
        self.array("B", "bool", type="flag")
        with self.assertRaises(SystemExit):
            self.find(image, {"one": {"sky": "SKY", "flag": "A"}, "two": {"sky": "SKY", "flag": "B"}})

    def test_another_images_flag_is_not_guessed_to_be_this_ones(self) -> None:
        image = self.array("SKY")
        self.array("FLAG_MODEL", "bool", type="flag")
        self.assertIsNone(self.find(image, {"model": {"sky": "MODEL", "flag": "FLAG_MODEL"}}))

    def test_the_committed_conversion_names_its_flag(self) -> None:
        fixture = Path(__file__).resolve().parent / "data" / "images" / "zarr" / "xradio" / "conformance_flagged"
        root = generate.read_metadata(fixture)
        image = generate.read_metadata(fixture / "SKY")
        self.assertEqual(generate.find_flag(fixture, "SKY", image, generate.member_names(fixture), root), "FLAG_SKY")


class IdentityTest(unittest.TestCase):
    """A rewrite is known by what it was made from, and that is the source's pixels, flags and
    coordinates as much as its metadata: a source that changes under the same path is another
    dataset, and reusing the old rewrite would measure pixels that are no longer there."""

    def setUp(self) -> None:
        self._directory = tempfile.TemporaryDirectory()
        self.source = Path(self._directory.name).resolve() / "source.zarr"
        (self.source / "SKY" / "c" / "0").mkdir(parents=True)
        (self.source / "l" / "c").mkdir(parents=True)
        (self.source / "zarr.json").write_text("{}")
        (self.source / "SKY" / "zarr.json").write_text(json.dumps({"shape": [4]}))
        (self.source / "SKY" / "c" / "0" / "0").write_bytes(b"pixel 11")
        (self.source / "l" / "c" / "0").write_bytes(b"coordinates")

    def tearDown(self) -> None:
        self._directory.cleanup()

    def identity(self) -> dict:
        return generate.identity(generate.parse_arguments(
            ["--source", str(self.source), "--output", "out", "--chunk", "l=1"]))

    def test_a_source_whose_pixels_change_is_another_dataset(self) -> None:
        before = self.identity()
        self.assertEqual(self.identity(), before, "the same source came to two identities")
        for changed in (self.source / "SKY" / "c" / "0" / "0", self.source / "l" / "c" / "0"):
            with self.subTest(changed=str(changed.relative_to(self.source))):
                was = self.identity()
                stat = changed.stat()
                changed.write_bytes(changed.read_bytes().replace(b"1", b"2"))
                # The same size, and a time a coarse filesystem could not tell apart from the old one
                # had it not been moved on: what is left to tell them apart is that it was written.
                os.utime(changed, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
                self.assertNotEqual(self.identity(), was, f"{changed.name} changed and the identity did not")

    def test_a_source_whose_linked_pixels_change_is_another_dataset(self) -> None:
        # The pixels are where the source's SKY links to: that they changed is what the rewrite has
        # to see, however the source reaches them.
        pixels = self.source.parent / "pixels"
        (self.source / "SKY").rename(pixels)
        (self.source / "SKY").symlink_to(pixels)
        changed = pixels / "c" / "0" / "0"
        was = self.identity()
        stat = changed.stat()
        changed.write_bytes(b"pixel 22")
        os.utime(changed, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
        self.assertNotEqual(self.identity(), was, "linked pixels changed and the identity did not")

    def test_a_template_whose_arrays_change_is_another_synthetic_dataset(self) -> None:
        # A synthetic cube takes its coordinates and their metadata from the template's arrays, not
        # only from its root document: a template rechunked along frequency writes another dataset.
        frequency = self.source / "frequency" / "zarr.json"
        frequency.parent.mkdir()
        frequency.write_text(json.dumps({"chunk_grid": {"configuration": {"chunk_shape": [4]}}}))

        def synthetic() -> dict:
            return generate.identity(generate.parse_arguments(
                ["--synthetic", "--shape", "time=1,frequency=4,polarization=1,l=8,m=8", "--template", str(self.source), "--output", "out", "--chunk", "l=1"]))

        was = synthetic()
        self.assertEqual(synthetic(), was, "the same template came to two identities")
        stat = frequency.stat()
        frequency.write_text(json.dumps({"chunk_grid": {"configuration": {"chunk_shape": [2]}}}))
        os.utime(frequency, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
        self.assertNotEqual(synthetic(), was, "a template array changed and the identity did not")

    def test_where_it_is_written_and_how_fast_are_not_its_identity(self) -> None:
        here = self.identity()
        elsewhere = generate.identity(generate.parse_arguments(
            ["--source", str(self.source), "--output", "elsewhere", "--chunk", "l=1", "--workers", "3"]))
        self.assertEqual(here, elsewhere)


class SyntheticCoordinatesTest(unittest.TestCase):
    """A synthetic cube's right ascension and declination are each as large as a plane -- 8 GiB apiece
    at 32768 square -- so they are written in blocks through fill(), as the pixels are, and not built
    whole in memory whatever --block-mib said."""

    def test_the_sky_direction_is_written_in_blocks_and_is_the_projection_inverted(self) -> None:
        import numpy as np
        import zarr

        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / "synthetic.zarr"
            args = generate.parse_arguments(
                ["--synthetic", "--shape", "frequency=1,polarization=1,l=1024,m=1024", "--output", str(out),
                 "--chunk", "l=128,m=128", "--block-mib", "1", "--workers", "1"])
            filled: dict[str, int] = {}
            real_fill = generate.fill

            def fill(job: Any, shape: list[int], unit: tuple[int, ...], itemsize: int, given: Any) -> None:
                filled[Path(job.target).name] = len(generate.plan_blocks(shape, unit, itemsize, given.block_mib << 20))
                real_fill(job, shape, unit, itemsize, given)

            generate.fill = fill
            try:
                generate.synthesize(args, out)
            finally:
                generate.fill = real_fill

            for name in ("right_ascension", "declination"):
                self.assertIn(name, filled, f"{name} was built whole rather than written through fill()")
                # 1024 x 1024 float64 is 8 MiB, so 1 MiB blocks take eight.
                self.assertGreater(filled[name], 1, f"{name} was written in one block")

            root = json.loads((Path(generate.DEFAULT_TEMPLATE) / "zarr.json").read_text())
            ra0, dec0 = root["attributes"]["coordinate_system_info"]["reference_direction"]["data"]
            l_axis = np.asarray(zarr.open_array(str(out / "l"), mode="r")[...])[:, None]
            m_axis = np.asarray(zarr.open_array(str(out / "m"), mode="r")[...])[None, :]
            n = np.sqrt(np.maximum(0.0, 1.0 - l_axis * l_axis - m_axis * m_axis))
            expected = {
                "declination": np.arcsin(m_axis * np.cos(dec0) + n * np.sin(dec0)),
                "right_ascension": ra0 + np.arctan2(l_axis, n * np.cos(dec0) - m_axis * np.sin(dec0)),
            }
            for name, values in expected.items():
                written = np.asarray(zarr.open_array(str(out / name), mode="r")[...])
                self.assertEqual(written.shape, (1024, 1024))
                np.testing.assert_allclose(written, values, rtol=0, atol=1e-15, err_msg=name)

if __name__ == "__main__":
    unittest.main()
