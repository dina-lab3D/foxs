#!/usr/bin/env python3
"""Standalone adaptation of IMP FoXS's basic 6lyz regression test."""

from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data"
FOXS = Path(sys.argv.pop(1) if len(sys.argv) > 1 else ROOT / "foxs").resolve()


class FoxsApplicationTest(unittest.TestCase):
    def test_simple(self):
        """Fit 6lyz against its experimental SAXS profile."""
        with tempfile.TemporaryDirectory(prefix="foxs-test-") as temp_name:
            temp = Path(temp_name)
            for name in ("6lyz.pdb", "lyzexp.dat"):
                shutil.copyfile(DATA / name, temp / name)

            result = subprocess.run(
                [str(FOXS), "-g", "6lyz.pdb", "lyzexp.dat"],
                cwd=temp,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)

            match = re.search(r"6lyz\.pdb.*Chi\^2\s+=\s+([\d.]+)", result.stdout)
            self.assertIsNotNone(match, msg="Chi output not found in " + result.stdout)
            self.assertAlmostEqual(float(match.group(1)), 0.20, delta=0.01)

            for name in (
                "6lyz.pdb.dat",
                "6lyz_lyzexp.dat",
                "6lyz_lyzexp.plt",
                "6lyz.plt",
                "6lyz_lyzexp.fit",
            ):
                self.assertTrue((temp / name).is_file(), msg=f"Missing output: {name}")


if __name__ == "__main__":
    unittest.main()
