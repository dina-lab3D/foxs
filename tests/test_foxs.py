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
    def run_6lyz(self, *options):
        temporary = tempfile.TemporaryDirectory(prefix="foxs-test-")
        temp = Path(temporary.name)
        for name in ("6lyz.pdb", "lyzexp.dat"):
            shutil.copyfile(DATA / name, temp / name)
        result = subprocess.run(
            [str(FOXS), *options, "6lyz.pdb", "lyzexp.dat"],
            cwd=temp,
            text=True,
            capture_output=True,
            check=False,
        )
        return temporary, temp, result

    def test_simple(self):
        """Fit 6lyz against its experimental SAXS profile."""
        temporary, temp, result = self.run_6lyz("-g")
        with temporary:
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

    def test_gpu_matches_simple_fit_when_available(self):
        """Exercise CUDA numerical parity on CUDA-capable test hosts."""
        temporary, _, result = self.run_6lyz("--gpu")
        with temporary:
            if result.returncode == 2 and "Cannot use --gpu" in result.stderr:
                self.skipTest(result.stderr.strip())
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            match = re.search(r"6lyz\.pdb.*Chi\^2\s+=\s+([\d.]+)", result.stdout)
            self.assertIsNotNone(match, msg="Chi output not found in " + result.stdout)
            self.assertAlmostEqual(float(match.group(1)), 0.20, delta=0.01)


if __name__ == "__main__":
    unittest.main()
