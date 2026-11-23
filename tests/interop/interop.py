#!/usr/bin/env python3
"""NumPy <-> cnpy interoperability test.

1. NumPy writes .npy/.npz files (all dtypes, C/Fortran order, big-endian, 0-d, empty,
   format version 2.0, savez and savez_compressed); cnpy_interop reads and checks them.
2. cnpy_interop writes files; NumPy loads them, checks the values and requires every .npy
   file to be byte-for-byte identical to what numpy.save produces for the same array.

Usage: interop.py --exe path/to/cnpy_interop --workdir scratch/dir
"""

from __future__ import annotations

import argparse
import io
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import numpy as np

TYPES = [
    "bool", "int8", "uint8", "int16", "uint16", "int32", "uint32", "int64", "uint64",
    "float32", "float64", "complex64", "complex128",
]


def values(dtype: str, n: int) -> np.ndarray:
    """Same formula as test::value<T> in tests/helpers.hpp."""
    i = np.arange(n, dtype=np.int64)
    dt = np.dtype(dtype)
    if dt.kind == "b":
        return i % 3 == 0
    if dt.kind in "iu":
        v = (i * 7 + 3) % 100
        return (v - 50 if dt.kind == "i" else v).astype(dt)
    if dt.kind == "f":
        return (i * 0.25 - 8).astype(dt)
    return ((i * 0.25 - 8) + 1j * (-i * 0.5)).astype(dt)


def cube(dtype: str) -> np.ndarray:
    return values(dtype, 60).reshape(3, 4, 5)


def generate(d: Path) -> None:
    for t in TYPES:
        a = cube(t)
        np.save(d / f"{t}_c.npy", a)
        np.save(d / f"{t}_f.npy", np.asfortranarray(a))
        np.save(d / f"{t}_be.npy", a.astype(a.dtype.newbyteorder(">")))
        np.save(d / f"{t}_1d.npy", values(t, 17))
        np.save(d / f"{t}_0d.npy", values(t, 1).reshape(()))
        np.save(d / f"{t}_empty.npy", values(t, 0).reshape(0, 3))
        with open(d / f"{t}_v2.npy", "wb") as fh:
            np.lib.format.write_array(fh, values(t, 6).reshape(2, 3), version=(2, 0))
        with open(d / f"{t}_v3.npy", "wb") as fh:
            np.lib.format.write_array(fh, values(t, 6).reshape(2, 3), version=(3, 0))
    np.savez(d / "savez.npz", **{t: cube(t) for t in TYPES})
    np.savez_compressed(d / "savez_compressed.npz", **{t: cube(t) for t in TYPES})
    np.save(d / "strings.npy", np.array(["hello", "wörld"], dtype="U5"))
    np.save(d / "structured.npy", np.zeros(2, dtype=[("a", "<i4"), ("b", "<f8")]))


def verify(d: Path) -> list[str]:
    errors: list[str] = []

    def same(name: str, got: np.ndarray, want: np.ndarray) -> None:
        if got.dtype != want.dtype or got.shape != want.shape or not np.array_equal(got, want):
            errors.append(f"{name}: got {got.dtype}{got.shape}, want {want.dtype}{want.shape} (or values differ)")

    for t in TYPES:
        a = cube(t)
        expected = {
            f"{t}_c.npy": a,
            f"{t}_f.npy": np.asfortranarray(a),
            f"{t}_1d.npy": values(t, 17),
            f"{t}_0d.npy": values(t, 1).reshape(()),
            f"{t}_empty.npy": values(t, 0).reshape(0, 3),
            f"{t}_appended.npy": a,
        }
        for name, want in expected.items():
            path = d / name
            got = np.load(path)
            same(name, got, want)
            if name.endswith("_f.npy") and not got.flags.f_contiguous:
                errors.append(f"{name}: not loaded as Fortran-ordered")
            ref = io.BytesIO()
            np.save(ref, want)
            if path.read_bytes() != ref.getvalue():
                errors.append(f"{name}: not byte-identical to numpy.save output")

    for archive, method in (("cnpy_savez.npz", zipfile.ZIP_STORED),
                            ("cnpy_savez_compressed.npz", zipfile.ZIP_DEFLATED)):
        path = d / archive
        with zipfile.ZipFile(path) as zf:
            bad = zf.testzip()  # verifies every CRC-32
            if bad is not None:
                errors.append(f"{archive}: CRC error in {bad}")
            if any(info.compress_type != method for info in zf.infolist()):
                errors.append(f"{archive}: unexpected compression method")
        with np.load(path) as z:
            if sorted(z.files) != sorted(TYPES):
                errors.append(f"{archive}: members {sorted(z.files)}")
            for t in TYPES:
                same(f"{archive}/{t}", z[t], cube(t))

    with np.load(d / "cnpy_appended.npz") as z:
        if z.files != ["a", "b", "c"]:
            errors.append(f"cnpy_appended.npz: members {z.files}")
        same("cnpy_appended.npz/a", z["a"], values("float64", 10))
        same("cnpy_appended.npz/b", z["b"], values("int32", 1000))
        same("cnpy_appended.npz/c", z["c"], values("float32", 6).reshape(2, 3))

    with zipfile.ZipFile(d / "cnpy_many.npz") as zf:
        if zf.testzip() is not None or len(zf.infolist()) != 65536 + 5:
            errors.append("cnpy_many.npz: bad zip64 archive")
    with np.load(d / "cnpy_many.npz") as z:
        if len(z.files) != 65536 + 5 or int(z["m65540"]) != 65540 or z["m7"].shape != ():
            errors.append("cnpy_many.npz: wrong content")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--exe", required=True, type=Path)
    parser.add_argument("--workdir", required=True, type=Path)
    args = parser.parse_args()

    shutil.rmtree(args.workdir, ignore_errors=True)
    from_numpy = args.workdir / "from_numpy"
    from_cnpy = args.workdir / "from_cnpy"
    from_numpy.mkdir(parents=True)
    from_cnpy.mkdir(parents=True)

    print(f"numpy {np.__version__}")
    generate(from_numpy)
    ok = subprocess.run([str(args.exe), "read", str(from_numpy)]).returncode == 0
    ok = subprocess.run([str(args.exe), "write", str(from_cnpy)]).returncode == 0 and ok

    try:
        errors = verify(from_cnpy)
    except Exception as exc:  # e.g. a file cnpy failed to write
        errors = [f"verify crashed: {exc!r}"]
    for e in errors:
        print("FAIL:", e, file=sys.stderr)
    if not errors:
        print("numpy verify: all checks passed (values and byte-identical .npy files)")
    return 0 if ok and not errors else 1


if __name__ == "__main__":
    sys.exit(main())
