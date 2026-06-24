"""LAMMPS dump streaming (vendored, self-contained).

Each dump column layout: `id type xu yu zu vx vy vz` (positions unwrapped,
velocities Å/ps). The reader parses the `ITEM: ATOMS ...` header per frame and
does not assume column order.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator
import numpy as np


@dataclass
class Frame:
    """One LAMMPS dump snapshot."""
    timestep: int
    natoms: int
    box: np.ndarray          # (3, 2)
    L: np.ndarray            # (3,)
    atom_id: np.ndarray      # (N,) int32
    atom_type: np.ndarray    # (N,) int8
    pos: np.ndarray          # (N, 3) float32
    vel: np.ndarray | None   # (N, 3) float32 or None


def iter_frames(path: str | Path, stride: int = 1, max_frames: int | None = None,
                want_vel: bool = True) -> Iterator[Frame]:
    """Stream frames from a LAMMPS dump, one at a time (flat memory)."""
    path = Path(path)
    with open(path) as f:
        frame_idx = 0
        yielded = 0
        while True:
            line = f.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            timestep = int(f.readline())
            f.readline()                          # "ITEM: NUMBER OF ATOMS"
            natoms = int(f.readline())
            f.readline()                          # "ITEM: BOX BOUNDS ..."
            box = np.array([[float(x) for x in f.readline().split()] for _ in range(3)])
            L = box[:, 1] - box[:, 0]
            header_tokens = f.readline().split()[2:]   # names after "ITEM: ATOMS"
            ncols = len(header_tokens)
            if frame_idx % stride != 0:
                for _ in range(natoms):
                    f.readline()
                frame_idx += 1
                continue
            block = "".join(f.readline() for _ in range(natoms))
            data = np.fromstring(block, sep=" ").reshape(natoms, ncols)
            ci = {name: header_tokens.index(name) for name in header_tokens}
            yield Frame(
                timestep=timestep, natoms=natoms, box=box, L=L,
                atom_id=data[:, ci["id"]].astype(np.int32),
                atom_type=data[:, ci["type"]].astype(np.int8),
                pos=data[:, [ci["xu"], ci["yu"], ci["zu"]]].astype(np.float32),
                vel=(data[:, [ci["vx"], ci["vy"], ci["vz"]]].astype(np.float32)
                     if want_vel and "vx" in ci else None),
            )
            yielded += 1
            frame_idx += 1
            if max_frames is not None and yielded >= max_frames:
                return


def count_frames(path: str | Path) -> int:
    """Count frames without parsing them."""
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith("ITEM: TIMESTEP"):
                n += 1
    return n
