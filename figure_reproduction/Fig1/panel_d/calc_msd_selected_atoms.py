#!/usr/bin/env python3
"""Compute MSD for selected atoms from a LAMMPS dump trajectory.

Configured for the 1500_0 trajectory in this directory:

* selected atom ids: 1..1440
* mass table: type 1 = Fe, type 2 = B
* coordinates: unwrapped xu/yu/zu
* drift correction: mass-weighted COM of the full system

Only the Python standard library is used.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Sequence, Set, Tuple


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_DUMP = SCRIPT_DIR.parents[2] / "two_phase_md" / "SUPERCELL.dump"
DEFAULT_OUT = SCRIPT_DIR / "regenerated_data" / "1500_0_msd.csv"

SELECTED_IDS = set(range(1, 1441))
TYPE_MASSES = {
    1: 55.845,  # Fe
    2: 10.81,  # B
}
REQUIRED_COLUMNS = ("id", "type", "xu", "yu", "zu")

Vector3 = Tuple[float, float, float]


@dataclass
class FrameData:
    index: int
    timestep: int
    natoms: int
    selected_positions: Dict[int, Vector3]
    selected_types: Dict[int, int]
    com: Vector3


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Stream-read a LAMMPS dump and compute first-frame-reference MSD "
            "for the hard-coded selected atom IDs."
        )
    )
    parser.add_argument("--dump", type=Path, default=DEFAULT_DUMP, help="LAMMPS dump file")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="CSV output path")
    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="Read at most this many frames, useful for quick smoke tests",
    )
    return parser.parse_args()


def parse_int_token(token: str, path: Path, line_number: int, column: str) -> int:
    try:
        return int(token)
    except ValueError:
        try:
            value = float(token)
        except ValueError as exc:
            raise ValueError(f"{path}:{line_number}: invalid integer in column {column!r}: {token!r}") from exc
        if value.is_integer():
            return int(value)
        raise ValueError(f"{path}:{line_number}: non-integer value in column {column!r}: {token!r}")


def parse_float_token(token: str, path: Path, line_number: int, column: str) -> float:
    try:
        return float(token)
    except ValueError as exc:
        raise ValueError(f"{path}:{line_number}: invalid float in column {column!r}: {token!r}") from exc


def require_columns(columns: Sequence[str], path: Path, line_number: int) -> Dict[str, int]:
    missing = [column for column in REQUIRED_COLUMNS if column not in columns]
    if missing:
        raise ValueError(
            f"{path}:{line_number}: missing required dump column(s): {', '.join(missing)}"
        )
    return {column: columns.index(column) for column in REQUIRED_COLUMNS}


def iter_lammps_frames(path: Path, selected_ids: Set[int]) -> Iterator[FrameData]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        line_number = 0
        frame_index = 0

        def read_required(context: str) -> str:
            nonlocal line_number
            line = handle.readline()
            if line == "":
                raise ValueError(f"{path}:{line_number + 1}: unexpected EOF while reading {context}")
            line_number += 1
            return line.rstrip("\r\n")

        while True:
            line = handle.readline()
            if line == "":
                return
            line_number += 1
            line = line.rstrip("\r\n")

            if not line.strip():
                continue
            if line != "ITEM: TIMESTEP":
                raise ValueError(f"{path}:{line_number}: expected 'ITEM: TIMESTEP', got {line!r}")

            timestep_line = read_required("timestep")
            try:
                timestep = int(timestep_line.strip())
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid timestep: {timestep_line!r}") from exc

            number_header = read_required("number-of-atoms header")
            if number_header != "ITEM: NUMBER OF ATOMS":
                raise ValueError(
                    f"{path}:{line_number}: expected 'ITEM: NUMBER OF ATOMS', got {number_header!r}"
                )

            natoms_line = read_required("number of atoms")
            try:
                natoms = int(natoms_line.strip())
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: invalid atom count: {natoms_line!r}") from exc

            box_header = read_required("box-bounds header")
            if not box_header.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(
                    f"{path}:{line_number}: expected 'ITEM: BOX BOUNDS', got {box_header!r}"
                )

            for axis in ("x", "y", "z"):
                read_required(f"{axis} box bounds")

            atoms_header = read_required("atoms header")
            if not atoms_header.startswith("ITEM: ATOMS "):
                raise ValueError(
                    f"{path}:{line_number}: expected 'ITEM: ATOMS ...', got {atoms_header!r}"
                )
            columns = atoms_header.split()[2:]
            column_index = require_columns(columns, path, line_number)

            id_idx = column_index["id"]
            type_idx = column_index["type"]
            x_idx = column_index["xu"]
            y_idx = column_index["yu"]
            z_idx = column_index["zu"]

            seen_ids: Set[int] = set()
            selected_positions: Dict[int, Vector3] = {}
            selected_types: Dict[int, int] = {}
            total_mass = 0.0
            weighted_x = 0.0
            weighted_y = 0.0
            weighted_z = 0.0

            for _ in range(natoms):
                atom_line = read_required("atom row")
                parts = atom_line.split()
                if len(parts) < len(columns):
                    raise ValueError(
                        f"{path}:{line_number}: atom row has {len(parts)} fields, "
                        f"expected at least {len(columns)}"
                    )

                atom_id = parse_int_token(parts[id_idx], path, line_number, "id")
                atom_type = parse_int_token(parts[type_idx], path, line_number, "type")
                x = parse_float_token(parts[x_idx], path, line_number, "xu")
                y = parse_float_token(parts[y_idx], path, line_number, "yu")
                z = parse_float_token(parts[z_idx], path, line_number, "zu")

                if atom_id in seen_ids:
                    raise ValueError(f"{path}:{line_number}: duplicated atom id {atom_id}")
                seen_ids.add(atom_id)

                try:
                    mass = TYPE_MASSES[atom_type]
                except KeyError as exc:
                    raise ValueError(
                        f"{path}:{line_number}: unknown atom type {atom_type}; "
                        "update TYPE_MASSES before computing COM drift"
                    ) from exc

                total_mass += mass
                weighted_x += mass * x
                weighted_y += mass * y
                weighted_z += mass * z

                if atom_id in selected_ids:
                    selected_positions[atom_id] = (x, y, z)
                    selected_types[atom_id] = atom_type

            if len(seen_ids) != natoms:
                raise ValueError(
                    f"{path}:{line_number}: frame {frame_index} has {len(seen_ids)} "
                    f"unique atom ids but natoms={natoms}"
                )
            if total_mass <= 0.0:
                raise ValueError(f"{path}:{line_number}: non-positive total mass")

            missing = selected_ids.difference(selected_positions)
            if missing:
                preview = ", ".join(str(atom_id) for atom_id in sorted(missing)[:10])
                if len(missing) > 10:
                    preview += ", ..."
                raise ValueError(
                    f"{path}:{line_number}: frame {frame_index} is missing selected atom id(s): {preview}"
                )

            yield FrameData(
                index=frame_index,
                timestep=timestep,
                natoms=natoms,
                selected_positions=selected_positions,
                selected_types=selected_types,
                com=(weighted_x / total_mass, weighted_y / total_mass, weighted_z / total_mass),
            )
            frame_index += 1


def compute_msd_rows(
    dump_path: Path,
    selected_ids: Set[int],
    max_frames: int | None,
) -> Iterator[List[object]]:
    if max_frames is not None and max_frames <= 0:
        raise ValueError("--max-frames must be a positive integer")

    selected_order = sorted(selected_ids)
    reference_positions: Dict[int, Vector3] | None = None
    reference_com: Vector3 | None = None

    for frame in iter_lammps_frames(dump_path, selected_ids):
        if max_frames is not None and frame.index >= max_frames:
            break

        if reference_positions is None:
            reference_positions = frame.selected_positions
            reference_com = frame.com

        assert reference_positions is not None
        assert reference_com is not None

        drift_x = frame.com[0] - reference_com[0]
        drift_y = frame.com[1] - reference_com[1]
        drift_z = frame.com[2] - reference_com[2]

        sum_x = 0.0
        sum_y = 0.0
        sum_z = 0.0
        for atom_id in selected_order:
            x, y, z = frame.selected_positions[atom_id]
            x0, y0, z0 = reference_positions[atom_id]

            dx = x - x0 - drift_x
            dy = y - y0 - drift_y
            dz = z - z0 - drift_z
            sum_x += dx * dx
            sum_y += dy * dy
            sum_z += dz * dz

        n_selected = len(selected_order)
        msd_x = sum_x / n_selected
        msd_y = sum_y / n_selected
        msd_z = sum_z / n_selected
        msd = msd_x + msd_y + msd_z

        yield [
            frame.index,
            frame.timestep,
            n_selected,
            f"{msd:.16g}",
            f"{msd_x:.16g}",
            f"{msd_y:.16g}",
            f"{msd_z:.16g}",
            f"{drift_x:.16g}",
            f"{drift_y:.16g}",
            f"{drift_z:.16g}",
        ]


def write_msd_csv(dump_path: Path, out_path: Path, max_frames: int | None) -> int:
    if not dump_path.exists():
        raise FileNotFoundError(f"dump file does not exist: {dump_path}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "frame",
        "timestep",
        "n_selected",
        "msd",
        "msd_x",
        "msd_y",
        "msd_z",
        "drift_x",
        "drift_y",
        "drift_z",
    ]

    rows_written = 0
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in compute_msd_rows(dump_path, SELECTED_IDS, max_frames):
            writer.writerow(row)
            rows_written += 1
            if rows_written % 500 == 0:
                print(f"processed {rows_written} frame(s)", file=sys.stderr)

    return rows_written


def main() -> int:
    args = parse_args()
    try:
        rows_written = write_msd_csv(args.dump, args.out, args.max_frames)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {rows_written} MSD row(s) to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
