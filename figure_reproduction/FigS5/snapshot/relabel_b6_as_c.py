#!/usr/bin/env python3
"""Relabel dynamically six-coordinated B atoms as type 3 in a LAMMPS dump."""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Dict, Iterator, List, Sequence, Tuple


DEFAULT_DUMP = Path("SUPERCELL-0.8-0.dump")
DEFAULT_OUT = Path("SUPERCELL-0.8-0_B6C.dump")
DEFAULT_CUTOFF = 2.37

FE_TYPE = 1
B_TYPE = 2
C_TYPE = 3
REQUIRED_COLUMNS = ("id", "type", "xu", "yu", "zu")
TOKEN_RE = re.compile(r"\S+")

Vector3 = Tuple[float, float, float]
Bounds3 = Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]


@dataclass
class AtomRow:
    raw: str
    atom_id: int
    atom_type: int
    position: Vector3
    type_span: Tuple[int, int]


@dataclass
class DumpFrame:
    index: int
    timestep: int
    natoms: int
    bounds: Bounds3
    header_lines: List[str]
    atom_rows: List[AtomRow]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dump", type=Path, default=DEFAULT_DUMP, help="input LAMMPS dump")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="output relabeled dump")
    parser.add_argument(
        "--summary",
        type=Path,
        default=None,
        help="per-frame CSV summary; defaults to OUT stem plus _counts.csv",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=DEFAULT_CUTOFF,
        help="inclusive B-B bond cutoff in Angstrom",
    )
    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="process at most this many frames, useful for testing",
    )
    parser.add_argument(
        "--verify-first-frame",
        action="store_true",
        help="compare the first-frame cell-list coordination against brute force",
    )
    args = parser.parse_args()
    if args.summary is None:
        args.summary = args.out.with_name(f"{args.out.stem}_counts.csv")
    return args


def parse_int(token: str, path: Path, line_number: int, column: str) -> int:
    try:
        return int(token)
    except ValueError:
        try:
            value = float(token)
        except ValueError as exc:
            raise ValueError(
                f"{path}:{line_number}: invalid integer in {column!r}: {token!r}"
            ) from exc
        if value.is_integer():
            return int(value)
        raise ValueError(f"{path}:{line_number}: non-integer value in {column!r}: {token!r}")


def parse_float(token: str, path: Path, line_number: int, column: str) -> float:
    try:
        return float(token)
    except ValueError as exc:
        raise ValueError(
            f"{path}:{line_number}: invalid float in {column!r}: {token!r}"
        ) from exc


def iter_dump_frames(path: Path) -> Iterator[DumpFrame]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        line_number = 0
        frame_index = 0

        def read_required(context: str) -> str:
            nonlocal line_number
            line = handle.readline()
            if line == "":
                raise ValueError(f"{path}:{line_number + 1}: unexpected EOF while reading {context}")
            line_number += 1
            return line

        while True:
            timestep_header = handle.readline()
            if timestep_header == "":
                return
            line_number += 1
            if timestep_header.rstrip("\r\n") != "ITEM: TIMESTEP":
                raise ValueError(
                    f"{path}:{line_number}: expected 'ITEM: TIMESTEP', "
                    f"got {timestep_header.rstrip()!r}"
                )

            timestep_line = read_required("timestep")
            try:
                timestep = int(timestep_line.strip())
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: invalid timestep {timestep_line.strip()!r}"
                ) from exc

            number_header = read_required("number-of-atoms header")
            if number_header.rstrip("\r\n") != "ITEM: NUMBER OF ATOMS":
                raise ValueError(
                    f"{path}:{line_number}: expected 'ITEM: NUMBER OF ATOMS', "
                    f"got {number_header.rstrip()!r}"
                )
            natoms_line = read_required("number of atoms")
            try:
                natoms = int(natoms_line.strip())
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: invalid atom count {natoms_line.strip()!r}"
                ) from exc
            if natoms <= 0:
                raise ValueError(f"{path}:{line_number}: atom count must be positive")

            box_header = read_required("box-bounds header")
            box_tokens = box_header.split()
            if box_tokens[:3] != ["ITEM:", "BOX", "BOUNDS"]:
                raise ValueError(
                    f"{path}:{line_number}: expected 'ITEM: BOX BOUNDS ...', "
                    f"got {box_header.rstrip()!r}"
                )
            if box_tokens[3:] != ["pp", "pp", "pp"]:
                raise ValueError(
                    f"{path}:{line_number}: only orthogonal fully periodic 'pp pp pp' boxes "
                    f"are supported, got {' '.join(box_tokens[3:])!r}"
                )

            bound_lines: List[str] = []
            parsed_bounds: List[Tuple[float, float]] = []
            for axis in ("x", "y", "z"):
                bound_line = read_required(f"{axis} box bounds")
                bound_lines.append(bound_line)
                parts = bound_line.split()
                if len(parts) != 2:
                    raise ValueError(
                        f"{path}:{line_number}: triclinic or malformed box bounds are unsupported"
                    )
                lower = parse_float(parts[0], path, line_number, f"{axis}lo")
                upper = parse_float(parts[1], path, line_number, f"{axis}hi")
                if upper <= lower:
                    raise ValueError(
                        f"{path}:{line_number}: {axis} upper bound must exceed lower bound"
                    )
                parsed_bounds.append((lower, upper))

            atoms_header = read_required("atoms header")
            if not atoms_header.startswith("ITEM: ATOMS "):
                raise ValueError(
                    f"{path}:{line_number}: expected 'ITEM: ATOMS ...', "
                    f"got {atoms_header.rstrip()!r}"
                )
            columns = atoms_header.split()[2:]
            missing = [column for column in REQUIRED_COLUMNS if column not in columns]
            if missing:
                raise ValueError(
                    f"{path}:{line_number}: missing required column(s): {', '.join(missing)}"
                )
            indices = {column: columns.index(column) for column in REQUIRED_COLUMNS}

            atom_rows: List[AtomRow] = []
            seen_ids = set()
            for _ in range(natoms):
                raw = read_required("atom row")
                matches = list(TOKEN_RE.finditer(raw))
                if len(matches) < len(columns):
                    raise ValueError(
                        f"{path}:{line_number}: atom row has {len(matches)} fields, "
                        f"expected at least {len(columns)}"
                    )
                tokens = [match.group(0) for match in matches]
                atom_id = parse_int(tokens[indices["id"]], path, line_number, "id")
                atom_type = parse_int(tokens[indices["type"]], path, line_number, "type")
                if atom_type not in (FE_TYPE, B_TYPE):
                    raise ValueError(
                        f"{path}:{line_number}: expected only type 1 (Fe) and type 2 (B), "
                        f"got type {atom_type}"
                    )
                if atom_id in seen_ids:
                    raise ValueError(f"{path}:{line_number}: duplicate atom id {atom_id}")
                seen_ids.add(atom_id)
                position = (
                    parse_float(tokens[indices["xu"]], path, line_number, "xu"),
                    parse_float(tokens[indices["yu"]], path, line_number, "yu"),
                    parse_float(tokens[indices["zu"]], path, line_number, "zu"),
                )
                atom_rows.append(
                    AtomRow(
                        raw=raw,
                        atom_id=atom_id,
                        atom_type=atom_type,
                        position=position,
                        type_span=matches[indices["type"]].span(),
                    )
                )

            if len(seen_ids) != natoms:
                raise ValueError(
                    f"{path}:{line_number}: frame {frame_index} has {len(seen_ids)} unique ids, "
                    f"expected {natoms}"
                )

            yield DumpFrame(
                index=frame_index,
                timestep=timestep,
                natoms=natoms,
                bounds=(parsed_bounds[0], parsed_bounds[1], parsed_bounds[2]),
                header_lines=[
                    timestep_header,
                    timestep_line,
                    number_header,
                    natoms_line,
                    box_header,
                    *bound_lines,
                    atoms_header,
                ],
                atom_rows=atom_rows,
            )
            frame_index += 1


def coordination_cell_list(
    frame: DumpFrame,
    cutoff: float,
) -> Tuple[List[AtomRow], List[int]]:
    lengths = tuple(upper - lower for lower, upper in frame.bounds)
    if cutoff >= min(lengths) / 2.0:
        raise ValueError(
            f"frame {frame.index}: cutoff {cutoff:g} must be less than half the shortest box length"
        )

    b_rows = [row for row in frame.atom_rows if row.atom_type == B_TYPE]
    lowers = tuple(bound[0] for bound in frame.bounds)
    cell_counts = tuple(max(1, int(length // cutoff)) for length in lengths)
    cell_widths = tuple(lengths[axis] / cell_counts[axis] for axis in range(3))
    wrapped: List[Vector3] = []
    cells: DefaultDict[Tuple[int, int, int], List[int]] = defaultdict(list)

    for index, row in enumerate(b_rows):
        position = tuple(
            (row.position[axis] - lowers[axis]) % lengths[axis] for axis in range(3)
        )
        wrapped.append((position[0], position[1], position[2]))
        key = tuple(
            min(cell_counts[axis] - 1, int(position[axis] / cell_widths[axis]))
            for axis in range(3)
        )
        cells[(key[0], key[1], key[2])].append(index)

    coordination = [0] * len(b_rows)
    cutoff_squared = cutoff * cutoff
    lx, ly, lz = lengths
    hx, hy, hz = lx / 2.0, ly / 2.0, lz / 2.0
    nx, ny, nz = cell_counts

    def count_pair(first: int, second: int) -> None:
        x1, y1, z1 = wrapped[first]
        x2, y2, z2 = wrapped[second]
        dx = abs(x2 - x1)
        if dx > hx:
            dx = lx - dx
        if dx > cutoff:
            return
        dy = abs(y2 - y1)
        if dy > hy:
            dy = ly - dy
        if dy > cutoff:
            return
        dz = abs(z2 - z1)
        if dz > hz:
            dz = lz - dz
        if dz > cutoff:
            return
        if dx * dx + dy * dy + dz * dz <= cutoff_squared:
            coordination[first] += 1
            coordination[second] += 1

    forward_offsets = [
        (dx, dy, dz)
        for dx in (-1, 0, 1)
        for dy in (-1, 0, 1)
        for dz in (-1, 0, 1)
        if (dx, dy, dz) > (0, 0, 0)
    ]

    for (ix, iy, iz), members in cells.items():
        for member_index, first in enumerate(members):
            for second in members[member_index + 1 :]:
                count_pair(first, second)

        seen_neighbor_keys = set()
        for dx, dy, dz in forward_offsets:
            neighbor_key = ((ix + dx) % nx, (iy + dy) % ny, (iz + dz) % nz)
            if neighbor_key in seen_neighbor_keys:
                continue
            seen_neighbor_keys.add(neighbor_key)
            neighbor_members = cells.get(neighbor_key)
            if not neighbor_members:
                continue
            for first in members:
                for second in neighbor_members:
                    count_pair(first, second)

    return b_rows, coordination


def coordination_brute_force(
    b_rows: Sequence[AtomRow],
    bounds: Bounds3,
    cutoff: float,
) -> List[int]:
    lengths = tuple(upper - lower for lower, upper in bounds)
    lowers = tuple(bound[0] for bound in bounds)
    wrapped = [
        tuple((row.position[axis] - lowers[axis]) % lengths[axis] for axis in range(3))
        for row in b_rows
    ]
    coordination = [0] * len(b_rows)
    cutoff_squared = cutoff * cutoff
    half_lengths = tuple(length / 2.0 for length in lengths)

    for first in range(len(b_rows)):
        x1, y1, z1 = wrapped[first]
        for second in range(first + 1, len(b_rows)):
            x2, y2, z2 = wrapped[second]
            dx = abs(x2 - x1)
            if dx > half_lengths[0]:
                dx = lengths[0] - dx
            if dx > cutoff:
                continue
            dy = abs(y2 - y1)
            if dy > half_lengths[1]:
                dy = lengths[1] - dy
            if dy > cutoff:
                continue
            dz = abs(z2 - z1)
            if dz > half_lengths[2]:
                dz = lengths[2] - dz
            if dz > cutoff:
                continue
            if dx * dx + dy * dy + dz * dz <= cutoff_squared:
                coordination[first] += 1
                coordination[second] += 1
    return coordination


def replace_type_token(row: AtomRow, new_type: int) -> str:
    start, end = row.type_span
    return f"{row.raw[:start]}{new_type}{row.raw[end:]}"


def histogram_text(coordination: Sequence[int]) -> str:
    histogram = Counter(coordination)
    return ";".join(f"{number}:{histogram[number]}" for number in sorted(histogram))


def output_temp_path(path: Path) -> Path:
    return path.with_name(f"{path.name}.tmp")


def validate_paths(dump_path: Path, out_path: Path, summary_path: Path) -> None:
    if not dump_path.is_file():
        raise FileNotFoundError(f"input dump does not exist: {dump_path}")
    resolved = [dump_path.resolve(), out_path.resolve(), summary_path.resolve()]
    if len(set(resolved)) != 3:
        raise ValueError("input dump, output dump, and summary CSV must be different paths")
    for path in (out_path, summary_path, output_temp_path(out_path), output_temp_path(summary_path)):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite existing file: {path}")


def process_dump(
    dump_path: Path,
    out_path: Path,
    summary_path: Path,
    cutoff: float,
    max_frames: int | None,
    verify_first_frame: bool,
) -> int:
    if cutoff <= 0.0:
        raise ValueError("--cutoff must be positive")
    if max_frames is not None and max_frames <= 0:
        raise ValueError("--max-frames must be a positive integer")
    validate_paths(dump_path, out_path, summary_path)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    out_temp = output_temp_path(out_path)
    summary_temp = output_temp_path(summary_path)
    started = time.monotonic()
    frames_written = 0

    try:
        with out_temp.open("x", encoding="utf-8", newline="") as dump_out, summary_temp.open(
            "x", encoding="utf-8", newline=""
        ) as summary_out:
            summary_writer = csv.writer(summary_out)
            summary_writer.writerow(
                (
                    "frame",
                    "timestep",
                    "n_atoms",
                    "n_fe",
                    "n_b",
                    "n_c",
                    "bb_coordination_histogram",
                )
            )

            for frame in iter_dump_frames(dump_path):
                if max_frames is not None and frame.index >= max_frames:
                    break

                b_rows, coordination = coordination_cell_list(frame, cutoff)
                if verify_first_frame and frame.index == 0:
                    brute_force = coordination_brute_force(b_rows, frame.bounds, cutoff)
                    if coordination != brute_force:
                        differences = sum(
                            first != second for first, second in zip(coordination, brute_force)
                        )
                        raise AssertionError(
                            f"first-frame cell-list/brute-force mismatch for {differences} B atoms"
                        )
                    print("first-frame cell-list result matches brute force", file=sys.stderr)

                c_ids = {
                    row.atom_id
                    for row, coordination_number in zip(b_rows, coordination)
                    if coordination_number == 6
                }
                n_fe = sum(row.atom_type == FE_TYPE for row in frame.atom_rows)
                n_original_b = len(b_rows)
                n_c = len(c_ids)
                n_b = n_original_b - n_c
                if n_fe + n_b + n_c != frame.natoms:
                    raise AssertionError(
                        f"frame {frame.index}: output type counts do not sum to {frame.natoms}"
                    )

                dump_out.writelines(frame.header_lines)
                for row in frame.atom_rows:
                    if row.atom_type == B_TYPE and row.atom_id in c_ids:
                        dump_out.write(replace_type_token(row, C_TYPE))
                    else:
                        dump_out.write(row.raw)

                summary_writer.writerow(
                    (
                        frame.index,
                        frame.timestep,
                        frame.natoms,
                        n_fe,
                        n_b,
                        n_c,
                        histogram_text(coordination),
                    )
                )
                frames_written += 1
                if frames_written == 1 or frames_written % 100 == 0:
                    elapsed = time.monotonic() - started
                    print(
                        f"processed {frames_written} frame(s); timestep={frame.timestep}; "
                        f"type3={n_c}; elapsed={elapsed:.1f}s",
                        file=sys.stderr,
                    )

        os.replace(out_temp, out_path)
        os.replace(summary_temp, summary_path)
    except Exception:
        for temp_path in (out_temp, summary_temp):
            try:
                temp_path.unlink()
            except FileNotFoundError:
                pass
        raise

    return frames_written


def main() -> int:
    args = parse_args()
    try:
        frames_written = process_dump(
            args.dump,
            args.out,
            args.summary,
            args.cutoff,
            args.max_frames,
            args.verify_first_frame,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Wrote {frames_written} frame(s) to {args.out}")
    print(f"Wrote per-frame counts to {args.summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
