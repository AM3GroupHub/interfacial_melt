#!/usr/bin/env python3
"""Compute average MSD for six atom-layer groups selected from a VASP structure."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, TextIO, Tuple

from calc_msd_selected_atoms import Vector3, iter_lammps_frames


SCRIPT_DIR = Path(__file__).resolve().parent
MD_DIR = SCRIPT_DIR.parents[2] / "two_phase_md"
DEFAULT_DUMP = MD_DIR / "SUPERCELL.dump"
DEFAULT_VASP = MD_DIR / "optimized_tot.vasp"
DEFAULT_OUT_DIR = SCRIPT_DIR / "regenerated_data"
DEFAULT_Z_TOLERANCE = 0.010
EXPECTED_ATOMS_PER_LAYER = 24

TYPE_BY_ELEMENT = {"Fe": 1, "B": 2}


@dataclass(frozen=True)
class LayerGroup:
    name: str
    system: str
    region: str
    element: str
    target_z: Tuple[float, float]

    @property
    def expected_type(self) -> int:
        return TYPE_BY_ELEMENT[self.element]

    @property
    def output_name(self) -> str:
        return f"{self.name}_layer_msd.csv"


@dataclass(frozen=True)
class StructureAtom:
    atom_id: int
    element: str
    fractional_z: float


@dataclass(frozen=True)
class LayerSelection:
    atom_id: int
    element: str
    fractional_z: float
    target_z: float
    delta_z: float


GROUPS = (
    LayerGroup("FeB_interface_B", "FeB", "interface", "B", (0.48506, 0.86386)),
    LayerGroup("FeB_interface_Fe", "FeB", "interface", "Fe", (0.47349, 0.87591)),
    LayerGroup("FeB_bulk_B", "FeB", "bulk", "B", (0.64905, 0.70184)),
    LayerGroup("FeB_bulk_Fe", "FeB", "bulk", "Fe", (0.66051, 0.69036)),
    LayerGroup("B_interface_B", "B", "interface", "B", (0.45169, 0.90469)),
    LayerGroup("B_bulk_B", "B", "bulk", "B", (0.10627, 0.16710)),
)

MSD_HEADER = (
    "frame",
    "timestep",
    "group",
    "system",
    "region",
    "element",
    "n_atoms",
    "target_z",
    "z_tolerance",
    "msd",
    "msd_x",
    "msd_y",
    "msd_z",
    "drift_x",
    "drift_y",
    "drift_z",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dump", type=Path, default=DEFAULT_DUMP, help="LAMMPS dump file")
    parser.add_argument("--vasp", type=Path, default=DEFAULT_VASP, help="VASP structure file")
    parser.add_argument(
        "--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="directory for output CSV files"
    )
    parser.add_argument(
        "--z-tolerance",
        type=float,
        default=DEFAULT_Z_TOLERANCE,
        help="maximum wrapped fractional-z distance from a target layer",
    )
    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="read at most this many frames, useful for a smoke test",
    )
    return parser.parse_args()


def parse_vector(line: str, context: str) -> Vector3:
    parts = line.split()
    if len(parts) < 3:
        raise ValueError(f"{context}: expected at least three numbers")
    try:
        return float(parts[0]), float(parts[1]), float(parts[2])
    except ValueError as exc:
        raise ValueError(f"{context}: invalid numeric vector {line!r}") from exc


def cross(a: Vector3, b: Vector3) -> Vector3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def dot(a: Vector3, b: Vector3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cartesian_to_fractional(
    position: Vector3,
    lattice: Tuple[Vector3, Vector3, Vector3],
) -> Vector3:
    a, b, c = lattice
    volume = dot(a, cross(b, c))
    if abs(volume) < 1e-15:
        raise ValueError("VASP lattice matrix is singular")
    return (
        dot(position, cross(b, c)) / volume,
        dot(position, cross(c, a)) / volume,
        dot(position, cross(a, b)) / volume,
    )


def read_vasp_atoms(path: Path) -> List[StructureAtom]:
    if not path.is_file():
        raise FileNotFoundError(f"VASP file does not exist: {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 9:
        raise ValueError(f"{path}: incomplete VASP structure")

    scale_tokens = lines[1].split()
    if len(scale_tokens) != 1:
        raise ValueError(f"{path}: only a single positive VASP scale factor is supported")
    try:
        scale = float(scale_tokens[0])
    except ValueError as exc:
        raise ValueError(f"{path}: invalid VASP scale factor {lines[1]!r}") from exc
    if scale <= 0.0:
        raise ValueError(f"{path}: VASP scale factor must be positive")

    raw_lattice = tuple(parse_vector(lines[index], f"{path}:{index + 1}") for index in range(2, 5))
    lattice = tuple(
        (vector[0] * scale, vector[1] * scale, vector[2] * scale)
        for vector in raw_lattice
    )

    elements = lines[5].split()
    try:
        counts = [int(token) for token in lines[6].split()]
    except ValueError as exc:
        raise ValueError(f"{path}: invalid element counts on line 7") from exc
    if len(elements) != len(counts):
        raise ValueError(f"{path}: element and count fields have different lengths")
    unknown = [element for element in elements if element not in TYPE_BY_ELEMENT]
    if unknown:
        raise ValueError(f"{path}: unsupported element(s): {', '.join(unknown)}")

    mode_index = 7
    if lines[mode_index].strip().lower().startswith("s"):
        mode_index += 1
    mode = lines[mode_index].strip().lower()
    direct = mode.startswith("d")
    cartesian = mode.startswith("c") or mode.startswith("k")
    if not direct and not cartesian:
        raise ValueError(f"{path}:{mode_index + 1}: unknown coordinate mode {lines[mode_index]!r}")

    atom_elements: List[str] = []
    for element, count in zip(elements, counts):
        atom_elements.extend([element] * count)

    coordinate_start = mode_index + 1
    if len(lines) < coordinate_start + len(atom_elements):
        raise ValueError(f"{path}: fewer coordinate rows than declared atoms")

    atoms: List[StructureAtom] = []
    for offset, element in enumerate(atom_elements):
        vector = parse_vector(
            lines[coordinate_start + offset],
            f"{path}:{coordinate_start + offset + 1}",
        )
        if direct:
            fractional = vector
        else:
            scaled_position = (vector[0] * scale, vector[1] * scale, vector[2] * scale)
            fractional = cartesian_to_fractional(scaled_position, lattice)
        atoms.append(
            StructureAtom(
                atom_id=offset + 1,
                element=element,
                fractional_z=fractional[2] % 1.0,
            )
        )
    return atoms


def wrapped_z_distance(a: float, b: float) -> float:
    difference = abs(a - b) % 1.0
    return min(difference, 1.0 - difference)


def select_layer_groups(
    atoms: Sequence[StructureAtom],
    tolerance: float,
) -> Dict[str, List[LayerSelection]]:
    if tolerance <= 0.0 or tolerance >= 0.5:
        raise ValueError("--z-tolerance must be between 0 and 0.5")

    selections: Dict[str, List[LayerSelection]] = {}
    atom_memberships: Dict[int, str] = {}
    for group in GROUPS:
        selected: List[LayerSelection] = []
        per_layer_counts = [0, 0]
        for atom in atoms:
            if atom.element != group.element:
                continue
            distances = [wrapped_z_distance(atom.fractional_z, z) for z in group.target_z]
            layer_index = 0 if distances[0] <= distances[1] else 1
            distance = distances[layer_index]
            if distance <= tolerance:
                selected.append(
                    LayerSelection(
                        atom_id=atom.atom_id,
                        element=atom.element,
                        fractional_z=atom.fractional_z,
                        target_z=group.target_z[layer_index],
                        delta_z=distance,
                    )
                )
                per_layer_counts[layer_index] += 1

        expected = [EXPECTED_ATOMS_PER_LAYER, EXPECTED_ATOMS_PER_LAYER]
        if per_layer_counts != expected:
            raise ValueError(
                f"{group.name}: selected {per_layer_counts[0]} and {per_layer_counts[1]} atoms "
                f"near z={group.target_z}; expected {expected[0]} per layer. "
                "Check the VASP file or --z-tolerance."
            )

        for selection in selected:
            previous = atom_memberships.get(selection.atom_id)
            if previous is not None:
                raise ValueError(
                    f"atom id {selection.atom_id} belongs to both {previous} and {group.name}"
                )
            atom_memberships[selection.atom_id] = group.name
        selections[group.name] = sorted(selected, key=lambda item: item.atom_id)
    return selections


def write_selection_csv(
    out_path: Path,
    selections: Dict[str, List[LayerSelection]],
    tolerance: float,
) -> None:
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            (
                "group",
                "system",
                "region",
                "element",
                "atom_id",
                "fractional_z",
                "target_z",
                "delta_z",
                "z_tolerance",
            )
        )
        for group in GROUPS:
            for selection in selections[group.name]:
                writer.writerow(
                    (
                        group.name,
                        group.system,
                        group.region,
                        group.element,
                        selection.atom_id,
                        f"{selection.fractional_z:.10g}",
                        f"{selection.target_z:.10g}",
                        f"{selection.delta_z:.10g}",
                        f"{tolerance:.10g}",
                    )
                )


def format_number(value: float) -> str:
    return f"{value:.16g}"


def write_layer_msd_files(
    dump_path: Path,
    vasp_path: Path,
    out_dir: Path,
    tolerance: float,
    max_frames: int | None,
) -> Tuple[int, List[Path], Dict[str, List[LayerSelection]]]:
    if not dump_path.is_file():
        raise FileNotFoundError(f"dump file does not exist: {dump_path}")
    if max_frames is not None and max_frames <= 0:
        raise ValueError("--max-frames must be a positive integer")

    atoms = read_vasp_atoms(vasp_path)
    selections = select_layer_groups(atoms, tolerance)
    selected_ids = {
        selection.atom_id
        for group_selections in selections.values()
        for selection in group_selections
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    selection_path = out_dir / "layer_atom_ids.csv"
    write_selection_csv(selection_path, selections, tolerance)
    output_paths = [out_dir / group.output_name for group in GROUPS]

    reference_positions: Dict[int, Vector3] | None = None
    reference_com: Vector3 | None = None
    rows_written = 0

    with ExitStack() as stack:
        handles: List[TextIO] = [
            stack.enter_context(path.open("w", encoding="utf-8", newline=""))
            for path in output_paths
        ]
        writers = [csv.writer(handle) for handle in handles]
        for writer in writers:
            writer.writerow(MSD_HEADER)

        for frame in iter_lammps_frames(dump_path, selected_ids):
            if max_frames is not None and frame.index >= max_frames:
                break

            if reference_positions is None:
                for group in GROUPS:
                    for selection in selections[group.name]:
                        actual_type = frame.selected_types[selection.atom_id]
                        if actual_type != group.expected_type:
                            raise ValueError(
                                f"atom id {selection.atom_id} is {group.element} in VASP but has "
                                f"dump type {actual_type}, expected type {group.expected_type}"
                            )
                reference_positions = frame.selected_positions.copy()
                reference_com = frame.com

            assert reference_positions is not None
            assert reference_com is not None
            drift = (
                frame.com[0] - reference_com[0],
                frame.com[1] - reference_com[1],
                frame.com[2] - reference_com[2],
            )

            for group, writer in zip(GROUPS, writers):
                group_selections = selections[group.name]
                sum_x = 0.0
                sum_y = 0.0
                sum_z = 0.0
                for selection in group_selections:
                    position = frame.selected_positions[selection.atom_id]
                    reference = reference_positions[selection.atom_id]
                    dx = position[0] - reference[0] - drift[0]
                    dy = position[1] - reference[1] - drift[1]
                    dz = position[2] - reference[2] - drift[2]
                    sum_x += dx * dx
                    sum_y += dy * dy
                    sum_z += dz * dz

                n_atoms = len(group_selections)
                msd_x = sum_x / n_atoms
                msd_y = sum_y / n_atoms
                msd_z = sum_z / n_atoms
                writer.writerow(
                    (
                        frame.index,
                        frame.timestep,
                        group.name,
                        group.system,
                        group.region,
                        group.element,
                        n_atoms,
                        ";".join(f"{z:.5f}" for z in group.target_z),
                        f"{tolerance:.10g}",
                        format_number(msd_x + msd_y + msd_z),
                        format_number(msd_x),
                        format_number(msd_y),
                        format_number(msd_z),
                        format_number(drift[0]),
                        format_number(drift[1]),
                        format_number(drift[2]),
                    )
                )

            rows_written += 1
            if rows_written % 500 == 0:
                print(f"processed {rows_written} frame(s)", file=sys.stderr)

    return rows_written, [selection_path, *output_paths], selections


def main() -> int:
    args = parse_args()
    try:
        rows_written, output_paths, selections = write_layer_msd_files(
            args.dump,
            args.vasp,
            args.out_dir,
            args.z_tolerance,
            args.max_frames,
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    for group in GROUPS:
        print(f"{group.name}: {len(selections[group.name])} atoms")
    print(f"Wrote {rows_written} MSD row(s) to each group CSV:")
    for path in output_paths:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
