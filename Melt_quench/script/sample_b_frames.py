#!/usr/bin/env python3
"""Sample B-network frames from a MD trajectory folder.

This script reuses the get_transport/get_avE logic to read a LAMMPS
trajectory folder, then selects 10 dump frames from the end at ~5 ps
intervals and computes the B-B bond matrix and B coordination for each
selected frame using a fixed cutoff (default 2.37 Å).

Outputs are written into one folder:
  <input_folder>/<output_name>/
"""

from __future__ import annotations

import argparse
import csv
import os
import time
from pathlib import Path

import numpy as np

from get_avE import parse_lammps_input, timestep_to_fs
from get_transport import parse_dump_positions, parse_element_labels_from_input, parse_type_labels


DEFAULT_CUTOFF = 2.37
TARGET_FRAME_COUNT = 10
TARGET_INTERVAL_PS = 5.0


def number_or_na(value):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "N/A"
    return f"{value:.6f}"


def format_time_tag(time_ps):
    return f"t{time_ps:.3f}ps".replace(".", "p")


def find_file_upwards(start_dir, filename, max_levels=3):
    current = start_dir.resolve()
    for _ in range(max_levels + 1):
        candidate = current / filename
        if candidate.exists():
            return candidate
        if current.parent == current:
            break
        current = current.parent
    return None


def list_target_dirs(root):
    targets = sorted([p for p in root.iterdir() if p.is_dir() and p.name.endswith("_dir")])
    return targets if targets else [root]


def parse_folder(folder, dump_name, input_name):
    input_path = find_file_upwards(folder, input_name)
    if input_path is None:
        raise FileNotFoundError(f"Missing input file: {folder / input_name}")
    settings = parse_lammps_input(str(input_path), dump_name)

    dump_file = settings.get("dump_file", dump_name)
    dump_path = find_file_upwards(folder, dump_file)
    if dump_path is None and dump_file != dump_name:
        dump_path = find_file_upwards(folder, dump_name)
    if dump_path is None:
        raise FileNotFoundError(f"Missing dump file under: {folder}")

    timesteps, _, types, pos_u, pos_w, boxes, volumes = parse_dump_positions(str(dump_path), stride=1, max_frames=None)
    timesteps = np.asarray(timesteps, dtype=float)
    if len(timesteps) < TARGET_FRAME_COUNT:
        raise RuntimeError(f"Need at least {TARGET_FRAME_COUNT} dump frames, got {len(timesteps)}")

    dt_step_fs = timestep_to_fs(settings["units"], settings["timestep"])
    time_ps = np.asarray((timesteps - timesteps[0]) * (dt_step_fs * 1.0e-3), dtype=float)

    return settings, dump_path, timesteps, types, pos_w, boxes, time_ps


def select_frame_indices(time_ps, count=TARGET_FRAME_COUNT, interval_ps=TARGET_INTERVAL_PS):
    selected = []
    cursor = len(time_ps) - 1
    for i in range(count):
        target = float(time_ps[-1] - i * interval_ps)
        candidates = np.arange(0, cursor + 1)
        if candidates.size == 0:
            raise RuntimeError("Not enough frames to satisfy the requested sampling pattern")
        local = candidates[np.argmin(np.abs(time_ps[candidates] - target))]
        selected.append(int(local))
        cursor = local - 1
    selected = sorted(set(selected))
    if len(selected) != count:
        raise RuntimeError(
            f"Could not select {count} unique frames at {interval_ps} ps spacing; got {len(selected)}"
        )
    return selected


def parse_b_types(input_path, settings, types, type_labels_text=None):
    unique_types = np.array(sorted(set(int(t) for t in types)), dtype=int)
    label_map = parse_element_labels_from_input(str(input_path), settings)
    type_labels, label_to_type = parse_type_labels(type_labels_text or "", unique_types, label_map)
    b_type = None
    for label, atom_type in label_to_type.items():
        if label.upper() == "B":
            b_type = int(atom_type)
            break
    if b_type is None:
        raise RuntimeError(
            f"Cannot identify B atom type from labels: {type_labels}. Use --type-labels if needed."
        )
    return b_type, type_labels


def wrap_minimum_image(delta, box_lengths):
    return delta - np.round(delta / box_lengths) * box_lengths


def build_b_bonds(frame_pos, box_lengths, b_mask, cutoff):
    b_pos = frame_pos[b_mask]
    n_b = len(b_pos)
    adjacency = np.zeros((n_b, n_b), dtype=np.uint8)
    distances = []
    for i in range(n_b):
        delta = b_pos[i + 1 :] - b_pos[i]
        delta = wrap_minimum_image(delta, box_lengths)
        dist = np.linalg.norm(delta, axis=1)
        hit = np.where(dist <= cutoff)[0]
        for offset in hit:
            j = i + 1 + int(offset)
            adjacency[i, j] = 1
            adjacency[j, i] = 1
            distances.append((i + 1, j + 1, float(dist[offset])))
    coordination = adjacency.sum(axis=1).astype(int)
    return b_pos, adjacency, coordination, distances


def write_matrix_csv(path, adjacency):
    n = adjacency.shape[0]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["B_index"] + [f"B{i + 1}" for i in range(n)])
        for i in range(n):
            writer.writerow([f"B{i + 1}"] + [int(x) for x in adjacency[i]])


def write_coordination_csv(path, b_positions, coordination):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["B_index", "x_A", "y_A", "z_A", "coordination_B"])
        for i, xyz in enumerate(b_positions):
            writer.writerow([i + 1, f"{xyz[0]:.6f}", f"{xyz[1]:.6f}", f"{xyz[2]:.6f}", int(coordination[i])])


def write_edges_csv(path, edges):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["B_index_i", "B_index_j", "distance_A"])
        writer.writerows(edges)


def main():
    parser = argparse.ArgumentParser(description="Sample 10 dump frames from the end and compute B-B bond matrices.")
    parser.add_argument("input_dir", type=Path, help="Parent folder containing *_dir temperature subfolders, or a single temperature folder")
    parser.add_argument("--dump", default="SUPERCELL.dump", help="Dump file name (default: SUPERCELL.dump)")
    parser.add_argument("--input", default="in.lammps", help="LAMMPS input file name (default: in.lammps)")
    parser.add_argument("--output-dir", type=Path, default=None, help="Output folder (default: <input_dir>/b_frame_analysis)")
    parser.add_argument("--cutoff", type=float, default=DEFAULT_CUTOFF, help="B-B cutoff in Angstrom")
    parser.add_argument(
        "--type-labels",
        default=None,
        help="Optional override for atom labels, e.g. 1:Fe,2:B; used to identify B type",
    )
    args = parser.parse_args()

    folder = args.input_dir.resolve()
    output_root = (args.output_dir.resolve() if args.output_dir is not None else folder / "b_frame_analysis")
    output_root.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    overall_rows = []
    for target_dir in list_target_dirs(folder):
        t_dir0 = time.perf_counter()
        settings, dump_path, timesteps, types, pos_w, boxes, time_ps = parse_folder(target_dir, args.dump, args.input)
        input_path = find_file_upwards(target_dir, args.input)
        b_type, type_labels = parse_b_types(input_path, settings, types, args.type_labels)
        b_mask = np.array(types == b_type)
        selected = select_frame_indices(time_ps, TARGET_FRAME_COUNT, TARGET_INTERVAL_PS)
        t_select = time.perf_counter()

        target_output = output_root / target_dir.name
        target_output.mkdir(parents=True, exist_ok=True)
        selection_rows = []
        summary_rows = []
        for frame_no, frame_index in enumerate(selected, start=1):
            frame_time = float(time_ps[frame_index])
            frame_tag = format_time_tag(frame_time)
            frame_dir = target_output / f"frame_{frame_no:02d}_{frame_tag}"
            frame_dir.mkdir(exist_ok=True)

            b_positions, adjacency, coordination, edges = build_b_bonds(
                pos_w[frame_index], boxes[frame_index], b_mask, args.cutoff
            )

            write_matrix_csv(frame_dir / "b_bond_matrix.csv", adjacency)
            write_coordination_csv(frame_dir / "b_coordination.csv", b_positions, coordination)
            write_edges_csv(frame_dir / "b_bond_edges.csv", edges)

            selection_rows.append((frame_no, frame_index + 1, int(timesteps[frame_index]), frame_time))
            summary_rows.append(
                {
                    "temperature_dir": target_dir.name,
                    "frame_no": frame_no,
                    "dump_frame_index": frame_index + 1,
                    "timestep": int(timesteps[frame_index]),
                    "time_ps": frame_time,
                    "b_atoms": int(len(b_positions)),
                    "b_bonds": int(adjacency.sum() // 2),
                    "avg_coordination": float(np.mean(coordination)),
                    "min_coordination": int(np.min(coordination)),
                    "max_coordination": int(np.max(coordination)),
                    "output_dir": str(frame_dir),
                }
            )

        with open(target_output / "selected_frames.csv", "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["frame_no", "dump_frame_index", "timestep", "time_ps"])
            writer.writerows(selection_rows)

        with open(target_output / "summary.csv", "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
            writer.writeheader()
            writer.writerows(summary_rows)

        total_dir_time = time.perf_counter() - t_dir0
        with open(target_output / "summary.txt", "w", encoding="utf-8") as handle:
            handle.write(f"Input folder: {target_dir}\n")
            handle.write(f"Dump file: {dump_path}\n")
            handle.write(f"B type: {b_type} ({type_labels.get(b_type, 'B')})\n")
            handle.write(f"B-B cutoff (A): {args.cutoff:.3f}\n")
            handle.write(f"Selected frames: {TARGET_FRAME_COUNT} every {TARGET_INTERVAL_PS:.1f} ps from the end\n")
            handle.write("\nSelected frame list:\n")
            for row in summary_rows:
                handle.write(
                    f"  frame {row['frame_no']:02d}: timestep={row['timestep']} time_ps={row['time_ps']:.6f} "
                    f"bonds={row['b_bonds']} avg_coord={row['avg_coordination']:.6f}\n"
                )
            handle.write("\nTiming (s):\n")
            handle.write(f"  parse_and_select: {t_select - t_dir0:.6f}\n")
            handle.write(f"  total: {total_dir_time:.6f}\n")

        overall_rows.extend(summary_rows)
        print(f"Processed: {target_dir}")
        print(f"  Output folder: {target_output}")
        print(f"  B type: {b_type} ({type_labels.get(b_type, 'B')})")
        for row in summary_rows:
            print(
                f"    frame {row['frame_no']:02d}: time={row['time_ps']:.6f} ps, "
                f"bonds={row['b_bonds']}, avg_coord={row['avg_coordination']:.6f}"
            )
        print(f"  Timing (s): parse_and_select={t_select - t_dir0:.6f}, total={total_dir_time:.6f}")

    with open(output_root / "overall_summary.csv", "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(overall_rows[0].keys()))
        writer.writeheader()
        writer.writerows(overall_rows)

    print(f"All outputs written to: {output_root}")
    print(f"Total elapsed (s): {time.perf_counter() - t0:.6f}")


if __name__ == "__main__":
    main()
