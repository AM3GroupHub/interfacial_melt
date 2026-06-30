#!/usr/bin/env python3
"""Density, Onsager transport, and RDF workflow for LAMMPS trajectories."""

import argparse
import os
import re
import numpy as np

from get_avE import block_average_stats, build_masses, commands_from_input, frame_dt_fs, parse_lammps_input
from get_avE import resolve_string_token

KB_EV_PER_K = 8.617333262145e-5
AMU_A3_TO_G_CM3 = 1.66053906660
A2_PER_PS_TO_CM2_PER_S = 1.0e-4
A3_TO_CM3 = 1.0e-24


def is_float(text):
    try:
        float(text.replace("D", "E"))
        return True
    except ValueError:
        return False


def norm_header(text):
    return text.strip().lower().replace("_", "")


def number_or_na(value):
    return "N/A" if value is None or np.isnan(value) else f"{value:.8g}"


def parse_element_labels_from_input(input_file, settings):
    """Read atom-type element labels from pair_coeff * * ... Element1 Element2 ..."""
    mass_map = settings.get("masses", {})
    if not mass_map:
        return {}
    n_types = max(int(atom_type) for atom_type in mass_map)
    variables = settings.get("variables", {})
    with open(input_file, "r", encoding="utf-8", errors="ignore") as handle:
        commands = commands_from_input(handle.readlines())
    for command in commands:
        parts = command.split()
        if len(parts) < 3 or parts[0].lower() != "pair_coeff":
            continue
        if parts[1] != "*" or parts[2] != "*":
            continue
        candidates = [resolve_string_token(token, variables) for token in parts[3:]]
        if len(candidates) < n_types:
            continue
        labels = candidates[-n_types:]
        if all(label.upper() == "NULL" or not is_float(label) for label in labels):
            return {i + 1: (label if label.upper() != "NULL" else f"type{i + 1}") for i, label in enumerate(labels)}
    return {}


def labels_for_atom_types(atom_types, label_map):
    labels = [label_map.get(int(atom_type), f"type{int(atom_type)}") for atom_type in atom_types]
    if len(set(labels)) == len(labels):
        return labels
    return [f"{label}_type{int(atom_type)}" for label, atom_type in zip(labels, atom_types)]


def parse_log_volumes(log_file):
    volumes = []
    in_thermo = False
    columns = {}
    with open(log_file, "r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            parts = line.split()
            lowered = [norm_header(p) for p in parts]
            if "step" in lowered and any(name in lowered for name in ("vol", "volume")):
                columns = {norm_header(name): i for i, name in enumerate(parts)}
                in_thermo = True
                continue
            if line.startswith(("Loop time", "WARNING", "ERROR")):
                in_thermo = False
                columns = {}
                continue
            if not in_thermo or not parts or not is_float(parts[0]):
                continue
            vol_i = columns.get("vol", columns.get("volume"))
            if vol_i is not None and vol_i < len(parts):
                try:
                    volumes.append(float(parts[vol_i].replace("D", "E")))
                except ValueError:
                    pass
    return volumes


def density_stats(log_file, total_mass_amu, skip_lines, block_size):
    volumes = parse_log_volumes(log_file)
    densities = [total_mass_amu * AMU_A3_TO_G_CM3 / v for v in volumes if v > 0]
    mean, err, samples, blocks = block_average_stats(densities, skip_lines, block_size)
    return mean, err, samples, blocks, len(volumes)


def read_box(handle):
    header = handle.readline()
    if not header.startswith("ITEM: BOX BOUNDS"):
        raise RuntimeError("Malformed dump: missing ITEM: BOX BOUNDS")
    bounds = []
    tilts = []
    for _ in range(3):
        values = [float(x) for x in handle.readline().split()]
        bounds.append(values[:2])
        tilts.append(values[2] if len(values) >= 3 else 0.0)
    lo = np.array([b[0] for b in bounds], dtype=float)
    hi = np.array([b[1] for b in bounds], dtype=float)
    lengths = hi - lo
    return lo, hi, lengths, np.array(tilts, dtype=float)


def detect_columns(col_to_idx, names):
    return [col_to_idx.get(name) for name in names]


def scaled_to_cart(values, lo, lengths):
    return lo + values * lengths


def parse_dump_positions(filename, stride=1, max_frames=None):
    timesteps, frames_u, frames_w, boxes, volumes = [], [], [], [], []
    ids_ref = types_ref = None
    frame_index = -1
    with open(filename, "r", encoding="utf-8", errors="ignore") as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            frame_index += 1
            ts = int(handle.readline().strip())
            handle.readline()
            n_atoms = int(handle.readline().strip())
            lo, _, lengths, tilts = read_box(handle)
            atom_header = handle.readline()
            if not atom_header.startswith("ITEM: ATOMS"):
                raise RuntimeError("Malformed dump: missing ITEM: ATOMS")
            cols = atom_header.split()[2:]
            col_to_idx = {name: i for i, name in enumerate(cols)}
            if "id" not in col_to_idx or "type" not in col_to_idx:
                raise RuntimeError("Dump must contain id and type columns")
            rows = [handle.readline().split() for _ in range(n_atoms)]
            if frame_index % stride != 0:
                continue
            if max_frames is not None and len(timesteps) >= max_frames:
                break
            ids = np.array([int(row[col_to_idx["id"]]) for row in rows], dtype=int)
            types = np.array([int(row[col_to_idx["type"]]) for row in rows], dtype=int)
            unwrapped, wrapped = positions_from_rows(rows, col_to_idx, lo, lengths)
            if ids_ref is None:
                ids_ref, types_ref = ids.copy(), types.copy()
                order = np.arange(n_atoms)
            else:
                id_to_idx = {int(atom_id): i for i, atom_id in enumerate(ids)}
                order = np.array([id_to_idx[int(atom_id)] for atom_id in ids_ref], dtype=int)
            timesteps.append(ts)
            frames_u.append(unwrapped[order])
            frames_w.append(wrapped[order])
            boxes.append(lengths.copy())
            volumes.append(float(np.prod(lengths)))
            if np.any(np.abs(tilts) > 1.0e-12):
                print("  - Warning: triclinic tilt detected; RDF uses orthogonal bounding lengths")
    if not frames_u:
        raise RuntimeError("No dump frames parsed")
    return timesteps, ids_ref, types_ref, np.stack(frames_u), np.stack(frames_w), np.stack(boxes), np.array(volumes)


def positions_from_rows(rows, col_to_idx, lo, lengths):
    n_atoms = len(rows)
    raw = lambda row, name: float(row[col_to_idx[name]].replace("D", "E"))
    if all(name in col_to_idx for name in ("xu", "yu", "zu")):
        unwrapped = np.array([[raw(row, "xu"), raw(row, "yu"), raw(row, "zu")] for row in rows])
    elif all(name in col_to_idx for name in ("xsu", "ysu", "zsu")):
        scaled = np.array([[raw(row, "xsu"), raw(row, "ysu"), raw(row, "zsu")] for row in rows])
        unwrapped = scaled_to_cart(scaled, lo, lengths)
    elif all(name in col_to_idx for name in ("x", "y", "z", "ix", "iy", "iz")):
        xyz = np.array([[raw(row, "x"), raw(row, "y"), raw(row, "z")] for row in rows])
        img = np.array([[raw(row, "ix"), raw(row, "iy"), raw(row, "iz")] for row in rows])
        unwrapped = xyz + img * lengths
    elif all(name in col_to_idx for name in ("xs", "ys", "zs", "ix", "iy", "iz")):
        scaled = np.array([[raw(row, "xs"), raw(row, "ys"), raw(row, "zs")] for row in rows])
        img = np.array([[raw(row, "ix"), raw(row, "iy"), raw(row, "iz")] for row in rows])
        unwrapped = scaled_to_cart(scaled + img, lo, lengths)
    elif all(name in col_to_idx for name in ("x", "y", "z")):
        print("  - Warning: using wrapped x/y/z as unwrapped positions for Onsager")
        unwrapped = np.array([[raw(row, "x"), raw(row, "y"), raw(row, "z")] for row in rows])
    elif all(name in col_to_idx for name in ("xs", "ys", "zs")):
        print("  - Warning: using wrapped xs/ys/zs as unwrapped positions for Onsager")
        scaled = np.array([[raw(row, "xs"), raw(row, "ys"), raw(row, "zs")] for row in rows])
        unwrapped = scaled_to_cart(scaled, lo, lengths)
    else:
        raise RuntimeError("Dump must contain xu/yu/zu or wrapped coordinates with image flags")
    wrapped = lo + np.mod(unwrapped - lo, lengths)
    return unwrapped.reshape(n_atoms, 3), wrapped.reshape(n_atoms, 3)


def fit_slope(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 2:
        return np.nan, np.nan
    coeff = np.polyfit(x, y, 1)
    pred = np.polyval(coeff, x)
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return float(coeff[0]), r2


def center_positions(pos, masses):
    m = masses.reshape(1, -1, 1)
    com = (pos * m).sum(axis=1, keepdims=True) / m.sum()
    return pos - com


def compute_onsager(pos, types, masses, volumes_a3, temperature, dt_ps, args):
    pos = center_positions(pos, masses)
    unique_types = np.array(sorted(set(int(t) for t in types)), dtype=int)
    masks = [types == atom_type for atom_type in unique_types]
    q = np.stack([pos[:, mask, :].sum(axis=1) for mask in masks], axis=1)
    max_lag = max(2, int((len(pos) - 1) * args.max_lag_fraction))
    lags = np.arange(args.corr_stride, max_lag + 1, args.corr_stride, dtype=int)
    if len(lags) < 2:
        lags = np.arange(1, len(pos), dtype=int)
    ns = len(unique_types)
    corr = np.zeros((len(lags), ns, ns), dtype=float)
    self_corr = np.zeros((len(lags), ns), dtype=float)
    for ilag, lag in enumerate(lags):
        origins = np.arange(0, len(pos) - lag, args.origin_stride, dtype=int)
        dq = q[origins + lag] - q[origins]
        for i in range(ns):
            for j in range(i, ns):
                value = np.mean(np.sum(dq[:, i, :] * dq[:, j, :], axis=1))
                corr[ilag, i, j] = corr[ilag, j, i] = value
            disp_i = pos[origins[:, None] + lag, masks[i], :] - pos[origins[:, None], masks[i], :]
            self_corr[ilag, i] = np.mean(np.sum(disp_i * disp_i, axis=(1, 2)))
    times_ps = lags * dt_ps
    fit_start, fit_end = fit_window(times_ps, args.onsager_fit_start_ps, args.onsager_fit_end_ps)
    fit_mask = (times_ps >= fit_start) & (times_ps <= fit_end)
    volume_cm3 = float(np.mean(volumes_a3) * A3_TO_CM3)
    prefactor = A2_PER_PS_TO_CM2_PER_S / (6.0 * KB_EV_PER_K * temperature * volume_cm3)
    lmat = np.zeros((ns, ns), dtype=float)
    r2 = np.zeros((ns, ns), dtype=float)
    for i in range(ns):
        for j in range(ns):
            slope, score = fit_slope(times_ps[fit_mask], corr[fit_mask, i, j])
            lmat[i, j] = slope * prefactor
            r2[i, j] = score
    lself = np.zeros(ns, dtype=float)
    dself = np.zeros(ns, dtype=float)
    r2self = np.zeros(ns, dtype=float)
    for i, mask in enumerate(masks):
        slope, score = fit_slope(times_ps[fit_mask], self_corr[fit_mask, i])
        lself[i] = slope * prefactor
        c_i = float(np.sum(mask)) / volume_cm3
        dself[i] = lself[i] * KB_EV_PER_K * temperature / c_i
        r2self[i] = score
    return unique_types, times_ps, corr, self_corr, lmat, r2, lself, dself, r2self, (fit_start, fit_end)


def fit_window(times_ps, start_ps, end_ps):
    tmax = float(times_ps[-1])
    start = 0.2 * tmax if start_ps is None else start_ps
    end = 0.8 * tmax if end_ps is None else end_ps
    if end <= start:
        raise ValueError("Invalid Onsager fit window")
    return start, end


def write_onsager(directory, labels, times, corr, self_corr, lmat, r2, lself, dself, r2self, fit_window_ps):
    header = "time_ps " + " ".join([f"C_{a}_{b}_A2" for a in labels for b in labels])
    header += " " + " ".join([f"Cself_{a}_A2" for a in labels])
    cols = [times]
    cols += [corr[:, i, j] for i in range(len(labels)) for j in range(len(labels))]
    cols += [self_corr[:, i] for i in range(len(labels))]
    np.savetxt(os.path.join(directory, "onsager_correlation_functions.dat"), np.column_stack(cols), header=header, comments="")
    with open(os.path.join(directory, "onsager_matrix.txt"), "w", encoding="utf-8") as handle:
        handle.write("# L_ij units: 1/(eV cm s), computed in center-of-mass frame\n")
        handle.write(f"# Fit window: {fit_window_ps[0]:.6g} to {fit_window_ps[1]:.6g} ps\n")
        handle.write("# labels: " + " ".join(labels) + "\n")
        handle.write("L_full\n")
        for row in lmat:
            handle.write(" ".join(f"{value:.8e}" for value in row) + "\n")
        handle.write("R2_full\n")
        for row in r2:
            handle.write(" ".join(f"{value:.6f}" for value in row) + "\n")
        handle.write("# type L_self D_self_cm2_s R2_self correlation_factor\n")
        for label, ls, ds, rs, lf in zip(labels, lself, dself, r2self, np.diag(lmat)):
            factor = lf / ls - 1.0 if ls != 0 else np.nan
            handle.write(f"{label} {ls:.8e} {ds:.8e} {rs:.6f} {factor:.8e}\n")


def parse_type_labels(label_text, unique_types, base_labels=None):
    labels = {}
    if base_labels:
        unique_names = labels_for_atom_types(unique_types, base_labels)
        labels.update({int(atom_type): name for atom_type, name in zip(unique_types, unique_names)})
    for atom_type in unique_types:
        labels.setdefault(int(atom_type), f"type{int(atom_type)}")
    if not label_text:
        return labels, {value: key for key, value in labels.items()}
    for item in label_text.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            key, value = item.split(":", 1)
        elif "=" in item:
            key, value = item.split("=", 1)
        else:
            raise ValueError("--type-labels entries must look like '1:Na,2:O'")
        labels[int(key.strip())] = value.strip()
    return labels, {value: key for key, value in labels.items()}


def resolve_type_token(token, label_to_type):
    token = token.strip()
    if re.fullmatch(r"\d+", token):
        return int(token)
    if token not in label_to_type:
        raise ValueError(f"Unknown RDF type label '{token}'")
    return label_to_type[token]


def requested_rdf_pairs(unique_types, args, label_to_type):
    if args.rdf_mode == "total":
        return None
    if args.rdf_mode == "all":
        return [(a, b) for ia, a in enumerate(unique_types) for b in unique_types[ia:]]
    if not args.rdf_pair:
        raise ValueError("--rdf-mode pair requires --rdf-pair, e.g. --rdf-pair 1:2")
    tokens = re.split(r"[:,\-]", args.rdf_pair)
    if len(tokens) != 2:
        raise ValueError("--rdf-pair must contain exactly two atom types, e.g. 1:2 or Na:O")
    a = resolve_type_token(tokens[0], label_to_type)
    b = resolve_type_token(tokens[1], label_to_type)
    if a not in unique_types or b not in unique_types:
        raise ValueError(f"RDF pair {a}:{b} not present in dump atom types {list(unique_types)}")
    return [(min(a, b), max(a, b))]


def compute_total_rdf(wrapped, boxes, centers, edges, shell, stride):
    hist = np.zeros(len(centers), dtype=float)
    norm = np.zeros(len(centers), dtype=float)
    for iframe in range(0, len(wrapped), stride):
        pos = wrapped[iframe]
        box = boxes[iframe]
        vol = float(np.prod(box))
        if len(pos) < 2:
            continue
        diff = pos[:, None, :] - pos[None, :, :]
        diff = diff[np.triu_indices(len(pos), k=1)]
        diff -= box * np.rint(diff / box)
        dist = np.linalg.norm(diff, axis=1)
        hist += np.histogram(dist, bins=edges)[0]
        norm += 0.5 * len(pos) * ((len(pos) - 1) / vol) * shell
    gr = np.divide(hist, norm, out=np.zeros_like(hist), where=norm > 0)
    return ["r_A", "g_total"], np.column_stack([centers, gr])


def compute_rdf(wrapped, types, boxes, args, type_label_map=None):
    unique_types = np.array(sorted(set(int(t) for t in types)), dtype=int)
    type_labels, label_to_type = parse_type_labels(args.type_labels, unique_types, type_label_map)
    rmax = args.rdf_rmax if args.rdf_rmax else 0.5 * float(np.min(np.mean(boxes, axis=0)))
    edges = np.arange(0.0, rmax + args.rdf_bin_width, args.rdf_bin_width)
    centers = 0.5 * (edges[:-1] + edges[1:])
    shell = 4.0 * np.pi / 3.0 * (edges[1:] ** 3 - edges[:-1] ** 3)
    if args.rdf_mode == "total":
        return compute_total_rdf(wrapped, boxes, centers, edges, shell, args.rdf_stride)

    pairs = requested_rdf_pairs(unique_types, args, label_to_type)
    hist = {pair: np.zeros(len(centers), dtype=float) for pair in pairs}
    norm = {pair: np.zeros(len(centers), dtype=float) for pair in pairs}
    for iframe in range(0, len(wrapped), args.rdf_stride):
        pos = wrapped[iframe]
        box = boxes[iframe]
        vol = float(np.prod(box))
        for a, b in pairs:
            ia = np.where(types == a)[0]
            ib = np.where(types == b)[0]
            if a == b:
                if len(ia) < 2:
                    continue
                diff = pos[ia][:, None, :] - pos[ia][None, :, :]
                triu = np.triu_indices(len(ia), k=1)
                diff = diff[triu]
                expected = 0.5 * len(ia) * ((len(ia) - 1) / vol) * shell
            else:
                diff = pos[ia][:, None, :] - pos[ib][None, :, :]
                diff = diff.reshape(-1, 3)
                expected = len(ia) * (len(ib) / vol) * shell
            diff -= box * np.rint(diff / box)
            dist = np.linalg.norm(diff, axis=1)
            hist[(a, b)] += np.histogram(dist, bins=edges)[0]
            norm[(a, b)] += expected
    gr = [centers]
    names = ["r_A"]
    for pair in pairs:
        gr.append(np.divide(hist[pair], norm[pair], out=np.zeros_like(hist[pair]), where=norm[pair] > 0))
        names.append(f"g_{type_labels[pair[0]]}_{type_labels[pair[1]]}")
    return names, np.column_stack(gr)


def process_directory(directory, args, output_root):
    input_path = os.path.join(directory, args.input)
    log_path = os.path.join(directory, args.filename)
    settings = parse_lammps_input(input_path, args.dump)
    temperature = settings["temperature"] or float(re.match(r"(\d+(?:\.\d+)?)_dir$", directory).group(1))
    result_dir = os.path.join(output_root, directory)
    os.makedirs(result_dir, exist_ok=True)
    dump_file = os.path.join(directory, settings["dump_file"])
    if not os.path.exists(dump_file):
        dump_file = os.path.join(directory, args.dump)
    timesteps, _, types, pos_u, pos_w, boxes, volumes = parse_dump_positions(dump_file, args.dump_stride, args.max_frames)
    skip_frames = args.skip if args.dump_skip_frames is None else args.dump_skip_frames
    skip_frames = min(skip_frames, len(pos_u) - 2) if len(pos_u) > 2 else 0
    masses = build_masses(types, settings["masses"])
    total_mass = float(np.sum(masses))
    density, density_err, n_samples, n_blocks, n_vol = density_stats(log_path, total_mass, args.skip, args.block_size)
    dt_ps = frame_dt_fs(settings, timesteps) * 1.0e-3 * args.dump_stride
    pos_u2, pos_w2, boxes2, volumes2 = pos_u[skip_frames:], pos_w[skip_frames:], boxes[skip_frames:], volumes[skip_frames:]
    labels, times, corr, scorr, lmat, r2, lself, dself, r2self, fwin = compute_onsager(
        pos_u2, types, masses, volumes2, temperature, dt_ps, args
    )
    element_labels = parse_element_labels_from_input(input_path, settings)
    label_names = labels_for_atom_types(labels, element_labels)
    write_onsager(result_dir, label_names, times, corr, scorr, lmat, r2, lself, dself, r2self, fwin)
    rdf_names, rdf_data = compute_rdf(pos_w2, types, boxes2, args, element_labels)
    np.savetxt(os.path.join(result_dir, args.rdf_output), rdf_data, header=" ".join(rdf_names), comments="")
    print(f"  - density: {number_or_na(density)} +/- {number_or_na(density_err)} g/cm3")
    print(f"  - Onsager/RDF outputs written in {result_dir}")
    return temperature, density, density_err, n_samples, n_blocks, args.block_size, n_vol, len(pos_u2), result_dir


def main():
    parser = argparse.ArgumentParser(description="Compute density, full Onsager matrix, and RDF for *_dir trajectories.")
    parser.add_argument("-s", "--skip", type=int, default=0, help="Skip first N log rows; also dump frames unless overridden")
    parser.add_argument("-f", "--filename", default="log.lammps", help="LAMMPS log file name")
    parser.add_argument(
        "-o",
        "--output",
        default="transport_results.txt",
        help="Summary output file; all results are written under OUTPUT_dir",
    )
    parser.add_argument("--input", default="in.lammps", help="LAMMPS input file name")
    parser.add_argument("--dump", default="SUPERCELL.dump", help="Default LAMMPS dump file name")
    parser.add_argument("--block-size", type=int, default=1000, help="Block size in log thermo rows for density error")
    parser.add_argument("--dump-skip-frames", type=int, default=None, help="Skip first N parsed dump frames")
    parser.add_argument("--dump-stride", type=int, default=1, help="Read every Nth dump frame")
    parser.add_argument("--max-frames", type=int, default=None, help="Maximum parsed dump frames")
    parser.add_argument("--corr-stride", type=int, default=50, help="Lag stride in frames for Onsager correlations")
    parser.add_argument("--origin-stride", type=int, default=10, help="Time-origin stride in frames for correlations")
    parser.add_argument("--max-lag-fraction", type=float, default=0.5, help="Maximum lag as fraction of trajectory length")
    parser.add_argument("--onsager-fit-start-ps", type=float, default=None, help="Fit start time in ps")
    parser.add_argument("--onsager-fit-end-ps", type=float, default=None, help="Fit end time in ps")
    parser.add_argument("--rdf-rmax", type=float, default=8.0, help="RDF maximum radius in Angstrom")
    parser.add_argument("--rdf-bin-width", type=float, default=0.02, help="RDF bin width in Angstrom")
    parser.add_argument("--rdf-stride", type=int, default=10, help="Use every Nth frame for RDF")
    parser.add_argument("--rdf-mode", choices=("all", "total", "pair"), default="all", help="RDF output mode")
    parser.add_argument("--rdf-pair", default=None, help="Atom pair for --rdf-mode pair, e.g. 1:2 or Na:O")
    parser.add_argument(
        "--type-labels",
        default=None,
        help="Optional override for atom labels, e.g. 1:Na,2:Si,3:B,4:O; otherwise read pair_coeff",
    )
    parser.add_argument("--rdf-output", default="rdf.dat", help="RDF output file name per directory")
    args = parser.parse_args()
    if args.block_size <= 0 or args.dump_stride <= 0 or args.corr_stride <= 0 or args.origin_stride <= 0:
        parser.error("block-size, dump-stride, corr-stride, and origin-stride must be positive")
    if args.rdf_stride <= 0 or args.rdf_bin_width <= 0:
        parser.error("rdf-stride and rdf-bin-width must be positive")
    output_root = f"{args.output}_dir"
    os.makedirs(output_root, exist_ok=True)
    summary_path = os.path.join(output_root, os.path.basename(args.output))
    results = []
    pattern = re.compile(r"(\d+(?:\.\d+)?)_dir$")
    for item in sorted(os.listdir(os.getcwd())):
        if not pattern.match(item) or not os.path.isdir(item):
            continue
        print(f"\nProcessing {item}")
        try:
            results.append(process_directory(item, args, output_root))
        except Exception as exc:
            print(f"  - Warning: failed to process {item}: {exc}")
    results.sort(key=lambda row: row[0])
    with open(summary_path, "w", encoding="utf-8") as handle:
        handle.write(
            "# Temperature(K)\tAverage_Density(g/cm3)\tDensity_err(g/cm3)\tN_samples\t"
            "N_blocks\tBlock_size\tN_log_volumes\tN_dump_frames\tOutput_Directory\n"
        )
        for row in results:
            handle.write(
                f"{row[0]:.2f}\t{number_or_na(row[1])}\t{number_or_na(row[2])}\t{row[3]}\t"
                f"{row[4]}\t{row[5]}\t{row[6]}\t{row[7]}\t{row[8]}\n"
            )
    print(f"\nSaved results directory: {output_root}")
    print(f"Saved summary: {summary_path}")


if __name__ == "__main__":
    main()
