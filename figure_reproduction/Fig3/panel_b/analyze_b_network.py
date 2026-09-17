import argparse
import csv
import time
from collections import Counter, defaultdict, deque
from pathlib import Path

import numpy as np


DEFAULT_CUTOFF = 2.37
MIN_RING = 3
MAX_RING = 8


def parse_poscar(path):
    lines = [line.rstrip() for line in path.read_text(encoding="utf-8", errors="ignore").splitlines() if line.strip()]
    scale = float(lines[1].split()[0])
    lattice = np.array([[float(value) for value in lines[index].split()] for index in range(2, 5)], dtype=float) * scale
    species = lines[5].split()
    counts = [int(value) for value in lines[6].split()]

    coord_mode_index = 7
    if lines[7].lower().startswith("s"):
        coord_mode_index = 8
    coord_mode = lines[coord_mode_index].lower()
    start = coord_mode_index + 1
    total_atoms = sum(counts)
    coords = np.array([[float(value) for value in lines[start + i].split()[:3]] for i in range(total_atoms)], dtype=float)

    if coord_mode.startswith("d"):
        frac = coords
        cart = frac @ lattice
    else:
        cart = coords
        frac = cart @ np.linalg.inv(lattice)

    species_per_atom = []
    poscar_indices = []
    offset = 0
    for symbol, count in zip(species, counts):
        for local_index in range(count):
            species_per_atom.append(symbol)
            poscar_indices.append(offset + local_index + 1)
        offset += count

    return lattice, frac, cart, species_per_atom, poscar_indices


def extract_b_atoms(lattice, frac, cart, species_per_atom, poscar_indices):
    indices = [index for index, symbol in enumerate(species_per_atom) if symbol == "B"]
    return {
        "lattice": lattice,
        "frac": frac[indices],
        "cart": cart[indices],
        "poscar_indices": [poscar_indices[index] for index in indices],
    }


def build_b_graph(frac, lattice, cutoff):
    n_atoms = len(frac)
    adjacency = [set() for _ in range(n_atoms)]
    edges = []
    distances = {}
    for i in range(n_atoms):
        delta_frac = frac[i + 1 :] - frac[i]
        delta_frac -= np.round(delta_frac)
        delta_cart = delta_frac @ lattice
        delta_dist = np.linalg.norm(delta_cart, axis=1)
        for offset in np.where(delta_dist <= cutoff)[0]:
            j = i + 1 + int(offset)
            distance = float(delta_dist[offset])
            adjacency[i].add(j)
            adjacency[j].add(i)
            distances[(i, j)] = distance
            edges.append((i, j, distance))
    return adjacency, edges, distances


def canonical_cycle(nodes):
    cycle = list(nodes)
    n_nodes = len(cycle)
    rotations = [tuple(cycle[i:] + cycle[:i]) for i in range(n_nodes)]
    reversed_cycle = list(reversed(cycle))
    rotations.extend(tuple(reversed_cycle[i:] + reversed_cycle[:i]) for i in range(n_nodes))
    return min(rotations)


def shortest_paths_between_edge_endpoints(adjacency, start, end, max_depth):
    distances = {start: 0}
    parents = defaultdict(list)
    queue = deque([start])

    while queue:
        node = queue.popleft()
        depth = distances[node]
        if depth >= max_depth:
            continue
        for neighbor in adjacency[node]:
            if (node == start and neighbor == end) or (node == end and neighbor == start):
                continue
            next_depth = depth + 1
            if neighbor not in distances:
                distances[neighbor] = next_depth
                parents[neighbor].append(node)
                queue.append(neighbor)
            elif distances[neighbor] == next_depth:
                parents[neighbor].append(node)

    if end not in distances:
        return None, []

    shortest_depth = distances[end]
    paths = []
    stack = [(end, [end])]
    while stack:
        node, reversed_path = stack.pop()
        if node == start:
            paths.append(list(reversed(reversed_path)))
            continue
        for parent in parents[node]:
            stack.append((parent, reversed_path + [parent]))
    return shortest_depth, paths


def enumerate_shortest_path_rings(adjacency, min_ring, max_ring):
    rings = set()
    for start, neighbors in enumerate(adjacency):
        for end in neighbors:
            if end <= start:
                continue
            shortest_depth, paths = shortest_paths_between_edge_endpoints(adjacency, start, end, max_ring - 1)
            if shortest_depth is None:
                continue
            ring_size = shortest_depth + 1
            if ring_size < min_ring or ring_size > max_ring:
                continue
            for path in paths:
                rings.add(canonical_cycle(path))
    return sorted(rings, key=lambda cycle: (len(cycle), cycle))


def per_atom_ring_counts(n_atoms, rings, min_ring, max_ring):
    counts = {size: np.zeros(n_atoms, dtype=int) for size in range(min_ring, max_ring + 1)}
    ring_size_counter = Counter()
    for ring in rings:
        ring_size = len(ring)
        ring_size_counter[ring_size] += 1
        for node in ring:
            counts[ring_size][node] += 1
    return counts, ring_size_counter


def write_adjacency_matrix(path, adjacency):
    n_atoms = len(adjacency)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        header = [""] + [f"B{i + 1}" for i in range(n_atoms)]
        writer.writerow(header)
        for i in range(n_atoms):
            row = [f"B{i + 1}"]
            row.extend(1 if j in adjacency[i] else 0 for j in range(n_atoms))
            writer.writerow(row)


def write_edge_list(path, edges, poscar_indices):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["B_index_i", "B_index_j", "POSCAR_index_i", "POSCAR_index_j", "distance_A"])
        for i, j, distance in edges:
            writer.writerow([i + 1, j + 1, poscar_indices[i], poscar_indices[j], f"{distance:.6f}"])


def write_per_atom_table(path, cart, poscar_indices, coordination, ring_counts, min_ring, max_ring):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        header = ["B_index", "POSCAR_index", "x_A", "y_A", "z_A", "coordination_B"]
        header.extend(f"n_ring_{size}" for size in range(min_ring, max_ring + 1))
        writer.writerow(header)
        for index, xyz in enumerate(cart):
            row = [index + 1, poscar_indices[index], f"{xyz[0]:.6f}", f"{xyz[1]:.6f}", f"{xyz[2]:.6f}", coordination[index]]
            row.extend(int(ring_counts[size][index]) for size in range(min_ring, max_ring + 1))
            writer.writerow(row)


def write_ring_list(path, rings, poscar_indices):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ring_size", "B_indices_1based", "POSCAR_indices_1based"])
        for ring in rings:
            writer.writerow([
                len(ring),
                " ".join(str(node + 1) for node in ring),
                " ".join(str(poscar_indices[node]) for node in ring),
            ])


def write_summary(path, poscar_path, cutoff, n_atoms, edges, coordination, ring_size_counter, timings):
    coordination_counter = Counter(int(value) for value in coordination)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f"POSCAR: {poscar_path}\n")
        handle.write(f"B-B cutoff (A): {cutoff:.3f}\n")
        handle.write(f"B atoms: {n_atoms}\n")
        handle.write(f"B-B bonds: {len(edges)}\n")
        handle.write(f"Average B coordination: {np.mean(coordination):.6f}\n")
        handle.write(f"Min/Max B coordination: {np.min(coordination)} / {np.max(coordination)}\n")
        handle.write("\nB coordination histogram:\n")
        for coordination_number in sorted(coordination_counter):
            handle.write(f"  coordination {coordination_number}: {coordination_counter[coordination_number]} B atoms\n")
        handle.write("\nUnique primitive shortest-path ring counts:\n")
        for size in range(MIN_RING, MAX_RING + 1):
            handle.write(f"  {size}-member rings: {ring_size_counter.get(size, 0)}\n")
        handle.write("\nTiming (seconds):\n")
        for key, value in timings.items():
            handle.write(f"  {key}: {value:.6f}\n")


def output_path(base_dir, stem, label, suffix):
    if label:
        return base_dir / f"{stem}_{label}{suffix}"
    return base_dir / f"{stem}{suffix}"


def main():
    parser = argparse.ArgumentParser(description="Analyze the B-only network in a POSCAR using a fixed B-B cutoff.")
    parser.add_argument("--poscar", type=Path, default=Path("POSCAR"), help="Path to POSCAR")
    parser.add_argument("--cutoff", type=float, default=DEFAULT_CUTOFF, help="B-B cutoff in Angstrom")
    parser.add_argument("--label", type=str, default="", help="Optional label suffix for output files")
    args = parser.parse_args()

    t0 = time.perf_counter()
    lattice, frac, cart, species_per_atom, poscar_indices = parse_poscar(args.poscar)
    b_data = extract_b_atoms(lattice, frac, cart, species_per_atom, poscar_indices)
    t_parse = time.perf_counter()

    adjacency, edges, distances = build_b_graph(b_data["frac"], b_data["lattice"], args.cutoff)
    coordination = np.array([len(neighbors) for neighbors in adjacency], dtype=int)
    t_graph = time.perf_counter()

    rings = enumerate_shortest_path_rings(adjacency, MIN_RING, MAX_RING)
    ring_counts, ring_size_counter = per_atom_ring_counts(len(b_data["frac"]), rings, MIN_RING, MAX_RING)
    t_rings = time.perf_counter()

    out_dir = args.poscar.resolve().parent
    write_adjacency_matrix(output_path(out_dir, "b_bond_matrix", args.label, ".csv"), adjacency)
    write_edge_list(output_path(out_dir, "b_bond_edges", args.label, ".csv"), edges, b_data["poscar_indices"])
    write_per_atom_table(
        output_path(out_dir, "b_network_per_atom", args.label, ".csv"),
        b_data["cart"],
        b_data["poscar_indices"],
        coordination,
        ring_counts,
        MIN_RING,
        MAX_RING,
    )
    write_ring_list(output_path(out_dir, "b_unique_rings", args.label, ".csv"), rings, b_data["poscar_indices"])
    timings = {
        "parse_poscar": t_parse - t0,
        "build_bond_graph": t_graph - t_parse,
        "enumerate_rings": t_rings - t_graph,
        "total": t_rings - t0,
    }
    write_summary(
        output_path(out_dir, "b_ring_summary", args.label, ".txt"),
        args.poscar.resolve(),
        args.cutoff,
        len(b_data["frac"]),
        edges,
        coordination,
        ring_size_counter,
        timings,
    )

    print(f"B atoms: {len(b_data['frac'])}")
    print(f"B-B bonds: {len(edges)}")
    for size in range(MIN_RING, MAX_RING + 1):
        print(f"{size}-member rings: {ring_size_counter.get(size, 0)}")
    print("Timing (s):")
    for key, value in timings.items():
        print(f"  {key}: {value:.6f}")
    print(f"Outputs written to: {out_dir}")


if __name__ == "__main__":
    main()
