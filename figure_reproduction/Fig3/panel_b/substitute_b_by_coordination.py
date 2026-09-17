import argparse
import csv
from pathlib import Path


def read_poscar(path):
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    species = lines[5].split()
    counts = [int(value) for value in lines[6].split()]

    coord_mode_index = 7
    if lines[7].strip().lower().startswith("s"):
        coord_mode_index = 8
    coord_mode = lines[coord_mode_index].strip()
    start = coord_mode_index + 1
    total_atoms = sum(counts)
    coords = [lines[start + i].rstrip() for i in range(total_atoms)]
    prefix = lines[:start]
    return prefix, species, counts, coord_mode, coords


def parse_mapping(text):
    mapping = {}
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        left, right = item.split(":", 1)
        mapping[int(left.strip())] = right.strip()
    return mapping


def read_substitution_map(path, mapping, only_zero_3ring=False):
    substitutions = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            coord = int(row["coordination_B"])
            poscar_index = int(row["POSCAR_index"])
            if only_zero_3ring and int(row["n_ring_3"]) != 0:
                continue
            if coord in mapping:
                substitutions[poscar_index] = mapping[coord]
    return substitutions


def species_per_atom(species, counts):
    out = []
    for symbol, count in zip(species, counts):
        out.extend([symbol] * count)
    return out


def regroup_coordinates(symbols, coords, species_order):
    grouped = {symbol: [] for symbol in species_order}
    for symbol, coord in zip(symbols, coords):
        grouped[symbol].append(coord)
    counts = [len(grouped[symbol]) for symbol in species_order]
    ordered_coords = []
    for symbol in species_order:
        ordered_coords.extend(grouped[symbol])
    return counts, ordered_coords


def write_poscar(path, prefix, species_order, counts, coord_mode, coords):
    lines = prefix[:]
    lines[5] = " ".join(species_order)
    lines[6] = " ".join(str(count) for count in counts)
    lines.extend(coords)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_species_order(original_species, mapping):
    species_order = list(original_species)
    for new_symbol in mapping.values():
        if new_symbol not in species_order:
            species_order.append(new_symbol)
    return species_order


def main():
    parser = argparse.ArgumentParser(description="Replace 3-coordinated B with N and 4-coordinated B with C in a POSCAR.")
    parser.add_argument("--poscar", type=Path, required=True, help="Input POSCAR path")
    parser.add_argument("--coordination-csv", type=Path, required=True, help="Per-atom coordination CSV")
    parser.add_argument("--output", type=Path, required=True, help="Output POSCAR path")
    parser.add_argument(
        "--mapping",
        type=str,
        default="3:N,4:C",
        help="Replacement rules as coordination:symbol pairs, e.g. '3:N,4:C' or '6:C'",
    )
    parser.add_argument(
        "--only-zero-3ring",
        action="store_true",
        help="Only replace B atoms with n_ring_3 = 0",
    )
    args = parser.parse_args()

    prefix, species, counts, coord_mode, coords = read_poscar(args.poscar)
    mapping = parse_mapping(args.mapping)
    substitutions = read_substitution_map(args.coordination_csv, mapping, only_zero_3ring=args.only_zero_3ring)
    atom_symbols = species_per_atom(species, counts)

    for poscar_index, new_symbol in substitutions.items():
        atom_symbols[poscar_index - 1] = new_symbol

    species_order = build_species_order(species, mapping)
    new_counts, new_coords = regroup_coordinates(atom_symbols, coords, species_order)
    write_poscar(args.output, prefix, species_order, new_counts, coord_mode, new_coords)

    n_c = sum(1 for symbol in atom_symbols if symbol == "C")
    n_n = sum(1 for symbol in atom_symbols if symbol == "N")
    n_b = sum(1 for symbol in atom_symbols if symbol == "B")
    print(f"Wrote {args.output}")
    print(f"Mapping: {mapping}")
    print(f"Total C: {n_c}")
    print(f"Total N: {n_n}")
    print(f"Remaining B: {n_b}")


if __name__ == "__main__":
    main()
