#!/usr/bin/env python3
import numpy as np
import argparse
import sys

PLANCK_eVs = 4.135667696e-15  # Planck constant [eV*s]

def parse_lammps_dump_velocity(filename):
    """
    Parse LAMMPS dump with 'ITEM:' headers.
    Returns:
      timesteps: list[int]
      ids_ref: (n_atoms,)
      types_ref: (n_atoms,)
      vel: (n_frames, n_atoms, 3) in the reference id order
    Automatically detects velocity columns: vx/vy/vz or v_x/v_y/v_z.
    """
    timesteps = []
    frames_vel = []
    ids_ref = None
    types_ref = None

    with open(filename, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue

            # TIMESTEP
            ts_line = f.readline()
            if not ts_line:
                break
            ts = int(ts_line.strip())
            timesteps.append(ts)

            # NUMBER OF ATOMS
            hdr = f.readline()
            if not hdr or not hdr.startswith("ITEM: NUMBER OF ATOMS"):
                raise RuntimeError("Malformed dump: missing 'ITEM: NUMBER OF ATOMS'")
            n_atoms = int(f.readline().strip())

            # BOX BOUNDS
            hdr = f.readline()
            if not hdr or not hdr.startswith("ITEM: BOX BOUNDS"):
                raise RuntimeError("Malformed dump: missing 'ITEM: BOX BOUNDS'")
            _ = f.readline(); _ = f.readline(); _ = f.readline()

            # ATOMS header
            hdr = f.readline()
            if not hdr or not hdr.startswith("ITEM: ATOMS"):
                raise RuntimeError("Malformed dump: missing 'ITEM: ATOMS'")
            cols = hdr.strip().split()[2:]
            col_to_idx = {name: i for i, name in enumerate(cols)}

            # Required columns
            if 'id' not in col_to_idx or 'type' not in col_to_idx:
                raise RuntimeError("Dump must contain 'id' and 'type' columns.")
            # Velocity columns
            vx_name = next((c for c in ('vx', 'v_x') if c in col_to_idx), None)
            vy_name = next((c for c in ('vy', 'v_y') if c in col_to_idx), None)
            vz_name = next((c for c in ('vz', 'v_z') if c in col_to_idx), None)
            if vx_name is None or vy_name is None or vz_name is None:
                raise RuntimeError(f"No velocity columns found. Got columns: {cols}")

            # Read atoms block
            ids_frame = np.empty(n_atoms, dtype=int)
            types_frame = np.empty(n_atoms, dtype=int)
            vel_frame = np.empty((n_atoms, 3), dtype=float)
            for i in range(n_atoms):
                parts = f.readline().strip().split()
                if len(parts) < len(cols):
                    raise RuntimeError(f"Atom line has fewer fields ({len(parts)}) than columns ({len(cols)}).")
                ids_frame[i] = int(parts[col_to_idx['id']])
                types_frame[i] = int(parts[col_to_idx['type']])
                vx = float(parts[col_to_idx[vx_name]])
                vy = float(parts[col_to_idx[vy_name]])
                vz = float(parts[col_to_idx[vz_name]])
                vel_frame[i, :] = (vx, vy, vz)

            if ids_ref is None:
                ids_ref = ids_frame.copy()
                types_ref = types_frame.copy()
                frames_vel.append(vel_frame)
            else:
                # Align current frame velocities to reference id order
                id_to_idx = {int(i): k for k, i in enumerate(ids_frame)}
                reorder = np.array([id_to_idx[int(i)] for i in ids_ref], dtype=int)
                frames_vel.append(vel_frame[reorder])

    if len(frames_vel) == 0:
        raise RuntimeError("No frames parsed from dump.")

    vel = np.stack(frames_vel, axis=0)
    return timesteps, ids_ref, types_ref, vel

def build_masses(types, mass_map_str=None):
    """
    Build mass array per atom type. If None, all ones.
    mass_map_str example: '1:6.94,2:35.45'
    """
    n_atoms = len(types)
    masses = np.ones(n_atoms, dtype=float)
    if mass_map_str:
        m = {}
        for e in mass_map_str.split(','):
            e = e.strip()
            if not e: continue
            k, v = e.split(':')
            m[int(k)] = float(v)
        for i in range(n_atoms):
            t = int(types[i])
            if t in m:
                masses[i] = m[t]
            else:
                raise ValueError(f"Atom type {t} missing in mass map.")
    return masses

def remove_com_drift(vel, masses):
    """
    Remove mass-weighted COM drift per frame.
    vel: (n_frames, n_atoms, 3)
    masses: (n_atoms,)
    """
    m = masses.reshape(1, -1, 1)
    v_com = (vel * m).sum(axis=1, keepdims=True) / m.sum(axis=1, keepdims=True)
    return vel - v_com

def compute_vdos(vel, dt_s, masses=None, use_hann=True):
    """
    Compute VDOS on energy axis (eV) with unit area.
    vel: (n_frames, n_atoms, 3) after COM removal
    dt_s: frame spacing in seconds
    masses: (n_atoms,), optional weighting
    use_hann: apply Hann window along time
    Returns:
      E_eV: (n_freq,)
      vdos: (n_freq,) normalized such that integral over E is 1
    """
    n_frames, n_atoms, _ = vel.shape
    if use_hann:
        w = np.hanning(n_frames).reshape(n_frames, 1, 1)
        vel_w = vel * w
    else:
        vel_w = vel

    # rFFT along time
    Vx = np.fft.rfft(vel_w[:, :, 0], axis=0)
    Vy = np.fft.rfft(vel_w[:, :, 1], axis=0)
    Vz = np.fft.rfft(vel_w[:, :, 2], axis=0)

    # Power spectrum per atom, sum over components
    P = (np.abs(Vx)**2 + np.abs(Vy)**2 + np.abs(Vz)**2)  # (n_freq, n_atoms)
    if masses is not None:
        P = P * masses.reshape(1, -1)

    # Average over atoms
    P_avg = P.mean(axis=1)  # (n_freq,)

    # Frequency axis and one-sided PSD scaling to per-Hz density
    freqs = np.fft.rfftfreq(n_frames, d=dt_s)  # Hz
    fs = 1.0 / dt_s
    S_f = P_avg / (n_frames * fs)  # per Hz
    # Double interior bins for one-sided spectrum
    if n_frames % 2 == 0:
        S_f[1:-1] *= 2.0
    else:
        S_f[1:] *= 2.0

    # Convert to energy axis: E = h f
    E_eV = PLANCK_eVs * freqs
    # Density per eV: S_E = S_f * (df/dE) = S_f / h
    S_E = S_f / PLANCK_eVs

    # Area normalize over energy axis
    area = np.trapz(S_E, E_eV)
    vdos = S_E / area if area > 0 else S_E

    return E_eV, vdos

def main():
    parser = argparse.ArgumentParser(description="Compute VDOS (energy axis in eV, area-normalized) from LAMMPS dump with velocities.")
    parser.add_argument("dump", help="LAMMPS dump file (text) with ATOMS including velocities")
    parser.add_argument("--dt-fs", type=float, default=10.0, help="Frame spacing in fs (e.g., output every 10 steps with 1 fs timestep -> 10)")
    parser.add_argument("--mass", type=str, default=None, help="Mass map 'type:mass,...' (optional), e.g. '1:6.94,2:35.45'")
    parser.add_argument("--no-window", action="store_true", help="Disable Hann window")
    parser.add_argument("--out", type=str, default="vdos_eV.dat", help="Output file (E[eV], VDOS)")
    args = parser.parse_args()

    dt_s = args.dt_fs * 1e-15

    try:
        timesteps, ids_ref, types_ref, vel = parse_lammps_dump_velocity(args.dump)
    except Exception as e:
        print(f"Error parsing dump: {e}", file=sys.stderr)
        sys.exit(1)

    n_frames, n_atoms, _ = vel.shape
    if n_frames < 2:
        print("Need at least 2 frames to compute spectrum.", file=sys.stderr)
        sys.exit(1)

    masses = build_masses(types_ref, args.mass)
    vel_corr = remove_com_drift(vel, masses)

    E_eV, vdos = compute_vdos(vel_corr, dt_s, masses=masses, use_hann=not args.no_window)
    np.savetxt(args.out, np.column_stack([E_eV, vdos]), header="E_eV  VDOS_normalized", comments='')

    print(f"Frames: {n_frames}, Atoms: {n_atoms}")
    print(f"dt = {args.dt_fs} fs -> Nyquist energy ~ {E_eV[-1]:.6f} eV")
    print(f"Saved: {args.out}")

if __name__ == "__main__":
    main()

