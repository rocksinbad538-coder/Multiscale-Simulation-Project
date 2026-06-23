from pathlib import Path
import argparse
import numpy as np


def read_cube(path):
    lines = path.read_text(errors="ignore").splitlines()

    nat = int(lines[2].split()[0])
    nat_abs = abs(nat)
    origin = np.array(list(map(float, lines[2].split()[1:4])))

    nx, vx = int(lines[3].split()[0]), np.array(list(map(float, lines[3].split()[1:4])))
    ny, vy = int(lines[4].split()[0]), np.array(list(map(float, lines[4].split()[1:4])))
    nz, vz = int(lines[5].split()[0]), np.array(list(map(float, lines[5].split()[1:4])))

    atoms = []
    for i in range(6, 6 + nat_abs):
        p = lines[i].split()
        atoms.append((int(float(p[0])), float(p[2]), float(p[3]), float(p[4])))

    values = []
    for line in lines[6 + nat_abs:]:
        values.extend(float(x) for x in line.split())

    expected = nx * ny * nz
    if len(values) > expected:
        values = values[:expected]
    if len(values) < expected:
        raise ValueError(f"{path}: expected {expected} cube values, found {len(values)}")

    rho = np.array(values).reshape((nx, ny, nz)) ** 2
    return atoms, rho, origin, vx, vy, vz


def infer_fragments_from_atoms(atoms):
    # Pyrene dimer convention in these hydrated-dimer cubes:
    # first 26 atoms = first chromophore, next 26 atoms = second chromophore,
    # remaining atoms = water shell.
    frag1 = np.array([[x, y, z] for _, x, y, z in atoms[:26]])
    frag2 = np.array([[x, y, z] for _, x, y, z in atoms[26:52]])
    water = np.array([[x, y, z] for _, x, y, z in atoms[52:]])
    return frag1, frag2, water


def min_dist2_to_fragment(r, frag):
    if len(frag) == 0:
        return np.inf
    d = frag - r
    return float(np.min(np.sum(d * d, axis=1)))


def classify_cube(cube_path, atom_radius=2.0):
    atoms, rho, origin, vx, vy, vz = read_cube(cube_path)
    frag1, frag2, water = infer_fragments_from_atoms(atoms)

    nx, ny, nz = rho.shape
    rcut2 = atom_radius ** 2

    sums = {"fragment1": 0.0, "fragment2": 0.0, "water": 0.0}

    for ix in range(nx):
        base_x = origin + ix * vx
        for iy in range(ny):
            base_xy = base_x + iy * vy
            for iz in range(nz):
                val = float(rho[ix, iy, iz])
                if val == 0.0:
                    continue

                r = base_xy + iz * vz

                d1 = min_dist2_to_fragment(r, frag1)
                d2 = min_dist2_to_fragment(r, frag2)
                dw = min_dist2_to_fragment(r, water)

                dmin = min(d1, d2, dw)

                # radial mask: assign only density within atom_radius of any atom.
                # If outside all masks, assign to nearest fragment to conserve norm.
                if dmin > rcut2:
                    if d1 <= d2 and d1 <= dw:
                        sums["fragment1"] += val
                    elif d2 <= d1 and d2 <= dw:
                        sums["fragment2"] += val
                    else:
                        sums["water"] += val
                else:
                    if d1 <= d2 and d1 <= dw:
                        sums["fragment1"] += val
                    elif d2 <= d1 and d2 <= dw:
                        sums["fragment2"] += val
                    else:
                        sums["water"] += val

    total = sum(sums.values())
    if total <= 0:
        raise ValueError(f"{cube_path}: zero cube density integral")

    return {k: v / total for k, v in sums.items()}


def find_cube(cube_dir, mo):
    matches = sorted(cube_dir.glob(f"*.mo{mo}.cube"))
    if not matches:
        raise FileNotFoundError(f"No cube file found for MO {mo} in {cube_dir}")
    if len(matches) > 1:
        raise RuntimeError(f"Multiple cube files found for MO {mo}: {matches}")
    return matches[0]


def infer_pair_label(cube_dir):
    name = cube_dir.name
    # Expected: PYR2_PYR3_water0p50_orbitals
    parts = name.split("_")
    pyrs = [p for p in parts if p.startswith("PYR")]
    if len(pyrs) >= 2:
        return pyrs[0], pyrs[1]
    return "fragment1", "fragment2"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cube-dir", required=True)
    ap.add_argument("--mos", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--atom-radius", type=float, default=2.0)
    args = ap.parse_args()

    cube_dir = Path(args.cube_dir)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    label1, label2 = infer_pair_label(cube_dir)

    rows = []
    for mo in args.mos:
        cube = find_cube(cube_dir, mo)
        frac = classify_cube(cube, atom_radius=args.atom_radius)
        rows.append((mo, frac["fragment1"], frac["fragment2"], frac["water"]))

    with out.open("w") as f:
        f.write(f"MO,{label1}_fraction,{label2}_fraction,WATER_fraction\n")
        for mo, f1, f2, fw in rows:
            f.write(f"{mo},{f1:.6f},{f2:.6f},{fw:.6f}\n")

    print(out)
    print(out.read_text())


if __name__ == "__main__":
    main()
