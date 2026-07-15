from pathlib import Path
import argparse
import numpy as np


def parse_cube_atom_lines(lines, nat_abs, origin, vx, vy, vz, nx, ny, nz):
    grid_min = origin
    grid_max = origin + (nx - 1) * vx + (ny - 1) * vy + (nz - 1) * vz
    lo = np.minimum(grid_min, grid_max) - 5.0
    hi = np.maximum(grid_min, grid_max) + 5.0

    candidates = []
    for mode in ("gaussian_5col", "orca_4or5col"):
        coords = []
        ok = True
        for i in range(6, 6 + nat_abs):
            p = lines[i].split()
            try:
                if mode == "gaussian_5col" and len(p) >= 5:
                    xyz = list(map(float, p[2:5]))
                else:
                    xyz = list(map(float, p[1:4]))
                coords.append(xyz)
            except Exception:
                ok = False
                break
        if not ok:
            continue
        arr = np.array(coords)
        score = np.mean(np.all((arr >= lo) & (arr <= hi), axis=1))
        candidates.append((score, arr))

    if not candidates:
        raise ValueError("Could not parse cube atom coordinates.")

    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]


def read_cube(path):
    lines = path.read_text(errors="ignore").splitlines()

    nat = int(lines[2].split()[0])
    nat_abs = abs(nat)
    origin = np.array(list(map(float, lines[2].split()[1:4])))

    nx, vx = int(lines[3].split()[0]), np.array(list(map(float, lines[3].split()[1:4])))
    ny, vy = int(lines[4].split()[0]), np.array(list(map(float, lines[4].split()[1:4])))
    nz, vz = int(lines[5].split()[0]), np.array(list(map(float, lines[5].split()[1:4])))

    atoms = parse_cube_atom_lines(lines, nat_abs, origin, vx, vy, vz, nx, ny, nz)

    values = []
    for line in lines[6 + nat_abs:]:
        values.extend(float(x) for x in line.split())

    expected = nx * ny * nz
    values = values[-expected:] if len(values) >= expected else values
    if len(values) != expected:
        raise ValueError(f"{path}: expected {expected} cube values, got {len(values)}")

    rho = np.array(values).reshape((nx, ny, nz)) ** 2
    return atoms, rho, origin, vx, vy, vz


def classify_cube(cube_path):
    atoms, rho, origin, vx, vy, vz = read_cube(cube_path)

    frag1 = atoms[:26]
    frag2 = atoms[26:52]
    water = atoms[52:]

    nx, ny, nz = rho.shape
    sums = {"fragment1": 0.0, "fragment2": 0.0, "water": 0.0}

    for ix in range(nx):
        rx = origin + ix * vx
        for iy in range(ny):
            rxy = rx + iy * vy
            for iz in range(nz):
                val = float(rho[ix, iy, iz])
                if val == 0.0:
                    continue
                r = rxy + iz * vz

                d1 = np.min(np.sum((frag1 - r) ** 2, axis=1))
                d2 = np.min(np.sum((frag2 - r) ** 2, axis=1))
                dw = np.min(np.sum((water - r) ** 2, axis=1)) if len(water) else np.inf

                if d1 <= d2 and d1 <= dw:
                    sums["fragment1"] += val
                elif d2 <= d1 and d2 <= dw:
                    sums["fragment2"] += val
                else:
                    sums["water"] += val

    total = sum(sums.values())
    return {k: v / total for k, v in sums.items()}


def infer_labels(cube_dir):
    parts = cube_dir.name.split("_")
    pyrs = [p for p in parts if p.startswith("PYR")]
    if len(pyrs) >= 2:
        return pyrs[0], pyrs[1]
    return "FRAG1", "FRAG2"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cube-dir", required=True)
    ap.add_argument("--mos", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cube_dir = Path(args.cube_dir)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    label1, label2 = infer_labels(cube_dir)

    with out.open("w") as f:
        f.write(f"MO,{label1}_fraction,{label2}_fraction,WATER_fraction\n")
        for mo in args.mos:
            matches = sorted(cube_dir.glob(f"*.mo{mo}.cube"))
            if len(matches) != 1:
                raise RuntimeError(f"Expected one cube for MO {mo}, found: {matches}")
            frac = classify_cube(matches[0])
            f.write(
                f"{mo},{frac['fragment1']:.6f},{frac['fragment2']:.6f},{frac['water']:.6f}\n"
            )

    print(out)
    print(out.read_text())


if __name__ == "__main__":
    main()
