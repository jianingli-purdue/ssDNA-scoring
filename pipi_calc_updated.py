#!/usr/bin/env python
"""
Count π–π stacking between nucleobase rings without relying on aromatic typing.
Adds:
  --strict-T        : vector alignment (true edge-to-face)
  --strict-T-line   : centroid vector must pass through face ring centroid (lateral tol)
  --lateral-tol     : max lateral offset of projected centroid on face ring plane (Å)

Usage:
  $SCHRODINGER/run pipi_calc.py input.(pdb|mae|maegz)
      [--chains A,B] [--csv stacks.csv] [--print-pairs]
      [--dmax 5.5] [--amax 30] [--tmin 60] [--tmax 120] [--tmaxdist 6.0] [--tprox 3.6]
      [--strict-T] [--strict-T-line] [--lateral-tol 1.2]
      [--strict-align-parallel 0.8] [--strict-align-perp 0.2]
"""

import argparse
import math
import csv
from collections import defaultdict
import numpy as np

from schrodinger.structure import StructureReader
from schrodinger import adapter

# SMARTS that do NOT require aromatic typing; match C/N 6-membered cycles
# NOTE: We intentionally do NOT match 5-member rings so that purines are
# treated as a single π system instead of two fused rings.
SMARTS_RING6_CN = "[#6,#7]1[#6,#7][#6,#7][#6,#7][#6,#7][#6,#7]1"

# Sugar anchor atom names that identify nucleic-acid residues robustly
SUGAR_ANCHORS = {"C1'", "C1*", "O4'", "O4*"}

# ----------------- Geometry helpers -----------------

def unit(v):
    n = np.linalg.norm(v)
    return v if n == 0 else v / n

def plane_fit(points):
    """Fit plane by PCA (SVD). Return (centroid, unit normal)."""
    pts = np.asarray(points, dtype=float)
    ctr = pts.mean(axis=0)
    X = pts - ctr
    _, _, vt = np.linalg.svd(X, full_matrices=False)
    normal = vt[-1]
    return ctr, unit(normal)

def angle_between_normals(n1, n2):
    """Angle in degrees between two normals; 0 == parallel (up to inversion)."""
    c = np.clip(np.dot(unit(n1), unit(n2)), -1.0, 1.0)
    ang = math.degrees(math.acos(c))
    return min(ang, 180.0 - ang)

def point_to_plane_distance(p, p0, n):
    return abs(np.dot((np.asarray(p) - np.asarray(p0)), unit(n)))

def lateral_offset_of_projection(c_face, n_face, c_edge):
    """
    Distance (in-plane) between the face ring centroid and the projection
    of the edge ring centroid onto the face ring plane.
    """
    r = c_edge - c_face
    n = unit(n_face)
    # subtract normal component -> in-plane vector
    in_plane = r - np.dot(r, n) * n
    return np.linalg.norm(in_plane)

def residue_key(a):
    """Tuple identifying a residue."""
    chain = getattr(a, "chain", "") or ""
    resnum = getattr(a, "resnum", None)
    ins = getattr(a, "insertion_code", "") or ""
    rname = (getattr(a, "pdbres", None) or getattr(a, "resname", "")).upper()
    return (rname, chain, resnum, ins)

def short_rid(rk):
    rname, chain, resnum, ins = rk
    return f"{rname}:{chain}{resnum or ''}{ins}"

# ----------------- Ring perception -----------------

def smart_matches(st, smarts, unique_sets=True):
    try:
        return adapter.evaluate_smarts(st, smarts, unique_sets=unique_sets) or []
    except Exception as e:
        raise RuntimeError(
            f"SMARTS evaluation failed for pattern '{smarts}'. "
            f"Use schrodinger.adapter.evaluate_smarts (2024+)."
        ) from e

def residue_has_sugar_anchor(atoms_in_res):
    for a in atoms_in_res:
        an = (getattr(a, "pdbname", None) or getattr(a, "name", "")).strip()
        if an in SUGAR_ANCHORS:
            return True
    return False

def find_aromatic_like_rings_in_nucleic_acids(st, allowed_chains=None):
    """
    Return list of ring dicts regardless of aromatic typing:
      {"rk": residue_key, "rid": 'XYZ:A15', "size": 6, "centroid": ..., "normal": ..., "coords": ...}

    We only keep 6-membered rings so that purines (fused 5+6 systems)
    are treated as a single π system represented by the 6-member ring.
    """
    idx_to_atom = {a.index: a for a in st.atom}

    res_to_atoms = defaultdict(list)
    for a in st.atom:
        res_to_atoms[residue_key(a)].append(a)

    def collect(smarts, size):
        rings = []
        for match in smart_matches(st, smarts, unique_sets=True):
            atoms = [idx_to_atom[i] for i in match]
            rks = {residue_key(a) for a in atoms}
            if len(rks) != 1:
                continue
            rk = next(iter(rks))
            rname, chain, _, _ = rk
            if allowed_chains and (chain not in allowed_chains):
                continue
            if not residue_has_sugar_anchor(res_to_atoms[rk]):
                continue
            coords = np.array([a.xyz for a in atoms], dtype=float)
            if coords.shape[0] != size:
                continue
            centroid, normal = plane_fit(coords)
            rings.append({
                "rk": rk,
                "rid": short_rid(rk),
                "size": size,
                "centroid": centroid,
                "normal": normal,
                "coords": coords
            })
        return rings

    # Only collect 6-membered rings (pyrimidines and the 6-ring of purines)
    rings6 = collect(SMARTS_RING6_CN, 6)

    # Deduplicate 6-membered rings if SMARTS produces overlaps
    seen, unique = set(), []
    for r in rings6:
        sig = (r["rid"], r["size"], tuple(np.round(r["centroid"], 3)))
        if sig not in seen:
            seen.add(sig)
            unique.append(r)
    return unique

# ----------------- Classification -----------------

def classify_ringpair(r1, r2,
                      dmax=5.5, amax=30.0,
                      tmin=60.0, tmax=120.0, tmaxdist=6.0, tprox=3.6,
                      strict_T=False, align_parallel=0.8, align_perp=0.2,
                      strict_T_line=False, lateral_tol=1.2):
    """
    Classify 'parallel' (face-to-face) or 'T' (edge-to-face), else None.

    strict_T:
        One ring's normal aligns with r12 (centroid vector), the other's is ~perpendicular:
           max(|n1·r12|, |n2·r12|) >= align_parallel
        and
           min(|n1·r12|, |n2·r12|) <= align_perp

    strict_T_line:
        Additionally require the centroid vector to pass (approximately) through
        the face ring centroid: the projection of the edge ring centroid onto the
        face ring plane must lie within 'lateral_tol' Å of the face ring centroid.
    """
    c1, n1, coords1 = r1["centroid"], r1["normal"], r1["coords"]
    c2, n2, coords2 = r2["centroid"], r2["normal"], r2["coords"]

    dcent = np.linalg.norm(c1 - c2)
    ang = angle_between_normals(n1, n2)

    # Parallel
    if dcent <= dmax and ang <= amax:
        return "parallel", {"d": dcent, "ang": ang, "prox12": None, "prox21": None, "lateral": None, "face": None}

    # Candidate T
    if not (tmin <= ang <= tmax and dcent <= tmaxdist):
        return None, {"d": dcent, "ang": ang, "prox12": None, "prox21": None, "lateral": None, "face": None}

    # Vector alignment (true-T) if requested
    r12 = unit(c2 - c1)
    if strict_T:
        a1 = abs(np.dot(unit(n1), r12))
        a2 = abs(np.dot(unit(n2), r12))
        # two possible orientations:
        # (A) ring1 face (aligned), ring2 edge (perp)
        face1 = (a1 >= align_parallel and a2 <= align_perp)
        # (B) ring2 face (aligned), ring1 edge (perp)
        face2 = (a2 >= align_parallel and a1 <= align_perp)
        if not (face1 or face2):
            return None, {"d": dcent, "ang": ang, "prox12": None, "prox21": None, "lateral": None, "face": None}
    else:
        # if not strict_T, allow either orientation
        face1 = True
        face2 = True

    # Proximity (edge atoms near face plane)
    prox12 = min(point_to_plane_distance(p, c2, n2) for p in coords1)
    prox21 = min(point_to_plane_distance(p, c1, n1) for p in coords2)
    if min(prox12, prox21) > tprox:
        return None, {"d": dcent, "ang": ang, "prox12": prox12, "prox21": prox21, "lateral": None, "face": None}

    # Centroid-in-line refinement
    if strict_T_line:
        # If ring1 is face, edge centroid (c2) projected onto plane1 must land near c1
        lat1 = lateral_offset_of_projection(c1, n1, c2) if face1 else np.inf
        # If ring2 is face, edge centroid (c1) projected onto plane2 must land near c2
        lat2 = lateral_offset_of_projection(c2, n2, c1) if face2 else np.inf

        passes1 = face1 and (lat1 <= lateral_tol)
        passes2 = face2 and (lat2 <= lateral_tol)

        if not (passes1 or passes2):
            return None, {
                "d": dcent,
                "ang": ang,
                "prox12": prox12,
                "prox21": prox21,
                "lateral": min(lat1, lat2),
                "face": None
            }

        # pick orientation that passes
        if passes1 and (lat1 <= lat2):
            lateral = lat1
            face = 1
        elif passes2:
            lateral = lat2
            face = 2
        else:
            lateral = min(lat1, lat2)
            face = 1 if lat1 <= lat2 else 2
    else:
        lateral, face = None, None

    return "T", {
        "d": dcent,
        "ang": ang,
        "prox12": prox12,
        "prox21": prox21,
        "lateral": lateral,
        "face": face
    }

# ----------------- Main -----------------

def main():
    ap = argparse.ArgumentParser(
        description="Count π–π stacks between nucleobase rings with strict T and centroid-line options."
    )
    ap.add_argument("input", help="Input structure (PDB/MAE/MAEGZ).")
    ap.add_argument("--chains", help="Comma-separated chain IDs to include (e.g., A,B).")
    ap.add_argument("--csv", help="Write interactions to CSV.")
    ap.add_argument("--print-pairs", action="store_true", help="Print detailed pairs.")

    # geometry cutoffs
    ap.add_argument("--dmax", type=float, default=5.5, help="Parallel centroid distance cutoff (Å).")
    ap.add_argument("--amax", type=float, default=30.0, help="Parallel angle cutoff (deg).")
    ap.add_argument("--tmin", type=float, default=60.0, help="T-shaped min angle (deg).")
    ap.add_argument("--tmax", type=float, default=120.0, help="T-shaped max angle (deg).")
    ap.add_argument("--tmaxdist", type=float, default=6.0, help="T-shaped centroid distance cutoff (Å).")
    ap.add_argument("--tprox", type=float, default=3.6, help="T-shaped min atom-to-plane proximity (Å).")

    # strict T options
    ap.add_argument("--strict-T", action="store_true", help="Require true edge-to-face vector alignment for T-type.")
    ap.add_argument(
        "--strict-align-parallel",
        type=float,
        default=0.8,
        help="|n·r12| threshold for the 'face' ring (default 0.8 ≈ 36°)."
    )
    ap.add_argument(
        "--strict-align-perp",
        type=float,
        default=0.2,
        help="|n·r12| threshold for the 'edge' ring (default 0.2 ≈ 78°)."
    )

    # centroid-line option
    ap.add_argument(
        "--strict-T-line",
        action="store_true",
        help="Additionally require centroid vector to pass through the face ring centroid (lateral tol)."
    )
    ap.add_argument(
        "--lateral-tol",
        type=float,
        default=1.2,
        help="Max lateral offset (Å) of projected centroid on face ring plane (default 1.2 Å)."
    )

    args = ap.parse_args()

    chains = None
    if args.chains:
        chains = set(c.strip() for c in args.chains.split(",") if c.strip())

    # Read first structure from file
    with StructureReader(args.input) as rdr:
        st = next(rdr)

    # Find nucleic-acid rings
    rings = find_aromatic_like_rings_in_nucleic_acids(st, allowed_chains=chains)
    if not rings:
        print(
            "No nucleic-acid rings found. If your sugar atom names are unusual, "
            "add them to SUGAR_ANCHORS or remove the anchor filter."
        )
        return

    # Group rings by residue; we only consider inter-residue pairs.
    by_res = defaultdict(list)
    for r in rings:
        by_res[r["rid"]].append(r)

    rows = []
    keys = list(by_res.keys())
    n_parallel = 0
    n_t = 0

    # NOTE: Using j in range(i+1, ...) ensures each residue pair is only
    # considered once (e.g., 5–6 but not 6–5).
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            ri, rj = keys[i], keys[j]
            for ring_i in by_res[ri]:
                for ring_j in by_res[rj]:
                    label, metrics = classify_ringpair(
                        ring_i,
                        ring_j,
                        dmax=args.dmax,
                        amax=args.amax,
                        tmin=args.tmin,
                        tmax=args.tmax,
                        tmaxdist=args.tmaxdist,
                        tprox=args.tprox,
                        strict_T=args.strict_T,
                        align_parallel=args.strict_align_parallel,
                        align_perp=args.strict_align_perp,
                        strict_T_line=args.strict_T_line,
                        lateral_tol=args.lateral_tol,
                    )
                    if label:
                        row = {
                            "type": label,
                            "res1": ri,
                            "ring1": f"{ring_i['size']}-member",
                            "res2": rj,
                            "ring2": f"{ring_j['size']}-member",
                            "dist_cc": round(metrics["d"], 3),
                            "angle_normals_deg": round(metrics["ang"], 2),
                            "prox12_planeÅ": round(metrics["prox12"], 3)
                            if metrics["prox12"] is not None
                            else "",
                            "prox21_planeÅ": round(metrics["prox21"], 3)
                            if metrics["prox21"] is not None
                            else "",
                        }
                        rows.append(row)
                        if label == "parallel":
                            n_parallel += 1
                        else:
                            n_t += 1
                        if args.print_pairs:
                            extra = ""
                            if args.strict_T or args.strict_T_line:
                                extra = (
                                    f"  (face={metrics.get('face')}, "
                                    f"lateral={None if metrics.get('lateral') is None else round(metrics['lateral'],3)} Å)"
                                )
                            print(
                                f"{label:8s}  {ri:<10s}({ring_i['size']}) — {rj:<10s}({ring_j['size']})"
                                f"  d={row['dist_cc']} Å  ang={row['angle_normals_deg']}°"
                                f"  prox12={row['prox12_planeÅ']} Å  prox21={row['prox21_planeÅ']} Å{extra}"
                            )

    print(
        f"[π–π stacking] total={len(rows)}  parallel={n_parallel}  T-shaped={n_t}"
        + (
            " (strict-T" + (", line" if args.strict_T_line else "") + ")"
            if args.strict_T or args.strict_T_line
            else ""
        )
    )

    if args.csv:
        fields = [
            "type",
            "res1",
            "ring1",
            "res2",
            "ring2",
            "dist_cc",
            "angle_normals_deg",
            "prox12_planeÅ",
            "prox21_planeÅ",
        ]
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        print(f"Wrote {len(rows)} interactions to {args.csv}")

if __name__ == "__main__":
    main()

