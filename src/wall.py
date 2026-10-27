"""STEP 9.1-9.4 — Regional LV wall thickness and aorto-septal angle from one mesh's geometry.

Input: the `_geom` dict from src.lvot.measure_mesh(..., keep_geometry=True).
- Surfaces: the three closed boundary components are LV cavity (lvot STEP 2.2), epicardium (largest enclosed
  volume) and RV cavity (the other one).
- Thickness at an LV endocardial point = distance to the nearest outer surface: the RV cavity surface for septal
  points (nearer to the RV than to the epicardium), the epicardium otherwise. Points near the valve rings dropped.
- Segments: level = fraction of LV length from the mitral centre towards the apex (basal / mid / apical);
  sector = angle around the long axis, 0 = septal direction (mean direction of RV-facing points), four 90-deg
  sectors septal, anterior (aortic side), lateral, inferior. Value = median thickness of the segment's points.
- Subaortic septum: septal points 3-15 mm below the aortic ring along the aortic axis (basal septal bulge).
- Aorto-septal angle (as in echo PLAX): line fitted to septal points in the PLAX plane (aortic axis + mitral
  centre) over the upper half of the LV; angle = 180 - angle(aortic axis into the LV, septal line). Straight
  continuation = 180; a sigmoid septum gives a smaller angle.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from src.config import load_config

W_ = load_config()["wall"]
SECTORS = ("septal", "anterior", "lateral", "inferior")


def wall_geometry(G: dict) -> dict:
    import pyvista as pv

    # STEP 9.1 — epicardium and RV cavity surfaces
    reg, lvr = G["regions"], G["lv_region"]
    rid = np.asarray(reg.cell_data["RegionId"])
    comps = []
    for r in np.unique(rid):
        if r == lvr:
            continue
        s = reg.extract_cells(rid == r).extract_surface()
        comps.append((abs(s.volume), s))
    comps.sort(key=lambda c: -c[0])
    if len(comps) < 2:
        raise ValueError("need epicardium and RV cavity components")
    epi, rvc = comps[0][1], comps[1][1]

    # STEP 9.2 — distance of every LV endo point to each surface
    endo = pv.PolyData(np.asarray(G["endo"].points))
    d_epi = np.abs(np.asarray(endo.compute_implicit_distance(epi)["implicit_distance"]))
    d_rv = np.abs(np.asarray(endo.compute_implicit_distance(rvc)["implicit_distance"]))
    X = np.asarray(endo.points)
    septal = d_rv < d_epi
    thick = np.where(septal, d_rv, d_epi)

    ao, mv = G["rings"]["ao"], G["rings"]["mv"]
    ring_pts = np.vstack([r["centre"] + r["xy"] @ np.vstack([r["e1"], r["e2"]]) for r in (ao, mv)])
    keep = cKDTree(ring_pts).query(X)[0] > W_["ring_exclude_mm"]

    # STEP 9.3 — long-axis frame and segments
    mvc = mv["centre"]
    la = G["apex"] - mvc
    L = np.linalg.norm(la); la /= L
    h = (X - mvc) @ la / L
    V = (X - mvc) - np.outer((X - mvc) @ la, la)  # radial component
    s_sel = septal & keep & (h > .2) & (h < .7)
    e_s = V[s_sel].mean(0); e_s -= (e_s @ la) * la; e_s /= np.linalg.norm(e_s)
    e_t = np.cross(la, e_s)
    if (ao["centre"] - mvc) @ e_t < 0:
        e_t = -e_t  # anterior = aortic side
    phi = np.degrees(np.arctan2(V @ e_t, V @ e_s))
    sector = np.select([np.abs(phi) < 45, (phi >= 45) & (phi < 135), np.abs(phi) >= 135], ["septal", "anterior", "lateral"], "inferior")
    out = {"wall_n_endo_pts": int(len(X)), "wall_septal_frac": float(septal[keep].mean()),
           "wall_n_regions": int(len(comps) + 1)}
    for lev, (lo, hi) in W_["levels"].items():
        for sc in SECTORS:
            m = keep & (h >= lo) & (h < hi) & (sector == sc)
            out[f"wt_{lev}_{sc}"] = float(np.median(thick[m])) if m.sum() >= 10 else np.nan
    for lev in ("basal", "mid"):
        out[f"sep_lat_ratio_{lev}"] = out[f"wt_{lev}_septal"] / out[f"wt_{lev}_lateral"]
    out["wt_max_p95"] = float(np.percentile(thick[keep & (h >= .05) & (h < .95)], 95))

    # STEP 9.3b — subaortic septum
    n, o = ao["normal"], G["origin"]
    z = (X - o) @ n
    lo, hi = W_["subaortic_depth_mm"]
    m = septal & keep & (z >= lo) & (z <= hi)
    out["wt_subaortic_septum"] = float(np.median(thick[m])) if m.sum() >= 10 else np.nan
    out["wt_subaortic_septum_p90"] = float(np.percentile(thick[m], 90)) if m.sum() >= 10 else np.nan
    out["sep_bulge_ratio"] = out["wt_subaortic_septum"] / out["wt_mid_septal"]

    # STEP 9.4 — aorto-septal angle in the PLAX plane (aortic axis + mitral centre)
    nA = np.cross(n, mvc - o); nA /= np.linalg.norm(nA)
    lo, hi = W_["aosa_height"]
    m = septal & keep & (np.abs((X - o) @ nA) < W_["plax_strip_mm"]) & (h >= lo) & (h <= hi)
    out["aosa_n_pts"] = int(m.sum())
    if m.sum() >= 10:
        Y = X[m] - X[m].mean(0)
        s = np.linalg.svd(Y, full_matrices=False)[2][0]
        s = s if s @ la > 0 else -s
        out["aosa_deg"] = 180 - float(np.degrees(np.arccos(np.clip(n @ s, -1, 1))))
    else:
        out["aosa_deg"] = np.nan
    out["_wall"] = {"epi": epi, "rvc": rvc, "X": X, "thick": thick, "septal": septal, "keep": keep, "h": h, "sector": sector,
                    "aosa_pts": X[m], "la": la, "e_s": e_s, "e_t": e_t, "mvc": mvc}
    return out
