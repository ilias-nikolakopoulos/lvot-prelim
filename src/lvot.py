"""Per-mesh LVOT / aortic annulus / mitral ring measurements (UKB representative meshes).

Geometry only (no UVCs): reads the .fiber.vtu once, coordinates um -> mm.
elemTag 0 RV, 1 LV, 4 mitral plug, 5 aortic plug, 6 tricuspid, 7 pulmonary.

Definitions (see README):
- Ring of a plug = LV nodes shared with that plug's tets (a ~2 mm band). Plane by SVD,
  normal oriented into the LV. Diameters from the ring nodes projected on that plane.
- LV cavity = the closed surface component bounded by LV tets and the mitral/aortic plugs
  (outer surface and RV cavity are the other components). Its enclosed volume = LV cavity
  volume (incl. LVOT up to the aortic plug).
- LV endocardium = cavity-surface triangles that come from LV tets (plug faces removed).
- LVOT sections: planes perpendicular to the aortic ring normal at depth d mm below the
  ventricular edge of the aortic ring band (d = 0 is where the plug meets the LV endo on the
  ventricular side). A section counts only if the plane cuts the LV endocardium in a closed loop
  that encircles the axis (closed tube). Depths are kept only while contiguous from d = 0.
"""
from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np

UM_TO_MM = 1e-3
from src.config import load_config

P_ = load_config()["measure"]  # all tunable parameters live in config.yaml
DEPTHS = tuple(float(d) for d in P_["depths_mm"])  # section planes below the ring
REPORT_DEPTHS = tuple(float(d) for d in P_["report_depths_mm"])  # depths written as columns
NAME_RE = re.compile(r"^(female|male)_age_(\d+)_bmi_(\d+)$")
warnings.filterwarnings("ignore", category=FutureWarning)
try:  # pyvista FutureWarning class
    from pyvista.core.errors import PyVistaFutureWarning
    warnings.filterwarnings("ignore", category=PyVistaFutureWarning)
except Exception:  # pragma: no cover
    pass


# ---------- 2D helpers ----------
def fit_ellipse(xy: np.ndarray) -> dict:
    """Direct least-squares ellipse fit (Fitzgibbon, Halir-Flusser form).
    Returns centre, full major/minor diameters (mm) and major-axis angle (rad)."""
    x, y = xy[:, 0] - xy[:, 0].mean(), xy[:, 1] - xy[:, 1].mean()
    s = max(x.std(), y.std())
    x, y = x / s, y / s
    D1 = np.column_stack([x * x, x * y, y * y])
    D2 = np.column_stack([x, y, np.ones_like(x)])
    S1, S2, S3 = D1.T @ D1, D1.T @ D2, D2.T @ D2
    T = -np.linalg.solve(S3, S2.T)
    M = S1 + S2 @ T
    M = np.array([M[2] / 2, -M[1], M[0] / 2])
    w, v = np.linalg.eig(M)
    v = np.real(v)
    cond = 4 * v[0] * v[2] - v[1] ** 2
    a1 = v[:, cond > 0][:, 0]
    A, B, C, Dc, E, F = np.concatenate([a1, T @ a1])
    # centre
    den = B * B - 4 * A * C
    x0 = (2 * C * Dc - B * E) / den
    y0 = (2 * A * E - B * Dc) / den
    Fc = A * x0 * x0 + B * x0 * y0 + C * y0 * y0 + Dc * x0 + E * y0 + F
    lam = np.linalg.eigvalsh(np.array([[A, B / 2], [B / 2, C]]))
    ax = np.sqrt(-Fc / lam)  # semi-axes (lam ascending -> ax descending)
    theta = 0.5 * np.arctan2(B, A - C) + (np.pi / 2 if A > C else 0.0)
    return {
        "cx": x0 * s + xy[:, 0].mean(), "cy": y0 * s + xy[:, 1].mean(),
        "dmax": 2 * ax.max() * s, "dmin": 2 * ax.min() * s, "theta": theta,
    }


def moment_ellipse(xy: np.ndarray) -> dict:
    """Ellipse with the same area centroid and second moments of area as the closed polygon xy.
    Exact for a true ellipse; robust for D-shaped or noisy outlines. Also Feret (caliper) diameters."""
    x, y = xy[:, 0], xy[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    A = cr.sum() / 2
    cx = ((x + x1) * cr).sum() / (6 * A)
    cy = ((y + y1) * cr).sum() / (6 * A)
    Ixx = ((x * x + x * x1 + x1 * x1) * cr).sum() / 12 - A * cx * cx
    Iyy = ((y * y + y * y1 + y1 * y1) * cr).sum() / 12 - A * cy * cy
    Ixy = ((x * y1 + 2 * x * y + 2 * x1 * y1 + x1 * y) * cr).sum() / 24 - A * cx * cy
    lam, vec = np.linalg.eigh(np.array([[Ixx, Ixy], [Ixy, Iyy]]) / A)
    # for an ellipse with semi-axes a, b: I/A = a^2/4, b^2/4
    dmin, dmax = 4 * np.sqrt(lam)
    th = np.arctan2(vec[1, 1], vec[0, 1])
    angs = np.linspace(0, np.pi, 180, endpoint=False)
    w = np.array([np.ptp(xy @ np.array([np.cos(t), np.sin(t)])) for t in angs])
    return {"cx": cx, "cy": cy, "dmax": dmax, "dmin": dmin, "theta": th,
            "feret_max": w.max(), "feret_min": w.min()}


def chord_through(xy: np.ndarray, c: np.ndarray, u: np.ndarray) -> float:
    """Length of the straight line through point c in direction u that lies across the closed polygon xy
    (edge-to-edge distance, as an echo caliper measures it)."""
    p, q = xy, np.roll(xy, -1, axis=0)
    w = q - p
    den = u[0] * w[:, 1] - u[1] * w[:, 0]
    ok = np.abs(den) > 1e-12
    dp = p - c
    t = (dp[:, 0] * w[:, 1] - dp[:, 1] * w[:, 0])[ok] / den[ok]   # position along the line
    sgm = (dp[:, 0] * u[1] - dp[:, 1] * u[0])[ok] / den[ok]      # position along the edge
    t = t[(sgm >= 0) & (sgm < 1)]
    return float(t.max() - t.min()) if len(t) >= 2 else np.nan


def echo_diams(xy: np.ndarray, u: np.ndarray) -> dict:
    """Diameters an echo would measure on this outline. u = unit in-plane direction of the parasternal
    long-axis (PLAX) cut, i.e. towards the mitral valve (septum -> aorto-mitral side).
    d_plax: line through the area centroid along u; d_perp: the orthogonal line; plax_vs_minor_deg:
    angle between u and the minor axis of the moment ellipse (0 = PLAX measures the minor axis)."""
    el = moment_ellipse(xy)
    c = np.array([el["cx"], el["cy"]])
    v = np.array([-u[1], u[0]])
    minor = el["theta"] + np.pi / 2
    ang = abs(np.degrees(np.angle(np.exp(2j * (np.arctan2(u[1], u[0]) - minor))))) / 2  # axis angle 0..90
    return {"d_plax": chord_through(xy, c, u), "d_perp": chord_through(xy, c, v), "plax_vs_minor_deg": ang}


def polygon_area_perim(xy: np.ndarray) -> tuple[float, float]:
    x, y = xy[:, 0], xy[:, 1]
    area = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    perim = np.linalg.norm(xy - np.roll(xy, -1, axis=0), axis=1).sum()
    return area, perim


def resample_loop(xy: np.ndarray, n: int = 180) -> np.ndarray:
    """Uniform arc-length resampling of a closed polyline (xy without repeated end point)."""
    c = np.vstack([xy, xy[:1]])
    seg = np.linalg.norm(np.diff(c, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)])
    t = np.linspace(0, s[-1], n, endpoint=False)
    return np.column_stack([np.interp(t, s, c[:, 0]), np.interp(t, s, c[:, 1])])


def point_in_poly(p, xy) -> bool:
    x, y = xy[:, 0], xy[:, 1]
    inside = False
    j = len(xy) - 1
    for i in range(len(xy)):
        if (y[i] > p[1]) != (y[j] > p[1]) and p[0] < (x[j] - x[i]) * (p[1] - y[i]) / (y[j] - y[i]) + x[i]:
            inside = not inside
        j = i
    return inside


def plane_basis(n):
    a = np.array([1.0, 0, 0]) if abs(n[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(n, a); e1 /= np.linalg.norm(e1)
    return e1, np.cross(n, e1)



def chain_loops(pts, lines, o, e1, e2):
    """Join 2-point line cells from a slice into polylines (merging duplicate points).
    Returns [(closed, xy)] in plane coordinates; closed loops have no repeated end point."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    key, inv = np.unique(np.round(pts, 6), axis=0, return_inverse=True)
    inv = inv.ravel()
    seg = lines.reshape(-1, 3)[:, 1:]
    seg = inv[seg]
    seg = seg[seg[:, 0] != seg[:, 1]]
    seg = np.unique(np.sort(seg, axis=1), axis=0)
    nn = len(key)
    A = coo_matrix((np.ones(len(seg)), (seg[:, 0], seg[:, 1])), shape=(nn, nn)).tocsr()
    A = A + A.T
    deg = np.asarray(A.sum(1)).ravel()
    ncomp, lab = connected_components(A, directed=False)
    nbr = [[] for _ in range(nn)]
    for a, b in seg:
        nbr[a].append(b); nbr[b].append(a)
    out = []
    for c in range(ncomp):
        nodes = np.flatnonzero(lab == c)
        if len(nodes) < 2:
            continue
        closed = bool(np.all(deg[nodes] == 2))
        ends = nodes[deg[nodes] == 1]
        start = nodes[0] if closed or len(ends) == 0 else ends[0]
        order, prev, cur = [start], -1, start
        while True:
            nxt = [q for q in nbr[cur] if q != prev]
            if not nxt or nxt[0] == start:
                break
            prev, cur = cur, nxt[0]
            order.append(cur)
            if len(order) > len(nodes):
                break
        X = key[order] - o
        out.append((closed, np.column_stack([X @ e1, X @ e2])))
    return out

# ---------- rings ----------
def ring_measure(X: np.ndarray, into_lv: np.ndarray, nbins: int = P_["ring_bins"]) -> dict:
    """STEPS 3.2-3.4. X: ring nodes (mm). into_lv: a point inside the LV used to orient the normal."""
    c = X.mean(0)
    _, sv, vt = np.linalg.svd(X - c)
    n = vt[2] if (into_lv - c) @ vt[2] > 0 else -vt[2]
    e1, e2 = plane_basis(n)
    d = (X - c) @ n
    xy = np.column_stack([(X - c) @ e1, (X - c) @ e2])
    el = fit_ellipse(xy)
    cen = np.array([el["cx"], el["cy"]])
    v = xy - cen
    ang = np.arctan2(v[:, 1], v[:, 0])
    r = np.linalg.norm(v, axis=1)
    b = np.floor((ang + np.pi) / (2 * np.pi) * nbins).astype(int).clip(0, nbins - 1)
    rb = np.array([np.median(r[b == k]) if np.any(b == k) else np.nan for k in range(nbins)])
    cover = np.isfinite(rb).mean()
    th = (np.arange(nbins) + 0.5) / nbins * 2 * np.pi - np.pi
    ok = np.isfinite(rb)
    poly = np.column_stack([rb[ok] * np.cos(th[ok]), rb[ok] * np.sin(th[ok])])
    area, perim = polygon_area_perim(poly)
    me = moment_ellipse(poly)
    me["cx"] += cen[0]; me["cy"] += cen[1]
    el = me
    centre3 = c + cen[0] * e1 + cen[1] * e2
    return {
        "centre": centre3, "normal": n, "e1": e1, "e2": e2, "n_nodes": len(X),
        "dmax": me["dmax"], "dmin": me["dmin"], "theta": me["theta"],
        "feret_max": me["feret_max"], "feret_min": me["feret_min"],
        "d_area": 2 * np.sqrt(area / np.pi), "d_perim": perim / np.pi,
        "area": area, "cover": cover, "planarity_rms": float(np.sqrt(np.mean(d ** 2))),
        "band_max_depth": float(d.max()), "xy": xy, "poly": poly + cen, "ell": el,
    }


# ---------- per-mesh ----------
def parse_name(mesh_id: str) -> dict:
    m = NAME_RE.match(mesh_id)
    return {"mesh_id": mesh_id, "sex": m[1][0].upper(), "age": int(m[2]), "bmi": int(m[3])}


def tet_volumes(P, conn):
    a, b, c, d = (P[conn[:, i]] for i in range(4))
    return np.abs(np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a)) / 6


def measure_mesh(folder: str | Path, depths=DEPTHS, keep_geometry: bool = False) -> dict:
    import pyvista as pv

    # STEP 1.1 — read the .vtu once, check tets only, convert um -> mm
    folder = Path(folder)
    mid = folder.name
    g = pv.read(folder / f"{mid}.fiber.vtu")
    tags = np.asarray(g.cell_data["elemTag"])
    conn = g.cells.reshape(-1, 5)
    if not np.all(conn[:, 0] == 4):
        raise ValueError("non-tetrahedral cells")
    conn = conn[:, 1:]
    P = np.asarray(g.points) * UM_TO_MM
    g.points = P
    out = parse_name(mid)
    flags = []

    # STEP 1.2 — LV myocardial volume (sum of tag-1 tet volumes)
    lv_nodes = np.unique(conn[tags == 1])
    vol = tet_volumes(P, conn)
    out["lv_myo_ml"] = vol[tags == 1].sum() / 1000

    # STEP 2.1 — boundary surface of all tets (myocardium + plugs) split into connected components
    surf = g.extract_surface(pass_cellid=True)
    ftag = tags[np.asarray(surf.cell_data["vtkOriginalCellIds"])]
    surf.cell_data["tag"] = ftag
    reg = surf.connectivity("all")
    rid = np.asarray(reg.cell_data["RegionId"])
    rtag = np.asarray(reg.cell_data["tag"])
    # STEP 2.2 — LV cavity = component made of LV + mitral/aortic plug faces only; its volume
    best, best_frac = None, 0
    for r in np.unique(rid):
        sel = rid == r
        frac = np.mean(np.isin(rtag[sel], (1, 4, 5)))
        has_plug = np.any(rtag[sel] == 5) and np.any(rtag[sel] == 4)
        # outer surface also contains LV faces; choose the component made only of LV + mitral/aortic
        if has_plug and frac > best_frac and not np.any(rtag[sel] == 0):
            best, best_frac = r, frac
    if best is None:
        raise ValueError("no closed LV cavity component found")
    cav = reg.extract_cells(rid == best).extract_surface()
    lv_region = int(best)
    out["lv_cavity_ml"] = abs(cav.volume) / 1000
    out["n_regions"] = int(len(np.unique(rid)))
    # STEP 2.3 — LV endocardium = cavity triangles from LV tets (plug faces removed)
    ctag = np.asarray(cav.cell_data["tag"])
    endo = cav.extract_cells(ctag == 1).extract_surface()
    endo_pts = np.asarray(endo.points)

    # STEP 3.1-3.4 — aortic (tag 5) and mitral (tag 4) rings: ring nodes, plane, diameters (ring_measure)
    lvc = endo_pts.mean(0)
    rings = {}
    for t, name in ((5, "ao"), (4, "mv")):
        ids = np.intersect1d(np.unique(conn[tags == t]), lv_nodes)
        rings[name] = ring_measure(P[ids], lvc)
        R = rings[name]
        for k in ("n_nodes", "dmax", "dmin", "feret_max", "feret_min", "d_area", "d_perim", "area", "cover", "planarity_rms"):
            out[f"{name}_{k}"] = R[k]
        out[f"{name}_ellip_ratio"] = R["dmax"] / R["dmin"]
        out[f"{name}_ellip_index"] = 1 - R["dmin"] / R["dmax"]
    # STEP 3.5 — aorto-mitral plane angle and centre distance
    ao, mv = rings["ao"], rings["mv"]
    out["ao_mv_centre_dist"] = np.linalg.norm(ao["centre"] - mv["centre"])
    out["ao_mv_angle"] = np.degrees(np.arccos(np.clip(ao["normal"] @ mv["normal"], -1, 1)))
    # STEP 3.7 — echo-style diameters of the annulus: along the PLAX direction (towards the mitral centre) and orthogonal
    dm = mv["centre"] - ao["centre"]
    u_ring = np.array([dm @ ao["e1"], dm @ ao["e2"]]); u_ring /= np.linalg.norm(u_ring)
    for k, v in echo_diams(ao["poly"], u_ring).items():
        out[f"ao_{k}"] = v

    # STEP 3.6 — LV long axis: mitral ring centre -> farthest LV endo point (apex); LV length; aortic axis angle
    dist = np.linalg.norm(endo_pts - mv["centre"], axis=1)
    apex = endo_pts[np.argmax(dist)]
    out["lv_length"] = dist.max()
    la = (apex - mv["centre"]) / dist.max()
    out["ao_axis_vs_lv_long_angle"] = np.degrees(np.arccos(np.clip(ao["normal"] @ la, -1, 1)))

    # STEP 4.1 — section planes perpendicular to the aortic-ring normal, depth 0 = ventricular edge of the ring band
    n, e1, e2 = ao["normal"], ao["e1"], ao["e2"]
    origin = ao["centre"] + ao["band_max_depth"] * n
    mv_dir = mv["centre"] - origin
    mv_ang = np.arctan2(mv_dir @ e2, mv_dir @ e1)  # in-plane direction of the mitral centre
    u_plax = np.array([np.cos(mv_ang), np.sin(mv_ang)])  # echo PLAX caliper direction in the section plane
    axis_pt = np.zeros(2)  # in-plane axis position = centre of the last strictly closed section
    rmax = P_["cc_radius_factor"] * ao["d_area"]  # chord-closed sections: wall points within this radius of the axis
    sections, cents, run = [], [], "waiting"
    for d in depths:
        rec = {"depth": d, "closed": False}
        o = origin + d * n
        # STEP 4.2 — exact plane x LV endocardium intersection, segments chained into loops
        sl = endo.slice(normal=n, origin=o)
        loops = chain_loops(np.asarray(sl.points), sl.lines, o, e1, e2) if sl.n_points else []
        enc = [xy for c, xy in loops if c and len(xy) > 8 and point_in_poly(axis_pt, xy)]
        # STEP 4.3 — strict: closed loop encircling the axis, not opened into the cavity
        if not loops:
            rec["reason"] = "no intersection"
        elif enc:
            xy = max(enc, key=lambda q: polygon_area_perim(q)[0])  # outermost
            area, perim = polygon_area_perim(xy)
            el = moment_ellipse(xy)
            rec.update(closed=True, area=area, d_area=2 * np.sqrt(area / np.pi), d_perim=perim / np.pi,
                       dmax=el["dmax"], dmin=el["dmin"], feret_max=el["feret_max"], feret_min=el["feret_min"],
                       xy=xy, ell=el, gap_deg=0.0, **echo_diams(xy, u_plax))
            if rec["d_area"] > P_["loop_too_large_factor"] * ao["d_area"]:
                rec.update(closed=False, reason=f"loop too large ({rec['d_area']:.1f} mm): opened into cavity")
            else:
                cents.append(o + xy.mean(0) @ np.vstack([e1, e2]))
                axis_pt = resample_loop(xy).mean(0)
        else:
            rec["reason"] = "open"
        # STEP 4.4 — chord-closed (secondary): find the largest empty angular sector, bridge it with a chord
        if not rec["closed"] and loops:
            # chord-closed: wall points near the axis, ordered by angle; the largest empty sector is the gap
            X = np.vstack([xy for _, xy in loops])
            v = X - axis_pt
            X = X[np.linalg.norm(v, axis=1) < rmax]
            if len(X) > 20:
                v = X - axis_pt
                ang = np.arctan2(v[:, 1], v[:, 0])
                nb = P_["cc_gap_bins"]
                occ = np.zeros(nb, bool)
                occ[np.floor((ang + np.pi) / (2 * np.pi) * nb).astype(int).clip(0, nb - 1)] = True
                if occ.all():
                    gap, gmid = 0.0, np.nan
                else:
                    e = np.concatenate([occ, occ])  # longest run of empty bins (circular)
                    best, cur, end = 0, 0, 0
                    for i, o_ in enumerate(e):
                        cur = 0 if o_ else cur + 1
                        if cur > best:
                            best, end = cur, i
                    best = min(best, nb)
                    gap = best * 360 / nb
                    gmid = ((end - best / 2 + 0.5) % nb) / nb * 2 * np.pi - np.pi
                rec["gap_deg"] = gap
                if np.isfinite(gmid):
                    rec["gap_vs_mitral_deg"] = abs(np.degrees(np.angle(np.exp(1j * (gmid - mv_ang)))))
                if gap <= P_["cc_max_gap_deg"]:
                    rel = np.mod(ang - (gmid if np.isfinite(gmid) else -np.pi), 2 * np.pi)
                    arc = X[np.argsort(rel)]
                    area, perim = polygon_area_perim(arc)
                    el = moment_ellipse(arc)
                    rec.update(cc=True, cc_area=area, cc_d_area=2 * np.sqrt(area / np.pi),
                               cc_dmax=el["dmax"], cc_dmin=el["dmin"], cc_feret_max=el["feret_max"],
                               cc_feret_min=el["feret_min"], cc_xy=arc, cc_ell=el,
                               **{f"cc_{k}": v for k, v in echo_diams(arc, u_plax).items()})
        sections.append(rec)
    # STEP 4.5 — strict run: contiguous closed sections from the first closed depth, which must be 0 or 1 mm
    closed = [s_["closed"] for s_ in sections]
    first = next((i for i, c in enumerate(closed) if c), None)
    ok = first is not None and sections[first]["depth"] <= 1
    for i, s_ in enumerate(sections):
        if ok and i >= first and not closed[i]:
            ok = False
        s_["kept"] = bool(ok and first is not None and i >= first)
    # STEP 4.6 — write strict and chord-closed measures at the report depths
    for s_ in sections:
        d = s_["depth"]
        if d not in REPORT_DEPTHS:
            continue
        dd = f"{d:g}".replace(".", "p")
        for k in ("dmax", "dmin", "feret_max", "feret_min", "d_area", "d_perim", "area",
                  "d_plax", "d_perp", "plax_vs_minor_deg"):
            out[f"lvot{dd}_{k}"] = s_.get(k, np.nan) if s_["kept"] else np.nan
        out[f"lvot{dd}_ellip_ratio"] = out[f"lvot{dd}_dmax"] / out[f"lvot{dd}_dmin"]
        out[f"lvot{dd}_ellip_index"] = 1 - out[f"lvot{dd}_dmin"] / out[f"lvot{dd}_dmax"]
        # chord-closed (secondary); equals the strict value where the section is closed
        src = (s_.get("d_area"), s_.get("dmax"), s_.get("dmin")) if s_["kept"] else \
            (s_.get("cc_d_area"), s_.get("cc_dmax"), s_.get("cc_dmin"))
        out[f"lvot{dd}cc_d_area"], out[f"lvot{dd}cc_dmax"], out[f"lvot{dd}cc_dmin"] = \
            (np.nan if v is None else v for v in src)
        out[f"lvot{dd}cc_ellip_ratio"] = out[f"lvot{dd}cc_dmax"] / out[f"lvot{dd}cc_dmin"]
        for k in ("area", "feret_max", "feret_min", "d_plax", "d_perp", "plax_vs_minor_deg"):
            out[f"lvot{dd}cc_{k}"] = s_.get(k, np.nan) if s_["kept"] else s_.get(f"cc_{k}", np.nan)
        out[f"lvot{dd}cc_gap_deg"] = s_.get("gap_deg", np.nan)
        out[f"lvot{dd}cc_gap_vs_mitral_deg"] = s_.get("gap_vs_mitral_deg", np.nan)
    # STEP 4.7 — closed-tube length and centreline tilt
    kept = [s_["depth"] for s_ in sections if s_["kept"]]
    out["lvot_tube_len"] = max(kept) if kept else np.nan
    out["lvot_n_sections"] = len(kept)
    if len(cents) >= 2:
        C = np.array(cents)
        _, _, vt = np.linalg.svd(C - C.mean(0))
        ax = vt[0] if vt[0] @ n > 0 else -vt[0]
        out["lvot_centreline_vs_ao_normal"] = np.degrees(np.arccos(np.clip(ax @ n, -1, 1)))
    else:
        out["lvot_centreline_vs_ao_normal"] = np.nan

    # STEP 5.1 — QC flags
    if ao["cover"] < P_["ring_cover_min"] or mv["cover"] < P_["ring_cover_min"]:
        flags.append("ring_gap")
    if not (P_["ao_plausible_mm"][0] <= ao["d_area"] <= P_["ao_plausible_mm"][1]):
        flags.append("ao_size_implausible")
    if not (P_["mv_plausible_mm"][0] <= mv["d_area"] <= P_["mv_plausible_mm"][1]):
        flags.append("mv_size_implausible")
    if len(kept) == 0:
        flags.append("no_closed_lvot_section")
    elif max(kept) < 2:
        flags.append("lvot_closed_lt2mm")
    if out["n_regions"] != 3:
        flags.append(f"surface_regions_{out['n_regions']}")
    out["qc_flags"] = ";".join(flags)
    out["section_reasons"] = " | ".join(f"{s['depth']:g}:{s.get('reason','ok')}" + (f"(gap {s['gap_deg']:.0f})" if 'gap_deg' in s and not s['closed'] else '') for s in sections)
    if keep_geometry:
        out["_geom"] = {"rings": rings, "sections": sections, "endo": endo, "apex": apex, "origin": origin,
                        "regions": reg, "lv_region": lv_region, "cavity": cav, "ao_ring_nodes": P[np.intersect1d(np.unique(conn[tags == 5]), lv_nodes)]}
    return out
