"""
Метрики качества полигональной (квад) сетки — только NumPy.

analyze(V, loop_v, poly_start, poly_size) -> dict:
  faces, quads, tris, ngons, quad_ratio
  valence_hist {'2','3','4','5','6+'} по внутренним вершинам, poles, poles_ratio
  angle_min/mean/max, angle_bad_ratio (вне 45..135°)
  edge_mean/min/max/std, edge_cv
  corner_sin_min/mean (1 = прямые углы)
  boundary_edges, nonmanifold_edges
  degenerate — индексы полигонов с углом < порога
  score — сводная оценка 0..100 (больше — лучше)
"""

import numpy as np


def _corner_angles(V, loop_v, poly_start, poly_size):
    nl = loop_v.shape[0]
    poly_of_loop = np.repeat(np.arange(poly_size.shape[0]), poly_size)
    idx = np.arange(nl)
    start = poly_start[poly_of_loop]
    size = poly_size[poly_of_loop]
    nxt = start + (idx - start + 1) % size
    prv = start + (idx - start - 1) % size
    p = V[loop_v]
    a = V[loop_v[prv]] - p
    b = V[loop_v[nxt]] - p
    cr = np.linalg.norm(np.cross(a, b), axis=1)
    dt = np.einsum("ij,ij->i", a, b)
    ang = np.degrees(np.arctan2(cr, dt))
    return ang, poly_of_loop, nxt


def analyze(V, loop_v, poly_start, poly_size, degenerate_angle_deg=20.0):
    V = np.asarray(V, dtype=np.float64).reshape(-1, 3)
    loop_v = np.asarray(loop_v, dtype=np.int64)
    poly_start = np.asarray(poly_start, dtype=np.int64)
    poly_size = np.asarray(poly_size, dtype=np.int64)
    r = {}
    nfaces = int(poly_size.shape[0])
    r["faces"] = nfaces
    r["quads"] = int(np.count_nonzero(poly_size == 4))
    r["tris"] = int(np.count_nonzero(poly_size == 3))
    r["ngons"] = int(np.count_nonzero(poly_size > 4))
    r["quad_ratio"] = r["quads"] / nfaces if nfaces else 0.0
    if nfaces == 0:
        r.update(valence_hist={}, poles=0, poles_ratio=0.0, angle_min=0.0, angle_mean=0.0, angle_max=0.0,
                 angle_bad_ratio=0.0, edge_mean=0.0, edge_min=0.0, edge_max=0.0, edge_std=0.0, edge_cv=0.0,
                 corner_sin_min=0.0, corner_sin_mean=0.0, boundary_edges=0, nonmanifold_edges=0,
                 degenerate=[], score=0.0)
        return r
    ang, poly_of_loop, nxt = _corner_angles(V, loop_v, poly_start, poly_size)
    # рёбра
    a = loop_v
    b = loop_v[nxt]
    lo = np.minimum(a, b)
    hi = np.maximum(a, b)
    key = lo * (V.shape[0] + 1) + hi
    uk, cnt = np.unique(key, return_counts=True)
    E = np.stack([uk // (V.shape[0] + 1), uk % (V.shape[0] + 1)], axis=1)
    r["boundary_edges"] = int(np.count_nonzero(cnt == 1))
    r["nonmanifold_edges"] = int(np.count_nonzero(cnt > 2))
    deg = np.bincount(E.ravel(), minlength=V.shape[0])
    boundary_v = np.zeros(V.shape[0], dtype=bool)
    if (cnt == 1).any():
        boundary_v[E[cnt == 1].ravel()] = True
    used = deg > 0
    interior = used & ~boundary_v
    d = deg[interior]
    hist = {"2": int(np.count_nonzero(d <= 2)), "3": int(np.count_nonzero(d == 3)),
            "4": int(np.count_nonzero(d == 4)), "5": int(np.count_nonzero(d == 5)),
            "6+": int(np.count_nonzero(d >= 6))}
    r["valence_hist"] = hist
    r["poles"] = int(d.shape[0] - hist["4"])
    r["poles_ratio"] = r["poles"] / max(int(d.shape[0]), 1)
    r["angle_min"] = float(ang.min())
    r["angle_max"] = float(ang.max())
    r["angle_mean"] = float(ang.mean())
    qmask_loop = poly_size[poly_of_loop] == 4
    qa = ang[qmask_loop]
    r["angle_bad_ratio"] = float(np.count_nonzero((qa < 45.0) | (qa > 135.0)) / max(qa.shape[0], 1))
    sins = np.sin(np.radians(qa)) if qa.size else np.zeros(1)
    r["corner_sin_min"] = float(sins.min())
    r["corner_sin_mean"] = float(sins.mean())
    el = np.linalg.norm(V[E[:, 0]] - V[E[:, 1]], axis=1)
    r["edge_mean"] = float(el.mean())
    r["edge_min"] = float(el.min())
    r["edge_max"] = float(el.max())
    r["edge_std"] = float(el.std())
    r["edge_cv"] = r["edge_std"] / max(r["edge_mean"], 1e-30)
    pmin = np.full(nfaces, 360.0)
    np.minimum.at(pmin, poly_of_loop, ang)
    r["degenerate"] = np.nonzero(pmin < degenerate_angle_deg)[0].tolist()
    # сводная оценка
    score = 100.0
    score -= 60.0 * (1.0 - r["quad_ratio"])
    score -= 40.0 * r["poles_ratio"]
    score -= 30.0 * r["angle_bad_ratio"]
    score -= 25.0 * (1.0 - r["corner_sin_mean"])
    score -= 10.0 * min(r["edge_cv"], 1.0)
    score -= 20.0 * min(1.0, r["nonmanifold_edges"] / max(nfaces, 1) * 100)
    score -= 20.0 * min(1.0, len(r["degenerate"]) / max(nfaces, 1) * 50)
    r["score"] = float(max(0.0, score))
    return r


def format_report(name, s, elapsed=None, engine=None, extra=None):
    L = [f"QuadForge · {name}"]
    L.append(f"Полигоны: {s['faces']} (квады {s['quad_ratio'] * 100:.1f}%, треуг. {s['tris']}, n-угольн. {s['ngons']})")
    h = s.get("valence_hist", {})
    if h:
        L.append(f"Валентности внутр.: 3→{h['3']} 4→{h['4']} 5→{h['5']} 6+→{h['6+']} (полюсов {s['poles_ratio'] * 100:.1f}%)")
    L.append(f"Углы: min {s['angle_min']:.1f}° · ср. {s['angle_mean']:.1f}° · max {s['angle_max']:.1f}° · "
             f"вне 45–135°: {s['angle_bad_ratio'] * 100:.1f}%")
    L.append(f"Рёбра: ср. {s['edge_mean']:.4g} · разброс (CV) {s['edge_cv']:.2f} · ортогональность {s['corner_sin_mean']:.3f}")
    if s.get("nonmanifold_edges"):
        L.append(f"Неманифольдных рёбер: {s['nonmanifold_edges']}")
    if s.get("degenerate"):
        L.append(f"Вырожденных полигонов (<20°): {len(s['degenerate'])}")
    L.append(f"Оценка качества: {s['score']:.0f}/100")
    tail = []
    if elapsed is not None:
        tail.append(f"время {elapsed:.2f} с")
    if engine:
        tail.append(f"движок: {engine}")
    if tail:
        L.append(" · ".join(tail))
    if extra:
        L.extend(extra)
    return "\n".join(L)
