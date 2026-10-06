"""Быстрый обмен мешами Blender <-> NumPy (foreach_get/foreach_set, без циклов Python)."""

import numpy as np


def vertices(me):
    n = len(me.vertices)
    V = np.empty(n * 3, dtype=np.float32)
    me.vertices.foreach_get("co", V)
    return V.reshape(n, 3)


def triangles(me):
    """Треугольники (m,3) через loop_triangles."""
    me.calc_loop_triangles()
    n = len(me.loop_triangles)
    T = np.empty(n * 3, dtype=np.int32)
    me.loop_triangles.foreach_get("vertices", T)
    return T.reshape(n, 3)


def edges(me):
    n = len(me.edges)
    E = np.empty(n * 2, dtype=np.int32)
    me.edges.foreach_get("vertices", E)
    return E.reshape(n, 2)


def polygons(me):
    """(loop_vertices, poly_start, poly_size)."""
    nl = len(me.loops)
    npoly = len(me.polygons)
    loop_v = np.empty(nl, dtype=np.int32)
    me.loops.foreach_get("vertex_index", loop_v)
    start = np.empty(npoly, dtype=np.int32)
    me.polygons.foreach_get("loop_start", start)
    size = np.empty(npoly, dtype=np.int32)
    me.polygons.foreach_get("loop_total", size)
    return loop_v, start, size


def edge_bool_attribute(me, name):
    """Булев атрибут рёбер (или None)."""
    attr = me.attributes.get(name)
    if attr is None or attr.domain != 'EDGE':
        return None
    n = len(me.edges)
    out = np.zeros(n, dtype=bool)
    try:
        if attr.data_type == 'BOOLEAN':
            attr.data.foreach_get("value", out)
        elif attr.data_type in {'FLOAT', 'INT', 'INT8'}:
            tmp = np.zeros(n, dtype=np.float32 if attr.data_type == 'FLOAT' else np.int32)
            attr.data.foreach_get("value", tmp)
            out = tmp > 0.5
        else:
            return None
    except (AttributeError, RuntimeError, TypeError):
        return None
    return out


def edge_float_attribute(me, name):
    attr = me.attributes.get(name)
    if attr is None or attr.domain != 'EDGE' or attr.data_type != 'FLOAT':
        return None
    out = np.zeros(len(me.edges), dtype=np.float32)
    try:
        attr.data.foreach_get("value", out)
    except (AttributeError, RuntimeError):
        return None
    return out


def edge_seams(me):
    n = len(me.edges)
    out = np.zeros(n, dtype=bool)
    try:
        me.edges.foreach_get("use_seam", out)
        return out
    except (AttributeError, RuntimeError, TypeError):
        pass
    a = edge_bool_attribute(me, "uv_seam")
    return a if a is not None else out


def vertex_group_weights(obj, me, name):
    """Веса группы вершин (0 для вершин вне группы)."""
    vg = obj.vertex_groups.get(name)
    n = len(me.vertices)
    w = np.zeros(n, dtype=np.float32)
    if vg is None:
        return w
    gi = vg.index
    for v in me.vertices:
        for g in v.groups:
            if g.group == gi:
                w[v.index] = g.weight
                break
    return w


def color_attribute_luma(me, name):
    """Яркость цветового атрибута на вершину (0..1); поддерживает POINT и CORNER домены."""
    attr = me.color_attributes.get(name) if hasattr(me, "color_attributes") else None
    if attr is None:
        return None
    n = len(attr.data)
    col = np.zeros(n * 4, dtype=np.float32)
    try:
        attr.data.foreach_get("color", col)
    except (AttributeError, RuntimeError):
        return None
    col = col.reshape(n, 4)
    luma = 0.2126 * col[:, 0] + 0.7152 * col[:, 1] + 0.0722 * col[:, 2]
    if attr.domain == 'POINT':
        return luma.astype(np.float32)
    if attr.domain == 'CORNER':
        loop_v = np.empty(len(me.loops), dtype=np.int32)
        me.loops.foreach_get("vertex_index", loop_v)
        acc = np.zeros(len(me.vertices), dtype=np.float64)
        cnt = np.zeros(len(me.vertices), dtype=np.float64)
        np.add.at(acc, loop_v, luma)
        np.add.at(cnt, loop_v, 1.0)
        return (acc / np.maximum(cnt, 1.0)).astype(np.float32)
    return None


def build_mesh(bpy_data, name, V, face_sizes, indices):
    """Создаёт bpy.types.Mesh из массивов (полигоны произвольной степени)."""
    V = np.ascontiguousarray(V, dtype=np.float32).reshape(-1, 3)
    face_sizes = np.ascontiguousarray(face_sizes, dtype=np.int32)
    indices = np.ascontiguousarray(indices, dtype=np.int32)
    me = bpy_data.meshes.new(name)
    me.vertices.add(V.shape[0])
    me.vertices.foreach_set("co", V.ravel())
    me.loops.add(indices.shape[0])
    me.loops.foreach_set("vertex_index", indices)
    me.polygons.add(face_sizes.shape[0])
    starts = np.zeros(face_sizes.shape[0], dtype=np.int32)
    if face_sizes.shape[0] > 1:
        starts[1:] = np.cumsum(face_sizes[:-1])
    me.polygons.foreach_set("loop_start", starts)
    try:  # Blender < 4.0
        me.polygons.foreach_set("loop_total", face_sizes)
    except (AttributeError, RuntimeError, TypeError):
        pass
    me.update(calc_edges=True)
    me.validate(clean_customdata=False)
    return me


def mirror_arrays(V, face_sizes, indices, axis, eps):
    """Отражение половины сетки по плоскости axis=0 со сшивкой вершин на плоскости."""
    V = np.array(V, dtype=np.float64)
    face_sizes = np.asarray(face_sizes, dtype=np.int64)
    indices = np.asarray(indices, dtype=np.int64)
    # сшиваем только граничные вершины, лежащие на плоскости
    starts0 = np.zeros(face_sizes.shape[0], dtype=np.int64)
    if face_sizes.shape[0] > 1:
        starts0[1:] = np.cumsum(face_sizes[:-1])
    poly_of_loop = np.repeat(np.arange(face_sizes.shape[0]), face_sizes)
    li = np.arange(indices.shape[0])
    nxt = starts0[poly_of_loop] + (li - starts0[poly_of_loop] + 1) % face_sizes[poly_of_loop]
    a = indices
    b = indices[nxt]
    key = np.minimum(a, b) * (V.shape[0] + 1) + np.maximum(a, b)
    uk, cnt = np.unique(key, return_counts=True)
    bkeys = uk[cnt == 1]
    boundary = np.zeros(V.shape[0], dtype=bool)
    boundary[bkeys // (V.shape[0] + 1)] = True
    boundary[bkeys % (V.shape[0] + 1)] = True
    on_plane = boundary & (np.abs(V[:, axis]) <= eps)
    V[on_plane, axis] = 0.0
    n = V.shape[0]
    Vm = V.copy()
    Vm[:, axis] *= -1.0
    # индексы отражённых вершин: вершины на плоскости переиспользуются
    new_index = np.arange(n, dtype=np.int64) + n
    new_index[on_plane] = np.nonzero(on_plane)[0]
    keep_new = ~on_plane
    remap = np.empty(n, dtype=np.int64)
    remap[on_plane] = np.nonzero(on_plane)[0]
    remap[keep_new] = n + np.arange(int(keep_new.sum()))
    V_out = np.vstack([V, Vm[keep_new]])
    # грани отражения — с обратной ориентацией
    face_sizes = np.asarray(face_sizes, dtype=np.int64)
    indices = np.asarray(indices, dtype=np.int64)
    starts = np.zeros(face_sizes.shape[0], dtype=np.int64)
    if face_sizes.shape[0] > 1:
        starts[1:] = np.cumsum(face_sizes[:-1])
    mir = np.empty_like(indices)
    for s, k in zip(starts, face_sizes):
        mir[s:s + k] = remap[indices[s:s + k][::-1]]
    return (V_out.astype(np.float32), np.concatenate([face_sizes, face_sizes]).astype(np.int32),
            np.concatenate([indices, mir]).astype(np.int32))
