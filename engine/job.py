"""
Задача ретопологии одного объекта.

Жизненный цикл:
  prepare()  — основной поток: читаем меш (с модификаторами), симметрию, фичери,
               направляющие, карту плотности; формируем параметры ядра;
  start()    — фоновый поток: нативное ядро (ctypes отпускает GIL, UI не блокируется);
  poll()     — прогресс/этап из общей структуры QFStatus;
  finish()   — основной поток: сборка объекта, зеркалирование, шейдинг, отчёт качества.
Если ядро недоступно или не справилось — резервный путь через Quadriflow Blender.
"""

import math
import threading
import time

import bmesh
import bpy
import numpy as np

from . import meshio, native_bridge, quality

AXES = "XYZ"


class RetopoError(RuntimeError):
    """Ошибка с понятным пользователю сообщением."""


def _mean_scale(obj):
    sc = obj.matrix_world.to_scale()
    return max((abs(sc.x) + abs(sc.y) + abs(sc.z)) / 3.0, 1e-12)


def _density_from_weights(w, mode):
    w = np.clip(np.asarray(w, dtype=np.float32), 0.0, 1.0)
    if mode == 'DECREASE':
        return 1.0 / (1.0 + 3.0 * w)
    return 1.0 + 3.0 * w


def _sym_axes(s):
    return [i for i, flag in enumerate((s.symmetry_x, s.symmetry_y, s.symmetry_z)) if flag]


class RetopoJob:
    def __init__(self, context, obj, settings, prefs=None):
        self.context = context
        self.obj = obj
        self.s = settings
        self.prefs = prefs
        self.status = native_bridge.QFStatus()
        self.thread = None
        self.result = None
        self.error = None
        self.done = False
        self.t0 = time.perf_counter()
        self.elapsed = 0.0
        self.engine_used = ""
        self.report = ""
        self.created = None
        self.native_log = ""
        self._inputs = None
        self._L_local = 0.0
        self._sym = []
        self._fallback_reason = ""

    # ------------------------------------------------------------------ подготовка
    def prepare(self):
        s = self.s
        obj = self.obj
        if obj is None or obj.type != 'MESH':
            raise RetopoError("Объект не является мешем")
        if obj.mode != 'OBJECT':
            raise RetopoError("Переключитесь в объектный режим")
        depsgraph = self.context.evaluated_depsgraph_get()
        owner = None
        tmp_mesh = None
        if s.apply_modifiers:
            owner = obj.evaluated_get(depsgraph)
            me = owner.to_mesh()
        else:
            me = obj.data
        try:
            if len(me.vertices) < 4 or len(me.polygons) == 0:
                raise RetopoError(f"«{obj.name}»: в меше нет геометрии")
            self._sym = _sym_axes(s)
            src = me
            if self._sym:
                bm = bmesh.new()
                bm.from_mesh(me)
                for ax in self._sym:
                    no = [0.0, 0.0, 0.0]
                    no[ax] = 1.0
                    geom = bm.verts[:] + bm.edges[:] + bm.faces[:]
                    bmesh.ops.bisect_plane(bm, geom=geom, dist=1e-7, plane_co=(0.0, 0.0, 0.0),
                                           plane_no=no, clear_inner=True)
                tmp_mesh = bpy.data.meshes.new("qf_tmp_half")
                bm.to_mesh(tmp_mesh)
                bm.free()
                src = tmp_mesh
                if len(src.polygons) == 0:
                    raise RetopoError("После разреза по плоскости симметрии геометрия пуста — "
                                      "проверьте, что центр объекта лежит на плоскости симметрии")
            V = meshio.vertices(src)
            F = meshio.triangles(src)
            E = meshio.edges(src)
            sharp = np.zeros(E.shape[0], dtype=bool)
            if s.use_marked_sharp:
                a = meshio.edge_bool_attribute(src, "sharp_edge")
                if a is not None:
                    sharp |= a
            if s.use_creases:
                c = meshio.edge_float_attribute(src, "crease_edge")
                if c is not None:
                    sharp |= c >= 0.5
            if s.use_seams:
                sharp |= meshio.edge_seams(src)
            guides = None
            if s.use_guides:
                g = meshio.edge_bool_attribute(src, "qf_guide")
                if g is None and src is not obj.data and len(obj.data.edges) == len(src.edges):
                    g = meshio.edge_bool_attribute(obj.data, "qf_guide")
                if g is not None and g.any():
                    guides = E[g]
            dens = None
            if s.density_source == 'VGROUP' and s.density_vgroup:
                w = meshio.vertex_group_weights(obj, src, s.density_vgroup)
                dens = _density_from_weights(w, s.density_mode)
            elif s.density_source == 'COLOR' and s.density_color:
                w = meshio.color_attribute_luma(src, s.density_color)
                if w is not None:
                    dens = _density_from_weights(w, s.density_mode)
            npoly_in = len(me.polygons)
        finally:
            if owner is not None:
                owner.to_mesh_clear()
        if tmp_mesh is not None:
            bpy.data.meshes.remove(tmp_mesh)

        # плотность результата
        sym_div = 2 ** len(self._sym)
        scale = _mean_scale(obj)
        p = None
        core = native_bridge.get_core()
        if core.available:
            p = core.default_params()
        else:
            p = native_bridge.QFParams()
        if s.density_mode_target == 'FACES':
            p.target_faces = max(8, int(s.target_faces / sym_div))
            p.target_edge_length = 0.0
        elif s.density_mode_target == 'EDGE':
            p.target_edge_length = float(s.target_edge_length) / scale
            p.target_faces = 0
        else:
            p.target_faces = max(8, int(npoly_in * s.target_ratio / sym_div))
            p.target_edge_length = 0.0
        area = self._area(V, F)
        if p.target_edge_length > 0:
            self._L_local = p.target_edge_length
        else:
            self._L_local = math.sqrt(area / max(1, p.target_faces))
        p.sharp_angle = float(s.sharp_angle)
        p.detect_sharp = int(bool(s.use_sharp_angle))
        p.preserve_boundary = 1 if (s.preserve_boundary or self._sym) else 0
        p.adaptive = float(s.adaptive)
        p.d_reserved[11] = float(getattr(s, "adaptive_gradation", 0.3))  # ограничение градиента размера
        p.curvature_align = float(s.curvature_align)
        p.guide_strength = float(s.guide_strength)
        p.density_strength = float(s.density_strength)
        p.seed = int(s.seed)
        p.threads = int(getattr(self.prefs, "threads", 0) or 0) if self.prefs else 0
        p.use_gpu = 0
        p.quality = int(s.quality)
        p.relax_iterations = int(s.relax_iterations) if s.relax_iterations > 0 else -1
        p.working_scale = float(s.working_scale)
        p.field_scale = 0.0
        p.verbose = 0
        self._inputs = dict(V=V, F=F, sharp=E[sharp] if sharp.any() else None, guides=guides,
                            dens=dens, params=p, npoly_in=npoly_in)

    @staticmethod
    def _area(V, F):
        a = V[F[:, 1]] - V[F[:, 0]]
        b = V[F[:, 2]] - V[F[:, 0]]
        return float(0.5 * np.linalg.norm(np.cross(a, b), axis=1).sum())

    # ------------------------------------------------------------------ расчёт
    def use_native(self):
        if self.s.engine == 'QUADRIFLOW':
            return False
        core = native_bridge.get_core()
        if not core.available:
            if self.s.engine == 'NATIVE':
                raise RetopoError("Нативное ядро не найдено: " + core.status_text())
            self._fallback_reason = "ядро недоступно"
            return False
        return True

    def start(self):
        if not self.use_native():
            self.done = True
            return
        core = native_bridge.get_core()
        inp = self._inputs

        def worker():
            try:
                self.result = core.run(inp["V"], inp["F"], inp["params"], sharp_edges=inp["sharp"],
                                       guide_edges=inp["guides"], vdensity=inp["dens"], status=self.status)
            except Exception as exc:  # не роняем поток молча
                self.error = exc
            finally:
                self.done = True

        self.thread = threading.Thread(target=worker, name="QuadForgeWorker", daemon=True)
        self.thread.start()

    def run_blocking(self):
        self.start()
        if self.thread is not None:
            self.thread.join()

    def poll(self):
        if self.thread is None:
            return 1.0, "Резервный движок", True
        return float(self.status.progress), self.status.text(), self.done

    def cancel(self):
        self.status.cancel = 1

    # ------------------------------------------------------------------ результат
    def finish(self, context):
        s = self.s
        obj = self.obj
        V = fs = idx = None
        if self.result is not None:
            self.native_log = self.result.log
            if self.result.ok:
                V, fs, idx = self.result.mesh()
                self.engine_used = "QuadForge (нативное ядро)"
            else:
                self._fallback_reason = self.result.error or "ядро не справилось"
            self.result.free()
        elif self.error is not None:
            self._fallback_reason = f"ошибка ядра: {self.error}"
        if self.status.cancel:
            raise RetopoError("Отменено")
        new_obj = None
        if V is not None and len(fs):
            for ax in reversed(self._sym):
                V, fs, idx = meshio.mirror_arrays(V, fs, idx, ax, max(1e-6, 0.05 * self._L_local))
            stats = quality.analyze(V, idx, self._starts(fs), fs)
            ok_quality = stats["quad_ratio"] >= 0.9 and stats["nonmanifold_edges"] == 0
            if not ok_quality and s.engine == 'AUTO' and s.fallback_compare:
                # сравним с резервным движком и возьмём лучшее
                alt = self._run_quadriflow(context)
                if alt is not None:
                    altV, altfs, altidx = alt
                    alt_stats = quality.analyze(altV, altidx, self._starts(altfs), altfs)
                    if alt_stats["score"] > stats["score"]:
                        V, fs, idx, stats = altV, altfs, altidx, alt_stats
                        self.engine_used = "Quadriflow (выбран по качеству)"
            new_obj = self._create_object(context, V, fs, idx)
        else:
            if s.engine == 'NATIVE':
                raise RetopoError(f"Ядро не справилось: {self._fallback_reason}")
            alt = self._run_quadriflow(context)
            if alt is None:
                raise RetopoError(f"Не удалось построить сетку ({self._fallback_reason})")
            V, fs, idx = alt
            stats = quality.analyze(V, idx, self._starts(fs), fs)
            self.engine_used = f"Quadriflow (резерв: {self._fallback_reason})"
            new_obj = self._create_object(context, V, fs, idx)
        self.elapsed = time.perf_counter() - self.t0
        extra = []
        if self.native_log and s.show_log_summary:
            for line in self.native_log.splitlines():
                if line.startswith("итог:") or "сингулярных" in line and "после" in line:
                    extra.append(line[:160])
        self.report = quality.format_report(new_obj.name, stats, self.elapsed, self.engine_used, extra)
        self.created = new_obj
        return new_obj

    @staticmethod
    def _starts(fs):
        fs = np.asarray(fs, dtype=np.int64)
        st = np.zeros(fs.shape[0], dtype=np.int64)
        if fs.shape[0] > 1:
            st[1:] = np.cumsum(fs[:-1])
        return st

    def _create_object(self, context, V, fs, idx):
        s = self.s
        src = self.obj
        name = src.name + s.suffix
        me = meshio.build_mesh(bpy.data, name, V, fs, idx)
        new_obj = bpy.data.objects.new(name, me)
        cols = list(src.users_collection) or [context.collection]
        for col in cols:
            col.objects.link(new_obj)
        new_obj.matrix_world = src.matrix_world.copy()
        if s.shade_smooth:
            me.polygons.foreach_set("use_smooth", np.ones(len(me.polygons), dtype=bool))
            self._mark_sharp_by_angle(me, s.sharp_angle if s.use_sharp_angle else math.radians(60.0))
        else:
            me.polygons.foreach_set("use_smooth", np.zeros(len(me.polygons), dtype=bool))
        me.update()
        if s.add_subdivision:
            mod = new_obj.modifiers.new(name="Subdivision", type='SUBSURF')
            mod.levels = int(s.subdivision_levels)
            mod.render_levels = int(s.subdivision_levels)
        new_obj["quadforge_version"] = "1.0.2"
        new_obj["quadforge_source"] = src.name
        new_obj["quadforge_engine"] = self.engine_used
        if s.hide_source:
            src.hide_set(True)
        for o in context.selected_objects:
            o.select_set(False)
        new_obj.select_set(True)
        context.view_layer.objects.active = new_obj
        return new_obj

    @staticmethod
    def _mark_sharp_by_angle(me, angle):
        """Помечает острые рёбра по углу между нормалями граней (атрибут sharp_edge)."""
        loop_v, start, size = meshio.polygons(me)
        n = len(me.polygons)
        if n == 0:
            return
        N = np.empty(n * 3, dtype=np.float32)
        me.polygons.foreach_get("normal", N)
        N = N.reshape(n, 3)
        E = meshio.edges(me)
        ne = E.shape[0]
        lek = np.empty(len(me.loops), dtype=np.int32)
        me.loops.foreach_get("edge_index", lek)
        poly_of_loop = np.repeat(np.arange(n), size)
        f0 = np.full(ne, -1, dtype=np.int64)
        f1 = np.full(ne, -1, dtype=np.int64)
        # первая и вторая грань ребра
        order = np.argsort(lek, kind="stable")
        le = lek[order]
        pf = poly_of_loop[order]
        first = np.ones(le.shape[0], dtype=bool)
        first[1:] = le[1:] != le[:-1]
        f0[le[first]] = pf[first]
        second = ~first
        f1[le[second]] = pf[second]
        sharp = np.zeros(ne, dtype=bool)
        both = (f0 >= 0) & (f1 >= 0)
        dots = np.einsum("ij,ij->i", N[f0[both]], N[f1[both]])
        sharp[both] = dots < math.cos(angle)
        if not sharp.any():
            return
        attr = me.attributes.get("sharp_edge")
        if attr is None:
            attr = me.attributes.new("sharp_edge", 'BOOLEAN', 'EDGE')
        attr.data.foreach_set("value", sharp)

    # ------------------------------------------------------------------ Quadriflow (резерв)
    def _run_quadriflow(self, context):
        s = self.s
        obj = self.obj
        if not hasattr(bpy.ops.object, "quadriflow_remesh"):
            return None
        depsgraph = context.evaluated_depsgraph_get()
        ev = obj.evaluated_get(depsgraph)
        me = bpy.data.meshes.new_from_object(ev) if s.apply_modifiers else obj.data.copy()
        tmp = bpy.data.objects.new(obj.name + "_qf_tmp", me)
        context.collection.objects.link(tmp)
        tmp.matrix_world = obj.matrix_world.copy()
        try:
            with context.temp_override(object=tmp, active_object=tmp, selected_objects=[tmp],
                                       selected_editable_objects=[tmp]):
                if s.quadriflow_voxel and hasattr(bpy.ops.object, "voxel_remesh"):
                    V = meshio.vertices(me)
                    diag = float(np.linalg.norm(V.max(0) - V.min(0))) or 1.0
                    me.remesh_voxel_size = max(self._L_local / 3.0, diag / 1500.0)
                    me.remesh_voxel_adaptivity = 0.0
                    bpy.ops.object.voxel_remesh()
                kw = dict(mode='FACES', target_faces=max(8, int(self._inputs["params"].target_faces or 1000)),
                          seed=int(s.seed), use_preserve_sharp=bool(s.use_sharp_angle),
                          use_preserve_boundary=bool(s.preserve_boundary),
                          use_mesh_symmetry=False, smooth_normals=False)
                avail = set(bpy.ops.object.quadriflow_remesh.get_rna_type().properties.keys())
                kw = {k: v for k, v in kw.items() if k in avail}
                bpy.ops.object.quadriflow_remesh(**kw)
            qme = tmp.data
            V = meshio.vertices(qme)
            loop_v, start, size = meshio.polygons(qme)
            return V, size, loop_v
        except RuntimeError:
            return None
        finally:
            data = tmp.data
            bpy.data.objects.remove(tmp)
            if data.users == 0:
                bpy.data.meshes.remove(data)
