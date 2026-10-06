"""Направляющие (guides): пометка рёбер в режиме редактирования + отрисовка во вьюпорте."""

import bmesh
import bpy
import gpu
from bpy.types import Operator
from gpu_extras.batch import batch_for_shader

ATTR = "qf_guide"


def _guide_layer(bm, create=False):
    lay = bm.edges.layers.bool.get(ATTR) if hasattr(bm.edges.layers, "bool") else None
    if lay is None and create:
        lay = bm.edges.layers.bool.new(ATTR)
    return lay


class QUADFORGE_OT_guides_mark(Operator):
    """Пометить выделенные рёбра как направляющие для потока рёбер"""

    bl_idname = "quadforge.guides_mark"
    bl_label = "Пометить направляющие"
    bl_options = {'REGISTER', 'UNDO'}

    clear: bpy.props.BoolProperty(name="Снять", default=False)

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == 'MESH' and obj.mode == 'EDIT'

    def execute(self, context):
        n = 0
        for obj in {o for o in context.objects_in_mode if o.type == 'MESH'}:
            bm = bmesh.from_edit_mesh(obj.data)
            lay = _guide_layer(bm, create=not self.clear)
            if lay is None:
                continue
            for e in bm.edges:
                if e.select:
                    e[lay] = not self.clear
                    n += 1
            bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)
        _invalidate()
        self.report({'INFO'}, f"{'Снято' if self.clear else 'Помечено'} рёбер: {n}")
        return {'FINISHED'}


class QUADFORGE_OT_guides_select(Operator):
    """Выделить рёбра-направляющие"""

    bl_idname = "quadforge.guides_select"
    bl_label = "Выделить направляющие"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == 'MESH' and obj.mode == 'EDIT'

    def execute(self, context):
        obj = context.active_object
        bm = bmesh.from_edit_mesh(obj.data)
        lay = _guide_layer(bm)
        if lay is None:
            return {'CANCELLED'}
        for e in bm.edges:
            if e[lay]:
                e.select = True
        bm.select_flush_mode()
        bmesh.update_edit_mesh(obj.data, loop_triangles=False, destructive=False)
        return {'FINISHED'}


class QUADFORGE_OT_guides_clear_all(Operator):
    """Удалить все направляющие с объекта"""

    bl_idname = "quadforge.guides_clear_all"
    bl_label = "Удалить все направляющие"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == 'MESH'

    def execute(self, context):
        obj = context.active_object
        if obj.mode == 'EDIT':
            bm = bmesh.from_edit_mesh(obj.data)
            lay = _guide_layer(bm)
            if lay is not None:
                bm.edges.layers.bool.remove(lay)
            bmesh.update_edit_mesh(obj.data)
        else:
            attr = obj.data.attributes.get(ATTR)
            if attr is not None:
                obj.data.attributes.remove(attr)
        _invalidate()
        return {'FINISHED'}


# ---------------------------------------------------------------------- отрисовка
_handle = None
_cache = {"key": None, "batch": None}


def _invalidate():
    _cache["key"] = None


def _draw():
    ctx = bpy.context
    s = getattr(ctx.scene, "quadforge", None)
    obj = ctx.active_object
    if s is None or not s.use_guides or obj is None or obj.type != 'MESH':
        return
    me = obj.data
    try:
        if obj.mode == 'EDIT':
            bm = bmesh.from_edit_mesh(me)
            lay = _guide_layer(bm)
            if lay is None:
                return
            mw = obj.matrix_world
            coords = []
            for e in bm.edges:
                if e[lay]:
                    coords.append(mw @ e.verts[0].co)
                    coords.append(mw @ e.verts[1].co)
        else:
            key = (obj.name, me.name, len(me.edges), me.attributes.get(ATTR) is not None)
            if _cache["key"] == key and _cache["batch"] is not None:
                coords = None
            else:
                attr = me.attributes.get(ATTR)
                if attr is None:
                    _cache["key"], _cache["batch"] = key, None
                    return
                mw = obj.matrix_world
                coords = []
                for i, d in enumerate(attr.data):
                    if d.value:
                        e = me.edges[i]
                        coords.append(mw @ me.vertices[e.vertices[0]].co)
                        coords.append(mw @ me.vertices[e.vertices[1]].co)
        shader = gpu.shader.from_builtin('UNIFORM_COLOR')
        if coords is not None:
            if not coords:
                return
            batch = batch_for_shader(shader, 'LINES', {"pos": coords})
            if obj.mode != 'EDIT':
                _cache["key"], _cache["batch"] = key, batch
        else:
            batch = _cache["batch"]
        gpu.state.line_width_set(3.0)
        gpu.state.depth_test_set('LESS_EQUAL')
        shader.uniform_float("color", (1.0, 0.55, 0.05, 1.0))
        batch.draw(shader)
        gpu.state.line_width_set(1.0)
        gpu.state.depth_test_set('NONE')
    except (ReferenceError, RuntimeError, ValueError):
        pass


_CLASSES = (QUADFORGE_OT_guides_mark, QUADFORGE_OT_guides_select, QUADFORGE_OT_guides_clear_all)


def register():
    global _handle
    for c in _CLASSES:
        bpy.utils.register_class(c)
    if not bpy.app.background:
        _handle = bpy.types.SpaceView3D.draw_handler_add(_draw, (), 'WINDOW', 'POST_VIEW')


def unregister():
    global _handle
    if _handle is not None:
        bpy.types.SpaceView3D.draw_handler_remove(_handle, 'WINDOW')
        _handle = None
    for c in reversed(_CLASSES):
        bpy.utils.unregister_class(c)
