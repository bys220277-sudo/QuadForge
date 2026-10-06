"""Анализ качества активного меша и выделение дефектов."""

from math import radians

import bpy
import numpy as np
from bpy.props import BoolProperty, FloatProperty
from bpy.types import Operator

from ..engine import meshio, quality


class QUADFORGE_OT_analyze(Operator):
    """Проанализировать качество квад-сетки активного объекта"""

    bl_idname = "quadforge.analyze"
    bl_label = "Анализ качества"
    bl_options = {'REGISTER', 'UNDO'}

    select_defects: BoolProperty(name="Выделить дефекты", default=False,
                                 description="Выделить не-квады и полигоны с острыми углами")
    degenerate_angle: FloatProperty(name="Порог угла", default=radians(20.0), min=0.0, max=radians(60.0),
                                    subtype='ANGLE')

    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return obj is not None and obj.type == 'MESH' and obj.mode == 'OBJECT'

    def execute(self, context):
        obj = context.active_object
        me = obj.data
        V = meshio.vertices(me)
        loop_v, start, size = meshio.polygons(me)
        s = quality.analyze(V, loop_v, start, size, degenerate_angle_deg=float(np.degrees(self.degenerate_angle)))
        context.scene.quadforge_report = quality.format_report(obj.name, s)
        if self.select_defects:
            bad = size != 4
            if s["degenerate"]:
                bad[np.asarray(s["degenerate"], dtype=np.int64)] = True
            me.polygons.foreach_set("select", bad)
            me.update()
        self.report({'INFO'}, f"Квады {s['quad_ratio'] * 100:.1f}%, оценка {s['score']:.0f}/100")
        return {'FINISHED'}


def register():
    bpy.utils.register_class(QUADFORGE_OT_analyze)


def unregister():
    bpy.utils.unregister_class(QUADFORGE_OT_analyze)
