"""N-панель QuadForge во вьюпорте."""

import bpy
from bpy.types import Panel

from .engine import native_bridge


class _Base:
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = "QuadForge"


class QUADFORGE_PT_main(_Base, Panel):
    bl_label = "QuadForge"

    def draw(self, context):
        layout = self.layout
        s = context.scene.quadforge
        layout.prop(s, "preset")
        col = layout.column(align=True)
        col.prop(s, "density_mode_target", text="")
        if s.density_mode_target == 'FACES':
            col.prop(s, "target_faces")
        elif s.density_mode_target == 'EDGE':
            col.prop(s, "target_edge_length")
        else:
            col.prop(s, "target_ratio")
        col.prop(s, "adaptive", slider=True)
        sub = col.row(align=True)
        sub.enabled = s.adaptive > 0.0
        sub.prop(s, "adaptive_gradation", slider=True)
        col.prop(s, "curvature_align", slider=True)
        row = layout.row(align=True)
        row.label(text="Симметрия:")
        row.prop(s, "symmetry_x", toggle=True)
        row.prop(s, "symmetry_y", toggle=True)
        row.prop(s, "symmetry_z", toggle=True)
        layout.prop(s, "batch_selected")
        row = layout.row(align=True)
        row.scale_y = 1.7
        row.operator("quadforge.remesh", icon='MOD_REMESH')
        report = context.scene.quadforge_report
        if report:
            box = layout.box()
            col = box.column(align=True)
            for line in report.splitlines()[:30]:
                col.label(text=line)
        layout.operator("quadforge.analyze", icon='VIEWZOOM')


class QUADFORGE_PT_features(_Base, Panel):
    bl_label = "Рёбра и фичери"
    bl_parent_id = "QUADFORGE_PT_main"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        s = context.scene.quadforge
        col = self.layout.column()
        row = col.row(align=True)
        row.prop(s, "use_sharp_angle")
        sub = row.row()
        sub.enabled = s.use_sharp_angle
        sub.prop(s, "sharp_angle", text="")
        col.prop(s, "use_marked_sharp")
        col.prop(s, "use_creases")
        col.prop(s, "use_seams")
        col.prop(s, "preserve_boundary")


class QUADFORGE_PT_guides(_Base, Panel):
    bl_label = "Направляющие и плотность"
    bl_parent_id = "QUADFORGE_PT_main"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        s = context.scene.quadforge
        layout = self.layout
        col = layout.column(align=True)
        col.prop(s, "use_guides")
        sub = col.column(align=True)
        sub.enabled = s.use_guides
        sub.prop(s, "guide_strength", slider=True)
        row = sub.row(align=True)
        row.operator("quadforge.guides_mark", text="Пометить").clear = False
        row.operator("quadforge.guides_mark", text="Снять").clear = True
        row = sub.row(align=True)
        row.operator("quadforge.guides_select", text="Выделить")
        row.operator("quadforge.guides_clear_all", text="Удалить все")
        layout.separator()
        layout.prop(s, "density_source")
        obj = context.active_object
        if s.density_source == 'VGROUP' and obj is not None:
            layout.prop_search(s, "density_vgroup", obj, "vertex_groups", text="Группа")
        elif s.density_source == 'COLOR' and obj is not None and obj.type == 'MESH':
            layout.prop_search(s, "density_color", obj.data, "color_attributes", text="Атрибут")
        if s.density_source != 'NONE':
            row = layout.row(align=True)
            row.prop(s, "density_mode", expand=True)
            layout.prop(s, "density_strength", slider=True)


class QUADFORGE_PT_advanced(_Base, Panel):
    bl_label = "Дополнительно"
    bl_parent_id = "QUADFORGE_PT_main"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        s = context.scene.quadforge
        col = self.layout.column()
        col.prop(s, "engine")
        col.prop(s, "quality")
        col.prop(s, "relax_iterations")
        col.prop(s, "working_scale")
        col.prop(s, "seed")
        col.prop(s, "apply_modifiers")
        col.prop(s, "fallback_compare")
        col.prop(s, "quadriflow_voxel")
        col.separator()
        col.prop(s, "suffix")
        col.prop(s, "hide_source")
        col.prop(s, "shade_smooth")
        row = col.row(align=True)
        row.prop(s, "add_subdivision")
        sub = row.row()
        sub.enabled = s.add_subdivision
        sub.prop(s, "subdivision_levels")
        col.prop(s, "show_log_summary")


class QUADFORGE_PT_accel(_Base, Panel):
    bl_label = "Ядро"
    bl_parent_id = "QUADFORGE_PT_main"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        col = self.layout.column(align=True)
        core = native_bridge.get_core()
        col.label(text="Ядро: " + core.status_text(), icon='CHECKMARK' if core.available else 'ERROR')
        col.label(text="Вычисления: многопоточный CPU")
        col.operator("quadforge.reload_core", icon='FILE_REFRESH')


class QUADFORGE_PT_log(_Base, Panel):
    bl_label = "Журнал ядра"
    bl_parent_id = "QUADFORGE_PT_main"
    bl_options = {'DEFAULT_CLOSED'}

    def draw(self, context):
        log = context.scene.quadforge_log
        col = self.layout.column(align=True)
        if not log:
            col.label(text="Пусто")
            return
        for line in log.splitlines()[-60:]:
            col.label(text=line[:150])


_CLASSES = (QUADFORGE_PT_main, QUADFORGE_PT_features, QUADFORGE_PT_guides, QUADFORGE_PT_advanced,
            QUADFORGE_PT_accel, QUADFORGE_PT_log)


def register():
    for c in _CLASSES:
        bpy.utils.register_class(c)


def unregister():
    for c in reversed(_CLASSES):
        bpy.utils.unregister_class(c)
