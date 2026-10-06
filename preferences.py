"""Настройки аддона: ядро, потоки CPU."""

import bpy
from bpy.props import IntProperty, StringProperty
from bpy.types import AddonPreferences

from .engine import native_bridge


def _reload_core(self, _context):
    native_bridge.get_core(preferred=self.native_lib_path or None, reload=True)


class QuadForgePreferences(AddonPreferences):
    bl_idname = __package__

    native_lib_path: StringProperty(
        name="Своё ядро",
        description="Путь к собственной сборке ядра (.dll/.so/.dylib). Пусто — встроенное",
        subtype='FILE_PATH', default="", update=_reload_core)
    threads: IntProperty(name="Потоки CPU", description="0 — все ядра процессора", default=0, min=0, max=256)

    def draw(self, context):
        layout = self.layout
        core = native_bridge.get_core(preferred=self.native_lib_path or None)
        box = layout.box()
        box.label(text="Ядро: " + core.status_text(), icon='CHECKMARK' if core.available else 'ERROR')
        if core.path:
            box.label(text=str(core.path))
        box.prop(self, "threads")
        box.prop(self, "native_lib_path")
        row = box.row()
        row.operator("quadforge.reload_core", icon='FILE_REFRESH')
        row.operator("quadforge.build_native", icon='TOOL_SETTINGS')


def register():
    bpy.utils.register_class(QuadForgePreferences)


def unregister():
    bpy.utils.unregister_class(QuadForgePreferences)
