"""Перезагрузка и (опционально) сборка нативного ядра из исходников."""

import subprocess
import sys
from pathlib import Path

import bpy
from bpy.types import Operator

from ..engine import native_bridge


class QUADFORGE_OT_reload_core(Operator):
    """Перезагрузить нативное ядро QuadForge"""

    bl_idname = "quadforge.reload_core"
    bl_label = "Перезагрузить ядро"

    def execute(self, context):
        addon = context.preferences.addons.get(__package__.rsplit(".", 1)[0])
        path = addon.preferences.native_lib_path if addon else ""
        core = native_bridge.get_core(preferred=path or None, reload=True)
        self.report({'INFO'} if core.available else {'WARNING'}, "Ядро: " + core.status_text())
        return {'FINISHED'}


class QUADFORGE_OT_build_native(Operator):
    """Собрать ядро из исходников (нужны CMake и компилятор C++17; CUDA не обязательна)"""

    bl_idname = "quadforge.build_native"
    bl_label = "Собрать ядро из исходников"

    def execute(self, context):
        root = Path(__file__).resolve().parent.parent
        script = root / "tools" / "build_native.py"
        if not script.is_file():
            self.report({'ERROR'}, "Не найден tools/build_native.py")
            return {'CANCELLED'}
        try:
            proc = subprocess.run([sys.executable, str(script), "--root", str(root)],
                                  capture_output=True, text=True, timeout=1800)
        except (OSError, subprocess.TimeoutExpired) as exc:
            self.report({'ERROR'}, f"Сборка не запустилась: {exc}")
            return {'CANCELLED'}
        log = (proc.stdout or "") + "\n" + (proc.stderr or "")
        context.scene.quadforge_log = log[-20000:]
        if proc.returncode != 0:
            self.report({'ERROR'}, "Сборка не удалась — см. «Журнал ядра»")
            return {'CANCELLED'}
        core = native_bridge.get_core(reload=True)
        self.report({'INFO'}, "Ядро собрано: " + core.status_text())
        return {'FINISHED'}


_CLASSES = (QUADFORGE_OT_reload_core, QUADFORGE_OT_build_native)


def register():
    for c in _CLASSES:
        bpy.utils.register_class(c)


def unregister():
    for c in reversed(_CLASSES):
        bpy.utils.unregister_class(c)
