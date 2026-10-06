"""
QuadForge — автоматическая ретопология в квады для Blender 4.2+ / 5.x.

Собственный солвер (нативное ядро C++17, многопоточное): изотропный ремешинг рабочей сетки,
кросс-поле 4-RoSy, целочисленная параметризация (MIQ), извлечение квадов,
ремонт и релаксация с проекцией на исходник. Резерв — Quadriflow Blender.
"""

bl_info = {
    "name": "QuadForge — Auto Retopology",
    "author": "QuadForge Project",
    "version": (1, 0, 0),
    "blender": (4, 2, 0),
    "location": "3D Viewport > N-панель > QuadForge",
    "description": "Автоматическая ретопология в чистые квады",
    "category": "Mesh",
}

from . import preferences, properties, panels  # noqa: E402
from .operators import op_analyze, op_guides, op_native, op_retopo  # noqa: E402

_MODULES = (preferences, properties, op_native, op_retopo, op_analyze, op_guides, panels)


def register():
    for m in _MODULES:
        m.register()


def unregister():
    for m in reversed(_MODULES):
        m.unregister()
