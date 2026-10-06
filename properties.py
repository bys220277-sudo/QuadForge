"""Параметры QuadForge и пресеты."""

from math import radians

import bpy
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    IntProperty,
    PointerProperty,
    StringProperty,
)

PRESETS = {
    'ORGANIC': dict(adaptive=0.5, adaptive_gradation=0.3, curvature_align=0.5, use_sharp_angle=False, sharp_angle=radians(45.0),
                    use_marked_sharp=True, preserve_boundary=True, quality='2', relax_iterations=3),
    'HARD_SURFACE': dict(adaptive=0.25, adaptive_gradation=0.3, curvature_align=0.3, use_sharp_angle=True, sharp_angle=radians(35.0),
                         use_marked_sharp=True, preserve_boundary=True, quality='2', relax_iterations=2),
    'CAD_LIKE': dict(adaptive=0.0, adaptive_gradation=0.3, curvature_align=0.2, use_sharp_angle=True, sharp_angle=radians(25.0),
                     use_marked_sharp=True, preserve_boundary=True, quality='3', relax_iterations=1),
    'GAME_ASSET': dict(adaptive=0.0, adaptive_gradation=0.3, curvature_align=0.5, use_sharp_angle=True, sharp_angle=radians(50.0),
                       use_marked_sharp=True, preserve_boundary=True, quality='1', relax_iterations=4),
}


def _apply_preset(self, _context):
    vals = PRESETS.get(self.preset)
    if not vals:
        return
    for k, v in vals.items():
        setattr(self, k, v)


def _vgroup_items(self, context):
    obj = context.active_object if context else None
    items = [('', "—", "")]
    if obj is not None and obj.type == 'MESH':
        items += [(vg.name, vg.name, "") for vg in obj.vertex_groups]
    return items


class QuadForgeSettings(bpy.types.PropertyGroup):
    """Все параметры ретопологии."""

    preset: EnumProperty(
        name="Пресет",
        description="Готовые наборы параметров под типовые задачи",
        items=[
            ('ORGANIC', "Органика", "Персонажи, скульпты, органические формы"),
            ('HARD_SURFACE', "Хард-серфейс", "Техника, оружие, продукт-дизайн: сохранение острых рёбер"),
            ('CAD_LIKE', "CAD", "Механические детали: точное следование фичерям, без адаптивности"),
            ('GAME_ASSET', "Игровой ассет", "Равномерная сетка под текстурирование и анимацию"),
            ('CUSTOM', "Свой", "Параметры задаются вручную"),
        ],
        default='ORGANIC',
        update=_apply_preset,
    )

    # --- плотность
    density_mode_target: EnumProperty(
        name="Плотность",
        items=[
            ('FACES', "Число квадов", "Целевое число четырёхугольников"),
            ('EDGE', "Длина ребра", "Целевая длина ребра в единицах сцены"),
            ('RATIO', "Доля", "Доля от числа полигонов исходника"),
        ],
        default='FACES',
    )
    target_faces: IntProperty(name="Квадов", description="Приблизительное число квадов в результате",
                              default=5000, min=16, max=5_000_000, soft_max=200_000)
    target_edge_length: FloatProperty(name="Длина ребра", default=0.05, min=1e-6, precision=4,
                                      subtype='DISTANCE', unit='LENGTH')
    target_ratio: FloatProperty(name="Доля", default=0.25, min=0.0001, max=1.0, subtype='FACTOR')
    adaptive: FloatProperty(
        name="Плотность по деталям",
        description="0 — равномерные квады. Больше — мельче там, где поверхность сильнее меняется "
                    "(складки, черты лица, мелкие формы), и крупнее на гладких участках (затылок, "
                    "спина). Общее число квадов сохраняется",
        default=0.5, min=0.0, max=1.0, subtype='FACTOR')
    adaptive_gradation: FloatProperty(
        name="Плавность перехода",
        description="Насколько быстро размер квадов может меняться по поверхности: меньше — "
                    "плавнее переход от плотных зон к разреженным",
        default=0.3, min=0.05, max=1.0, subtype='FACTOR')

    # --- поток рёбер
    curvature_align: FloatProperty(
        name="Следование кривизне",
        description="Выравнивать рёбра по главным направлениям кривизны там, где это надёжно "
                    "(трубы, цилиндры, тела вращения)",
        default=0.5, min=0.0, max=1.0, subtype='FACTOR')
    use_guides: BoolProperty(
        name="Направляющие",
        description="Учитывать рёбра, помеченные как направляющие (QuadForge ▸ Направляющие)",
        default=True)
    guide_strength: FloatProperty(name="Сила направляющих", default=0.7, min=0.0, max=1.0, subtype='FACTOR')

    # --- фичери
    use_sharp_angle: BoolProperty(name="Острые рёбра по углу", default=False,
                                  description="Автоматически находить острые рёбра по углу между гранями")
    sharp_angle: FloatProperty(name="Угол", default=radians(35.0), min=radians(1.0), max=radians(179.0),
                               subtype='ANGLE')
    use_marked_sharp: BoolProperty(name="Помеченные острые", default=True,
                                   description="Рёбра, помеченные как Sharp, становятся линиями сетки")
    use_creases: BoolProperty(name="Складки (Crease ≥ 0.5)", default=True)
    use_seams: BoolProperty(name="UV-швы", default=False)
    preserve_boundary: BoolProperty(name="Сохранять границы", default=True,
                                    description="Открытые границы меша остаются рёбрами результата")

    # --- карта плотности
    density_source: EnumProperty(
        name="Карта плотности",
        items=[('NONE', "Нет", ""), ('VGROUP', "Группа вершин", ""), ('COLOR', "Цвет вершин", "")],
        default='NONE')
    density_vgroup: StringProperty(name="Группа вершин", default="")
    density_color: StringProperty(name="Цветовой атрибут", default="")
    density_mode: EnumProperty(
        name="Режим",
        items=[('INCREASE', "Сгущать", "Окрашенные места получают более мелкие квады"),
               ('DECREASE', "Разрежать", "Окрашенные места получают более крупные квады")],
        default='INCREASE')
    density_strength: FloatProperty(name="Сила", default=0.7, min=0.0, max=1.0, subtype='FACTOR')

    # --- симметрия
    symmetry_x: BoolProperty(name="X", default=False, description="Симметрия по локальной оси X объекта")
    symmetry_y: BoolProperty(name="Y", default=False)
    symmetry_z: BoolProperty(name="Z", default=False)

    # --- движок и качество
    engine: EnumProperty(
        name="Движок",
        items=[('AUTO', "Авто", "Нативное ядро QuadForge; при сбое — Quadriflow"),
               ('NATIVE', "QuadForge", "Только собственное ядро"),
               ('QUADRIFLOW', "Quadriflow", "Встроенный Quadriflow Blender")],
        default='AUTO')
    quality: EnumProperty(
        name="Качество",
        items=[('0', "Быстро", "Минимум попыток"), ('1', "Норма", ""),
               ('2', "Высокое", "Повторные попытки при сложной геометрии"), ('3', "Максимум", "")],
        default='2')
    relax_iterations: IntProperty(name="Релаксация", description="Итерации релаксации результата "
                                  "с проекцией на исходник (0 — выключено)", default=3, min=0, max=30)
    working_scale: FloatProperty(
        name="Детализация расчёта",
        description="Во сколько раз рабочая треугольная сетка мельче квадов (0 — авто). "
                    "Больше — надёжнее на сложных фичерях, но медленнее",
        default=0.0, min=0.0, max=8.0)
    seed: IntProperty(name="Seed", default=1, min=0, max=999999)
    apply_modifiers: BoolProperty(name="С модификаторами", default=True,
                                  description="Брать геометрию с применёнными модификаторами")
    fallback_compare: BoolProperty(name="Сравнивать с Quadriflow", default=True,
                                   description="Если качество ниже порога — построить и Quadriflow-вариант "
                                               "и оставить лучший")
    quadriflow_voxel: BoolProperty(name="Воксельная подготовка для Quadriflow", default=True)

    # --- вывод
    suffix: StringProperty(name="Суффикс", default="_QF")
    hide_source: BoolProperty(name="Скрыть исходник", default=False)
    shade_smooth: BoolProperty(name="Гладкий шейдинг", default=True)
    add_subdivision: BoolProperty(name="Subdivision", default=False)
    subdivision_levels: IntProperty(name="Уровни", default=1, min=1, max=4)
    batch_selected: BoolProperty(name="Все выделенные", default=False,
                                 description="Обработать все выделенные меши по очереди")
    show_log_summary: BoolProperty(name="Сводка ядра в отчёте", default=True)


def register():
    bpy.utils.register_class(QuadForgeSettings)
    bpy.types.Scene.quadforge = PointerProperty(type=QuadForgeSettings)
    bpy.types.Scene.quadforge_report = StringProperty(name="Отчёт QuadForge", default="")
    bpy.types.Scene.quadforge_log = StringProperty(name="Лог ядра QuadForge", default="")


def unregister():
    del bpy.types.Scene.quadforge_log
    del bpy.types.Scene.quadforge_report
    del bpy.types.Scene.quadforge
    bpy.utils.unregister_class(QuadForgeSettings)
