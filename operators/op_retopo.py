"""Оператор ретопологии: фоновый поток для ядра, прогресс, отмена по Esc, батч."""

import bpy
from bpy.types import Operator

from ..engine.job import RetopoError, RetopoJob


def _prefs(context):
    addon = context.preferences.addons.get(__package__.rsplit(".", 1)[0])
    return addon.preferences if addon else None


def _targets(context, s):
    if s.batch_selected:
        return [o for o in context.selected_objects if o.type == 'MESH']
    obj = context.active_object
    return [obj] if obj is not None and obj.type == 'MESH' else []


class QUADFORGE_OT_remesh(Operator):
    """Автоматическая ретопология в квады (QuadForge)"""

    bl_idname = "quadforge.remesh"
    bl_label = "QuadForge: Ретопология"
    bl_options = {'REGISTER', 'UNDO'}

    @classmethod
    def poll(cls, context):
        s = getattr(context.scene, "quadforge", None)
        if s is None or context.mode != 'OBJECT':
            return False
        return bool(_targets(context, s))

    # ---- синхронный режим (скрипты, фон, тесты)
    def execute(self, context):
        s = context.scene.quadforge
        reports, failed = [], 0
        for obj in _targets(context, s):
            job = RetopoJob(context, obj, s, _prefs(context))
            try:
                job.prepare()
                job.run_blocking()
                job.finish(context)
                reports.append(job.report)
                if job.native_log:
                    context.scene.quadforge_log = job.native_log
            except RetopoError as exc:
                failed += 1
                reports.append(f"ОШИБКА ({obj.name}): {exc}")
        context.scene.quadforge_report = "\n\n".join(reports)
        if failed:
            self.report({'WARNING'}, f"QuadForge: ошибок {failed}")
            return {'FINISHED'} if failed < len(reports) else {'CANCELLED'}
        self.report({'INFO'}, "QuadForge: готово")
        return {'FINISHED'}

    # ---- интерактивный режим
    def invoke(self, context, event):
        s = context.scene.quadforge
        self._queue = _targets(context, s)
        if not self._queue:
            self.report({'ERROR'}, "Нет выделенных меш-объектов")
            return {'CANCELLED'}
        self._index = 0
        self._job = None
        self._reports = []
        self._failed = 0
        self._cancelling = False
        wm = context.window_manager
        wm.progress_begin(0.0, 100.0)
        self._timer = wm.event_timer_add(0.1, window=context.window)
        wm.modal_handler_add(self)
        return {'RUNNING_MODAL'}

    def _start_job(self, context):
        s = context.scene.quadforge
        obj = self._queue[self._index]
        self._job = RetopoJob(context, obj, s, _prefs(context))
        self._job.prepare()
        self._job.start()

    def modal(self, context, event):
        if event.type == 'ESC' and event.value == 'PRESS':
            self._cancelling = True
            if self._job is not None:
                self._job.cancel()
            return {'RUNNING_MODAL'}
        if event.type != 'TIMER':
            return {'PASS_THROUGH'}
        try:
            if self._job is None:
                if self._cancelling:
                    return self._end(context, cancelled=True)
                self._start_job(context)
            frac, msg, done = self._job.poll()
            total = 100.0 * (self._index + min(max(frac, 0.0), 1.0)) / len(self._queue)
            context.window_manager.progress_update(total)
            if context.area:
                context.area.header_text_set(
                    f"QuadForge [{self._index + 1}/{len(self._queue)}] {msg or '...'} — {total:.0f}%  (Esc — отмена)")
            if not done:
                return {'RUNNING_MODAL'}
            job = self._job
            self._job = None
            if self._cancelling:
                return self._end(context, cancelled=True)
            job.finish(context)
            self._reports.append(job.report)
            if job.native_log:
                context.scene.quadforge_log = job.native_log
        except RetopoError as exc:
            self._failed += 1
            name = self._queue[self._index].name if self._index < len(self._queue) else "?"
            self._reports.append(f"ОШИБКА ({name}): {exc}")
            self._job = None
            if self._cancelling:
                return self._end(context, cancelled=True)
        except Exception as exc:  # неожиданное: завершаем аккуратно
            self._reports.append(f"Внутренняя ошибка: {exc}")
            self._failed += 1
            self._job = None
            return self._end(context, cancelled=True)
        self._index += 1
        if self._index >= len(self._queue):
            return self._end(context)
        return {'RUNNING_MODAL'}

    def _end(self, context, cancelled=False):
        wm = context.window_manager
        wm.event_timer_remove(self._timer)
        wm.progress_end()
        if context.area:
            context.area.header_text_set(None)
        if self._reports:
            context.scene.quadforge_report = "\n\n".join(self._reports)
        if cancelled:
            self.report({'WARNING'}, "QuadForge: отменено")
            return {'CANCELLED'}
        ok = len(self._queue) - self._failed
        if self._failed:
            self.report({'WARNING'}, f"QuadForge: готово {ok}, ошибок {self._failed}")
        else:
            self.report({'INFO'}, f"QuadForge: готово ({ok})")
        return {'FINISHED'}


def register():
    bpy.utils.register_class(QUADFORGE_OT_remesh)


def unregister():
    bpy.utils.unregister_class(QUADFORGE_OT_remesh)
