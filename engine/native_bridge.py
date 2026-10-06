"""
Мост к нативному ядру QuadForge (C++17) через ctypes.

Ядро поставляется собранным (lib/<платформа>/…) и не требует компилятора на
машине пользователя. Python-версия Blender не важна: ABI — чистый C.
"""

import ctypes
import os
import platform
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

API_VERSION = 4


class QFParams(ctypes.Structure):
    _pack_ = 8
    _fields_ = [
        ("struct_size", ctypes.c_int32),
        ("target_faces", ctypes.c_int32),
        ("target_edge_length", ctypes.c_double),
        ("sharp_angle", ctypes.c_double),
        ("detect_sharp", ctypes.c_int32),
        ("preserve_boundary", ctypes.c_int32),
        ("adaptive", ctypes.c_double),
        ("curvature_align", ctypes.c_double),
        ("guide_strength", ctypes.c_double),
        ("density_strength", ctypes.c_double),
        ("seed", ctypes.c_int32),
        ("threads", ctypes.c_int32),
        ("use_gpu", ctypes.c_int32),
        ("quality", ctypes.c_int32),
        ("relax_iterations", ctypes.c_int32),
        ("debug", ctypes.c_int32),
        ("verbose", ctypes.c_int32),
        ("pad0", ctypes.c_int32),
        ("working_scale", ctypes.c_double),
        ("field_scale", ctypes.c_double),
        ("i_reserved", ctypes.c_int32 * 16),
        ("d_reserved", ctypes.c_double * 16),
    ]


class QFInput(ctypes.Structure):
    _pack_ = 8
    _fields_ = [
        ("struct_size", ctypes.c_int32),
        ("nv", ctypes.c_int32),
        ("V", ctypes.c_void_p),
        ("nf", ctypes.c_int32),
        ("n_sharp", ctypes.c_int32),
        ("F", ctypes.c_void_p),
        ("sharp_edges", ctypes.c_void_p),
        ("n_guide", ctypes.c_int32),
        ("pad0", ctypes.c_int32),
        ("guide_edges", ctypes.c_void_p),
        ("vdensity", ctypes.c_void_p),
        ("reserved", ctypes.c_int32 * 8),
    ]


class QFStatus(ctypes.Structure):
    _pack_ = 8
    _fields_ = [
        ("cancel", ctypes.c_int32),
        ("progress", ctypes.c_float),
        ("stage", ctypes.c_int32),
        ("done", ctypes.c_int32),
        ("message", ctypes.c_char * 192),
    ]

    def text(self):
        return self.message.decode("utf-8", "replace")


def _platform_dir():
    sysname = sys.platform
    machine = platform.machine().lower()
    if sysname.startswith("win"):
        return "windows-x64"
    if sysname == "darwin":
        return "macos-arm64" if machine in ("arm64", "aarch64") else "macos-x64"
    return "linux-arm64" if machine in ("arm64", "aarch64") else "linux-x64"


_LIB_NAMES = {
    "windows": ("quadforge_native.dll",),
    "darwin": ("libquadforge_native.dylib",),
    "linux": ("libquadforge_native.so",),
}


def _lib_names():
    if sys.platform.startswith("win"):
        return _LIB_NAMES["windows"]
    if sys.platform == "darwin":
        return _LIB_NAMES["darwin"]
    return _LIB_NAMES["linux"]


def _candidate_paths(preferred=None):
    if preferred:
        yield Path(preferred)
    env = os.environ.get("QUADFORGE_NATIVE_LIB")
    if env:
        yield Path(env)
    root = Path(__file__).resolve().parent.parent
    for name in _lib_names():
        yield root / "lib" / _platform_dir() / name
        yield root / "lib" / name
        yield root / "native" / "build" / name
        yield root / "native" / "build" / "Release" / name


class RunResult:
    """Результат одного запуска ядра."""

    def __init__(self, core, handle):
        self._core = core
        self._h = handle
        self.ok = bool(core.lib.qf_result_ok(handle))
        self.error = core.lib.qf_result_error(handle).decode("utf-8", "replace")
        self.log = core.lib.qf_result_log(handle).decode("utf-8", "replace")
        self._mesh = None

    def mesh(self):
        """(V float32 (n,3), face_sizes int32, indices int32)."""
        if self._mesh is None:
            lib = self._core.lib
            pv = ctypes.POINTER(ctypes.c_float)()
            pfs = ctypes.POINTER(ctypes.c_int32)()
            pidx = ctypes.POINTER(ctypes.c_int32)()
            nv = ctypes.c_int()
            nfc = ctypes.c_int()
            nidx = ctypes.c_int()
            lib.qf_result_mesh(self._h, ctypes.byref(pv), ctypes.byref(nv), ctypes.byref(pfs),
                               ctypes.byref(nfc), ctypes.byref(pidx), ctypes.byref(nidx))
            V = (np.ctypeslib.as_array(pv, shape=(nv.value * 3,)).reshape(-1, 3).copy()
                 if nv.value else np.zeros((0, 3), np.float32))
            fs = np.ctypeslib.as_array(pfs, shape=(nfc.value,)).copy() if nfc.value else np.zeros(0, np.int32)
            idx = np.ctypeslib.as_array(pidx, shape=(nidx.value,)).copy() if nidx.value else np.zeros(0, np.int32)
            self._mesh = (V, fs, idx)
        return self._mesh

    def array(self, name):
        data = ctypes.c_void_p()
        n = ctypes.c_int64()
        dt = ctypes.c_int()
        rc = self._core.lib.qf_result_array(self._h, name.encode(), ctypes.byref(data),
                                            ctypes.byref(n), ctypes.byref(dt))
        if rc != 0 or n.value == 0:
            return None
        if dt.value == 0:
            ptr = ctypes.cast(data, ctypes.POINTER(ctypes.c_double))
            return np.ctypeslib.as_array(ptr, shape=(n.value,)).copy()
        ptr = ctypes.cast(data, ctypes.POINTER(ctypes.c_int32))
        return np.ctypeslib.as_array(ptr, shape=(n.value,)).copy()

    def free(self):
        if self._h:
            self._core.lib.qf_result_free(self._h)
            self._h = None

    def __del__(self):
        try:
            self.free()
        except Exception:
            pass


def _shadow_copy(path):
    """Windows блокирует загруженную DLL: обновить/удалить аддон было бы нельзя.
    Поэтому грузим копию из временной папки (имя зависит от размера и времени файла)."""
    if not sys.platform.startswith("win"):
        return path
    try:
        st = path.stat()
        dst_dir = Path(tempfile.gettempdir()) / "quadforge_core"
        dst_dir.mkdir(parents=True, exist_ok=True)
        dst = dst_dir / f"{path.stem}_{st.st_size}_{st.st_mtime_ns}{path.suffix}"
        if not dst.is_file() or dst.stat().st_size != st.st_size:
            tmp = dst.with_suffix(".tmp")
            shutil.copy2(path, tmp)
            os.replace(tmp, dst)
        return dst
    except OSError:
        return path


class NativeCore:
    def __init__(self):
        self.lib = None
        self.path = None
        self.error = None

    @property
    def available(self):
        return self.lib is not None

    def load(self, preferred=None):
        self.lib = None
        self.path = None
        self.error = None
        last = None
        for cand in _candidate_paths(preferred):
            if not cand.is_file():
                continue
            try:
                lib = ctypes.CDLL(str(_shadow_copy(cand)))
                self._bind(lib)
                ver = lib.qf_api_version()
                if ver != API_VERSION:
                    raise OSError(f"версия ядра {ver} != {API_VERSION}")
                self.lib = lib
                self.path = cand
                return True
            except (OSError, AttributeError) as exc:
                last = exc
        if last is not None:
            self.error = str(last)
        return False

    def _bind(self, lib):
        lib.qf_api_version.restype = ctypes.c_int
        lib.qf_default_params.argtypes = [ctypes.POINTER(QFParams)]
        lib.qf_default_params.restype = ctypes.c_int
        lib.qf_set_threads.argtypes = [ctypes.c_int]
        lib.qf_get_threads.restype = ctypes.c_int
        lib.qf_run.restype = ctypes.c_void_p
        lib.qf_run.argtypes = [ctypes.POINTER(QFParams), ctypes.POINTER(QFInput), ctypes.POINTER(QFStatus)]
        lib.qf_result_ok.argtypes = [ctypes.c_void_p]
        lib.qf_result_ok.restype = ctypes.c_int
        lib.qf_result_error.argtypes = [ctypes.c_void_p]
        lib.qf_result_error.restype = ctypes.c_char_p
        lib.qf_result_log.argtypes = [ctypes.c_void_p]
        lib.qf_result_log.restype = ctypes.c_char_p
        lib.qf_result_mesh.argtypes = [
            ctypes.c_void_p, ctypes.POINTER(ctypes.POINTER(ctypes.c_float)),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.POINTER(ctypes.c_int32)),
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.POINTER(ctypes.c_int32)),
            ctypes.POINTER(ctypes.c_int)]
        lib.qf_result_array.argtypes = [
            ctypes.c_void_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_void_p),
            ctypes.POINTER(ctypes.c_int64), ctypes.POINTER(ctypes.c_int)]
        lib.qf_result_free.argtypes = [ctypes.c_void_p]

    def default_params(self):
        p = QFParams()
        self.lib.qf_default_params(ctypes.byref(p))
        return p

    def run(self, V, F, params, sharp_edges=None, guide_edges=None, vdensity=None, status=None):
        """Синхронный запуск. V (n,3) float32, F (m,3) int32 (треугольники).

        sharp_edges / guide_edges — массивы пар индексов вершин; vdensity — плотность
        на вершину (1.0 = обычная). Вызов отпускает GIL (ctypes), поэтому его можно
        выполнять в фоновом потоке, опрашивая status из основного.
        """
        V = np.ascontiguousarray(V, dtype=np.float32).reshape(-1, 3)
        F = np.ascontiguousarray(F, dtype=np.int32).reshape(-1, 3)
        keep = [V, F]
        inp = QFInput()
        inp.struct_size = ctypes.sizeof(QFInput)
        inp.nv = V.shape[0]
        inp.V = V.ctypes.data
        inp.nf = F.shape[0]
        inp.F = F.ctypes.data
        if sharp_edges is not None and len(sharp_edges):
            se = np.ascontiguousarray(sharp_edges, dtype=np.int32).reshape(-1, 2)
            keep.append(se)
            inp.sharp_edges = se.ctypes.data
            inp.n_sharp = se.shape[0]
        if guide_edges is not None and len(guide_edges):
            ge = np.ascontiguousarray(guide_edges, dtype=np.int32).reshape(-1, 2)
            keep.append(ge)
            inp.guide_edges = ge.ctypes.data
            inp.n_guide = ge.shape[0]
        if vdensity is not None:
            vd = np.ascontiguousarray(vdensity, dtype=np.float32).reshape(-1)
            if vd.shape[0] == V.shape[0]:
                keep.append(vd)
                inp.vdensity = vd.ctypes.data
        st = status if status is not None else QFStatus()
        params.struct_size = ctypes.sizeof(QFParams)
        h = self.lib.qf_run(ctypes.byref(params), ctypes.byref(inp), ctypes.byref(st))
        del keep
        return RunResult(self, h)

    def set_threads(self, n):
        if self.available:
            self.lib.qf_set_threads(int(n))

    def gpu_text(self):
        if not self.available:
            return "—"
        try:
            fn = self.lib.qf_gpu_status
        except AttributeError:
            return "не поддерживается этой сборкой ядра"
        fn.restype = ctypes.c_char_p
        raw = fn()
        return raw.decode("utf-8", "replace") if raw else "нет данных"

    def status_text(self):
        if not self.available:
            if self.error:
                return f"недоступно: {self.error[:80]}"
            return "не найдено"
        return f"загружено (ядро v{self.lib.qf_api_version()}, потоков {self.lib.qf_get_threads()})"


_core = None


def get_core(preferred=None, reload=False):
    global _core
    if _core is None or reload:
        _core = NativeCore()
        _core.load(preferred)
    return _core
