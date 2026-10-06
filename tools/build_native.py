"""
Сборка нативного ядра QuadForge под Windows x64.

Способы (по порядку):
  1. CMake + Visual Studio 2019/2022 (cl.exe) — обычный путь на Windows;
  2. zig (pip install ziglang) — кросс-компиляция, работает и с Linux.

Результат: <addon>/lib/windows-x64/quadforge_native.dll

    python build_native.py [--root <папка аддона>] [--method cmake|zig]
"""

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def run(cmd, cwd=None):
    print("$", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], cwd=cwd, check=True)


def build_cmake(src_dir, out_dll):
    cmake = shutil.which("cmake")
    if not cmake:
        raise RuntimeError("CMake не найден")
    build = Path(tempfile.mkdtemp(prefix="quadforge_build_"))
    try:
        run([cmake, "-S", src_dir, "-B", build, "-A", "x64"])
        run([cmake, "--build", build, "--config", "Release", "--parallel"])
        cands = list(build.rglob("quadforge_native.dll"))
        if not cands:
            raise RuntimeError("CMake не создал quadforge_native.dll")
        shutil.copy2(cands[0], out_dll)
    finally:
        shutil.rmtree(build, ignore_errors=True)


def build_zig(src_dir, out_dll):
    try:
        import ziglang  # noqa: F401
    except ImportError as exc:
        raise RuntimeError("zig не установлен (pip install ziglang)") from exc
    sources = sorted((Path(src_dir) / "src").glob("*.cpp"))
    tmp = Path(tempfile.mkdtemp(prefix="quadforge_zig_"))
    try:
        dll = tmp / "quadforge_native.dll"
        run([sys.executable, "-m", "ziglang", "c++", "-target", "x86_64-windows-gnu", "-std=c++17", "-O3",
             "-fno-math-errno", "-DNDEBUG", "-shared", "-o", dll, *sources])
        shutil.copy2(dll, out_dll)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument("--method", choices=("auto", "cmake", "zig"), default="auto")
    a = ap.parse_args()
    root = Path(a.root)
    src_dir = root / "native"
    out_dir = root / "lib" / "windows-x64"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_dll = out_dir / "quadforge_native.dll"
    methods = [a.method] if a.method != "auto" else (["cmake", "zig"] if sys.platform.startswith("win") else ["zig"])
    errors = []
    for m in methods:
        try:
            (build_cmake if m == "cmake" else build_zig)(src_dir, out_dll)
            print("Готово:", out_dll)
            return 0
        except (RuntimeError, subprocess.CalledProcessError, OSError) as exc:
            errors.append(f"{m}: {exc}")
            print("Не удалось:", errors[-1], flush=True)
    print("Сборка не удалась:\n  " + "\n  ".join(errors), file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
