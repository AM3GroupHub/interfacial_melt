#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from pathlib import Path

def main():
    try:
        from pymatgen.core import Structure
        from pymatgen.io.vasp import Poscar
    except ImportError:
        print("未找到 pymatgen。请先安装：pip install pymatgen")
        sys.exit(1)

    cwd = Path(".").resolve()
    cif_files = sorted(list(cwd.glob("*.cif")) + list(cwd.glob("*.CIF")))

    if not cif_files:
        print("当前目录未发现 .cif 文件。")
        return

    success, failed = 0, 0
    for cif_path in cif_files:
        stem = cif_path.stem  # 文件名不含扩展名
        out_path = cif_path.with_suffix(".vasp")  # 输出为同名 .vasp（POSCAR 格式）

        try:
            structure = Structure.from_file(str(cif_path))
            # 写为 VASP POSCAR 格式（VASP5 标准包含元素行）
            Poscar(structure).write_file(str(out_path))
            print(f"[OK] {cif_path.name} -> {out_path.name}")
            success += 1
        except Exception as e:
            print(f"[FAIL] {cif_path.name}: {e}")
            failed += 1

    print(f"\n完成：成功 {success} 个，失败 {failed} 个。输出文件为同名 .vasp（POSCAR）格式。")

if __name__ == "__main__":
    main()
