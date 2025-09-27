#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
import shutil
import bisect
from pymatgen.io.vasp.outputs import Vasprun


def main():
    # 固定的 work 根目录（按需修改）
    work_root = "/data/home/ycshen/openMX_wannier/multi_test/work"
    wtin_dir = "/data/home/ycshen/openMX_wannier/multi_test/src/wt.in"

    ap = argparse.ArgumentParser(
        description="先用 BandStructure.is_metal() 判 EF 是否被能带穿越；有 DOSCAR 时同时检查 DOS(EF)；为金属则收集 WR/wannier90_hr.dat 到 work/hrs/"
    )
    ap.add_argument("--no-rename", action="store_true",
                    help="复制到 hrs/ 时不重命名（可能覆盖）")
    ap.add_argument("--dos-threshold", type=float, default=1,
                    help="DOS(EF) 判金属的阈值（默认 1e-3 states/eV）")
    args = ap.parse_args()

    hrs_dir = os.path.join(work_root, "hrs")
    if not os.path.isdir(work_root):
        raise SystemExit(f"工作目录不存在：{work_root}")

    subdirs = [os.path.join(work_root, d) for d in os.listdir(work_root)
               if os.path.isdir(os.path.join(work_root, d)) and d != "hrs"]

    total = len(subdirs)
    copied = 0
    skipped = 0
    failed  = 0

    print(f"[开始] 扫描 {total} 个材料目录；输出到：{hrs_dir}")
    for mat_dir in sorted(subdirs):
        name = os.path.basename(mat_dir)

        # ---- 1) BandStructure 判据（核心） ----
        # vrp = find_vasprun(mat_dir)
        # if vrp is None:
        #     band_flag, band_detail = None, "未找到 vasprun.xml(.gz)"
        # else:
        #     band_flag, band_detail = band_is_metal(vrp)
        band_flag = True
        # ---- 2) DOSCAR 判据（辅助/提示） ----
        dos_flag, dos_detail = doscar_judgement(mat_dir, args.dos_threshold)

        # ---- 3) 最终判定：二者之一为 True 即判为金属 ----
        if dos_flag is True:
            # 矛盾提示：Band 金属但 DOS≈0，或 DOS 金属但 Band 非金属
            if band_flag is True and dos_flag is False:
                print(f"[警告] {dos_detail}显示为非金属")
            if band_flag is False and dos_flag is True:
                print(f"[警告] {name}: DOS 显示金属（{dos_detail}），但 Band 非金属（{band_detail}）。检查路径：是否用的是Γ线而非网格？")

            ok, info = copy_hr(mat_dir, hrs_dir, rename=(not args.no_rename))
            ok, info = write_wt_in_from_poscar(mat_dir, wtin_dir, hrs_dir)
            if ok:
                src_hint = "Band" if band_flag is True else "DOS"
                print(f"[复制] {name}: 判定金属（依据：{src_hint}）。{(' '+dos_detail) if isinstance(dos_detail,str) else ''} 已复制到 {info}")
                copied += 1
            else:
                print(f"[缺失] {name}: 判定金属，但 {info}")
                failed += 1

        elif (dos_flag is False):
            # 二者都认定为绝缘/半导体
            print(f"[绝缘] {name}: 非金属。")
            skipped += 1

        else:
            # 无法判定（至少一个为 None，另一个也不是 True）
            details = []
            if band_flag is None: details.append(f"Band 不可用（{band_detail}）")
            else: details.append(f"Band 判定={band_flag}")
            if dos_flag is None: details.append(f"DOS 不可用（{dos_detail}）")
            else: details.append(f"DOS 判定={dos_flag}（{dos_detail}）")
            print(f"[跳过] {name}: 无法可靠判定：{'; '.join(details)}")
            failed += 1

    print("\n[完成] 汇总：")
    print(f"  材料总数      : {total}")
    print(f"  已复制（金属）: {copied}")
    print(f"  已跳过（绝缘）: {skipped}")
    print(f"  失败/缺文件   : {failed}")





# ===================== 基本工具 =====================


def _replace_section(text, header, new_block):
    """用 new_block 替换以 header 开头的区块；没有则在末尾追加。"""
    title_re = re.compile(r"^[A-Z_ ]+\s*$", re.M)
    m = re.search(rf"(?m)^{re.escape(header)}\s*$", text)
    if not m:
        if not text.endswith("\n"):
            text += "\n"
        if not new_block.endswith("\n"):
            new_block += "\n"
        return text + "\n" + new_block

    start = m.start()
    m2 = title_re.search(text, pos=m.end()+1)
    end = m2.start() if m2 else len(text)

    before, after = text[:start], text[end:]
    if not before.endswith("\n"):
        before += "\n"
    if not new_block.endswith("\n"):
        new_block += "\n"
    return before + new_block + after

def _ensure_extra_blank_line_between_blocks(text):
    """确保关键区块之间至少有一行空白。"""
    # LATTICE 和 ATOM_POSITIONS
    text = re.sub(r'(?ms)(^LATTICE.*?\n)(^ATOM_POSITIONS\b)', r'\1\n\2', text)
    # ATOM_POSITIONS 和 PROJECTORS
    text = re.sub(r'(?ms)(^ATOM_POSITIONS.*?\n)(^PROJECTORS\b)', r'\1\n\2', text)
    # PROJECTORS 和 KCUBE_BULK（若模板含该块）
    text = re.sub(r'(?ms)(^PROJECTORS.*?\n)(^KCUBE_BULK\b)', r'\1\n\2', text)
    return text

def _replace_key_value_lines(text, hr_filename, efermi):
    """替换 Hrfile 与 E_FERMI（保留原行注释）。"""
    text = re.sub(r"(?m)^\s*Hrfile\s*=\s*'.*?'\s*$",
                  f"Hrfile = '{hr_filename}'", text)
    def _repl_ef(m):
        parts = m.group(0).split("!")
        new = f"E_FERMI = {efermi:.6f}"
        if len(parts) > 1:
            new += "    !" + "!".join(parts[1:])
        return new
    text = re.sub(r"(?m)^\s*E_FERMI\s*=\s*.*$", _repl_ef, text)
    return text

# ===== POSCAR 解析 & 格式化 =====

def _parse_poscar(poscar_path):
    """返回: (expanded_symbols, counts, coord_type, lattice(3x3), coords(Nx3))"""
    with open(poscar_path, "r") as f:
        lines = [l.rstrip() for l in f if l.strip() != "" or True]

    if len(lines) < 8:
        raise ValueError(f"POSCAR太短: {poscar_path}")

    scale = float(lines[1].split()[0])

    a = [float(x) for x in lines[2].split()]
    b = [float(x) for x in lines[3].split()]
    c = [float(x) for x in lines[4].split()]
    lattice = [[e*scale for e in a],
               [e*scale for e in b],
               [e*scale for e in c]]

    def _is_int_list(s):
        try:
            _ = [int(x) for x in s.split()]
            return True
        except Exception:
            return False

    if _is_int_list(lines[5]):     # 老格式
        species = []
        counts = [int(x) for x in lines[5].split()]
        start_idx = 6
    else:                          # VASP5
        species = lines[5].split()
        counts = [int(x) for x in lines[6].split()]
        start_idx = 7

    # Selective dynamics
    coord_type_line = lines[start_idx].strip().lower()
    if coord_type_line.startswith('s'):
        start_idx += 1
        coord_type_line = lines[start_idx].strip().lower()

    coord_type = 'Direct' if coord_type_line.startswith('d') else 'Cartesian'
    start_idx += 1
    natoms = sum(counts)

    coords = []
    for i in range(natoms):
        toks = lines[start_idx + i].split()
        coords.append([float(toks[0]), float(toks[1]), float(toks[2])])

    if not species:
        species = [f"X{i+1}" for i in range(len(counts))]

    expanded_symbols = []
    for sym, n in zip(species, counts):
        expanded_symbols.extend([sym]*n)

    return expanded_symbols, counts, coord_type, lattice, coords

def _format_lattice_block(lattice):
    lines = ["LATTICE", "Angstrom"]
    for vec in lattice:
        lines.append("   " + "    ".join(f"{v:.16f}" for v in vec))
    lines.append("")  # 至少留一空行
    return "\n".join(lines)

def _format_positions_block(symbols, coord_type, coords):
    lines = ["ATOM_POSITIONS",
             f"{len(symbols)}                               ! number of atoms for projectors",
             "Direct" if coord_type.lower().startswith("d") else "Cartesian"]
    for sym, (x, y, z) in zip(symbols, coords):
        lines.append(f"{sym:<3} {x:.16f}    {y:.16f}    {z:.16f} ")
    lines.append("")  # 结束空行
    return "\n".join(lines)

# ===== win 投影解析 & 格式化 =====

def _parse_win_projections(win_path):
    """返回 per_atom_orbitals: List[List[str]]，顺序: s | px py pz | dz2 dxz dyz dx2-y2 dxy"""
    if not os.path.isfile(win_path):
        raise FileNotFoundError(f"未找到 {win_path}")

    with open(win_path, "r", encoding="utf-8", errors="ignore") as f:
        txt = f.read()

    m = re.search(r"(?is)Begin\s+Projections(.*?)End\s+Projections", txt)
    if not m:
        raise ValueError("未在 wannier90.win 中找到 Begin/End Projections")

    block = m.group(1)
    lines = [ln.strip() for ln in block.splitlines() if ln.strip()]

    p_orbs = ["px", "py", "pz"]
    d_orbs = ["dz2", "dxz", "dyz", "dx2-y2", "dxy"]

    per_atom_orbitals = []
    for ln in lines:
        ls = re.findall(r"l\s*=\s*([0-2])", ln)
        chosen = []
        if "0" in ls: chosen += ["s"]
        if "1" in ls: chosen += p_orbs
        if "2" in ls: chosen += d_orbs
        if not chosen:
            chosen = ["s"]
        per_atom_orbitals.append(chosen)
    return per_atom_orbitals

def _format_projectors_block(expanded_symbols, per_atom_orbitals):
    n = min(len(expanded_symbols), len(per_atom_orbitals))

    counts = []
    for orbs in per_atom_orbitals[:n]:
        c = 0
        if "s" in orbs: c += 1
        if any(o in orbs for o in ("px","py","pz")): c += 3
        if any(o in orbs for o in ("dz2","dxz","dyz","dx2-y2","dxy")): c += 5
        counts.append(c)
    head = " ".join(str(x) for x in counts)

    lines = ["PROJECTORS", f" {head}                  ! number of projectors"]
    for sym, orbs in zip(expanded_symbols[:n], per_atom_orbitals[:n]):
        row = [sym]
        if "s" in orbs: row.append("s")
        if any(o in orbs for o in ("px","py","pz")): row += ["px","py","pz"]
        if any(o in orbs for o in ("dz2","dxz","dyz","dx2-y2","dxy")): row += ["dz2","dxz","dyz","dx2-y2","dxy"]
        lines.append(" " + " ".join(row))
    lines.append("")
    return "\n".join(lines)

# ===== 总入口（一次性完成所有操作） =====

def write_wt_in_from_poscar(material_dir, template_wt_in, out_root):
    """
    - 若 out_root/<材料名> 目录不存在：返回 False, "<mat> 没有 hr"
    - 否则：
        * 复制模板 wt.in 到该目录
        * 写入/替换 LATTICE, ATOM_POSITIONS, PROJECTORS
        * 从 vasprun 读取 E_FERMI；同时替换 Hrfile='<mat>_hr.dat'
    返回: (ok, 输出文件路径 | 错误信息)
    """
    mat = os.path.basename(material_dir.rstrip(os.sep))
    target_dir = os.path.join(out_root, mat)

    # 你要求：目录不存在则直接返回，不创建
    if not os.path.exists(target_dir):
        return False, f"{mat} 没有 hr"

    wt_out = os.path.join(target_dir, "wt.in")
    shutil.copy2(template_wt_in, wt_out)

    # 读取模板
    with open(wt_out, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    # 读取 POSCAR/CONTCAR → LATTICE & ATOM_POSITIONS
    poscar_path = None
    for name in ("POSCAR", "CONTCAR"):
        p = os.path.join(material_dir, name)
        if os.path.isfile(p):
            poscar_path = p
            break
    if poscar_path is None:
        return False, f"未找到 POSCAR/CONTCAR 于 {material_dir}"

    try:
        symbols, counts, coord_type, lattice, coords = _parse_poscar(poscar_path)
    except Exception as e:
        return False, f"解析POSCAR失败: {e}"

    lattice_block = _format_lattice_block(lattice)
    positions_block = _format_positions_block(symbols, coord_type, coords)
    content = _replace_section(content, "LATTICE", lattice_block)
    content = _replace_section(content, "ATOM_POSITIONS", positions_block)

    # 从 vasprun 读取费米能
    vasprun_path = None
    for fn in ("vasprun.xml", "vasprun.xml.gz"):
        p = os.path.join(material_dir, p := os.path.join(material_dir, fn))
        # 修正写法
    vasprun_path = None
    for fn in ("vasprun.xml", "vasprun.xml.gz"):
        p = os.path.join(material_dir, fn)
        if os.path.isfile(p):
            vasprun_path = p
            break
    if vasprun_path is None:
        return False, "未找到 vasprun.xml 或 vasprun.xml.gz"
    vr = Vasprun(vasprun_path, parse_projected_eigen=False)
    efermi = vr.efermi

    # Hrfile / E_FERMI
    content = _replace_key_value_lines(content, f"{mat}_hr.dat", efermi)

    # PROJECTORS：从 wannier90.win 解析
    win_path = os.path.join(material_dir, "WR", "wannier90.win")
    per_atom_orbitals = _parse_win_projections(win_path)
    # 用 POSCAR 展开的元素名对齐
    expanded_symbols = symbols  # 已经是展开后的
    proj_block = _format_projectors_block(expanded_symbols, per_atom_orbitals)
    content = _replace_section(content, "PROJECTORS", proj_block)

    # 关键块间加空行
    content = _ensure_extra_blank_line_between_blocks(content)

    with open(wt_out, "w", encoding="utf-8") as f:
        f.write(content)

    return True, wt_out

def find_vasprun(material_dir):
    for fn in ("vasprun.xml", "vasprun.xml.gz"):
        p = os.path.join(material_dir, fn)
        if os.path.isfile(p):
            return p
    return None

def copy_hr(material_dir, hrs_dir, rename=True):
    src = os.path.join(material_dir, "WR", "wannier90_hr.dat")
    src_wosoc = os.path.join(material_dir, "WR", "spin_wo_soc_hr.dat")
    os.makedirs(hrs_dir, exist_ok=True)


    if (os.path.isfile(src) ) :
        mat = os.path.basename(material_dir.rstrip(os.sep))
        # 为每个材料新建一个子目录
        target_dir = os.path.join(hrs_dir, mat)
        os.makedirs(target_dir, exist_ok=True)
        dst = os.path.join(target_dir, f"{mat}_hr.dat" if rename else "wannier90_hr.dat")
        shutil.copy2(src, dst)
    elif ( os.path.isfile(src_wosoc) ):
        mat = os.path.basename(material_dir.rstrip(os.sep))
        # 为每个材料新建一个子目录
        target_dir = os.path.join(hrs_dir, mat)
        os.makedirs(target_dir, exist_ok=True)
        dst = os.path.join(target_dir, f"{mat}_hr.dat" if rename else "wannier90_hr.dat")
        shutil.copy2(src_wosoc, dst)    
    else :
        return False, f"缺少 {src}"
    

    return True, dst

# ===================== DOSCAR 判据（优先用于提示/佐证） =====================

def get_fermi_from_outcar(outcar_path):
    """从 OUTCAR 抓取第一处 E-fermi 数值；失败返回 None"""
    try:
        with open(outcar_path, "r", errors="ignore") as f:
            for line in f:
                if "E-fermi" in line:
                    toks = line.strip().split()
                    try:
                        return float(toks[2])
                    except Exception:
                        m = re.search(r"E-fermi\s*:\s*([\-0-9.]+)", line)
                        if m:
                            return float(m.group(1))
                        return None
    except Exception:
        return None
    return None

def parse_total_dos_from_doscar(doscar_path):
    """
    读取 DOSCAR 的总 DOS（合并自旋），假定第6行第3列就是 NEDOS。
    返回 energies(list), total_dos(list)
    """
    with open(doscar_path, "r", errors="ignore") as f:
        lines = [ln for ln in f if ln.strip()]

    if len(lines) < 6:
        raise RuntimeError("DOSCAR 太短")

    header = lines[5].split()
    ngrid = int(float(header[2]))  # 第6行第3列 = NEDOS


    data = lines[6:6+ngrid]

    energies, total_dos = [], []
    for ln in data:
        toks = ln.split()
        if not toks:
            continue
        E = float(toks[0])
        if len(toks) >= 5:
            # 自旋极化
            DOS = float(toks[1]) + float(toks[2])
        elif len(toks) >= 3:
            # 非自旋
            DOS = float(toks[1])
        else:
            continue
        energies.append(E)
        total_dos.append(DOS)

    return energies, total_dos

def dos_at_efermi(energies, total_dos, efermi):
    """
    线性插值估算 DOS(E_F)。如果 E_F 超出网格范围，则取最近端点。
    """
    if not energies:
        return None
    if energies[0] > energies[-1]:

        energies = energies[::-1]; total_dos = total_dos[::-1]
    idx = bisect.bisect_left(energies, efermi)
    if idx == 0:
        
        return total_dos[0]
    if idx >= len(energies):
        return total_dos[-1]
    x0, x1 = energies[idx-1], energies[idx]
    y0, y1 = total_dos[idx-1], total_dos[idx]
    if x1 == x0:

        return 0.5 * (y0 + y1)
    t = (efermi - x0) / (x1 - x0)
    
    return y0 + t * (y1 - y0)

def doscar_judgement(material_dir, dos_threshold):
    """
    返回 (is_metal_by_dos: True/False/None, detail: str)
    None 表示缺少文件或解析失败。
    """
    doscar = os.path.join(material_dir, "DOSCAR")
    outcar = os.path.join(material_dir, "OUTCAR")
    if not (os.path.isfile(doscar) and os.path.isfile(outcar)):
        return None, "缺少 DOSCAR/OUTCAR"

    efermi = get_fermi_from_outcar(outcar)
    # print(efermi)
    if efermi is None:
        return None, "OUTCAR 中未找到 E-fermi"

    try:
        energies, tdos = parse_total_dos_from_doscar(doscar)
        d_ef = dos_at_efermi(energies, tdos, efermi)
        print(d_ef)
        if d_ef is None:
            return None, "DOSCAR 解析失败"
        return (d_ef > dos_threshold), f"DOS(EF)={d_ef:.6g} states/eV"
    except Exception as e:
        return None, f"DOSCAR 解析异常：{e}"

# ===================== BandStructure 判据（核心是否穿越 EF） =====================

def band_is_metal(vasprun_path):
    """
    使用 band structure 判据：有无能带跨越 EF（pymatgen 的 BandStructure.is_metal）
    返回 (True/False/None, detail: str)
    """
    try:
        vr = Vasprun(vasprun_path, parse_dos=False, parse_projected_eigen=False)
        bs = vr.get_band_structure()
        return bs.is_metal(), "BandStructure.is_metal()"
    except Exception as e:
        return None, f"BandStructure 构建失败：{e}"

# ===================== 主流程 =====================


if __name__ == "__main__":
    main()