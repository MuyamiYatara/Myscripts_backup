import os
import re
import subprocess as sp


#-----------------------提取磁矩-----------------------#
def get_magmom(poscar_file):
    """
    从 POSCAR（可为超胞或超胞）读取每个原子的磁矩并按元素分类返回字典。
    返回值示例：
      {
        "Co": [["0", "0", "2.6"], ["0", "0", "2.6"]],
        "F" : [["0", "0", "0"], ["0", "0", "0"], ...],
        ...
      }
    注：本函数不对磁矩做数值变换，仅按原文件读取后三列（若存在）。
    """
    with open(poscar_file, "r") as f:
        # 去掉空行（保持每行的内部空格）
        lines = [line.strip() for line in f if line.strip()]

    # 找到 Direct/Cartesian 行的位置（i 为该行索引）
    direct_idx = None
    for i, line in enumerate(lines):
        low = line.lower()
        if low.startswith("direc") or low.startswith("cartesian"):
            direct_idx = i
            break
    if direct_idx is None:
        raise ValueError("POSCAR 文件中找不到 Direct 或 Cartesian 行")

    # 从 Direct 行上方向上寻找 "numbers" 行（只含数字和空白的行）
    numbers_idx = None
    for j in range(direct_idx - 1, -1, -1):
        # 若该行至少含有一个数字，并且不含字母（排除元素名行）
        if re.search(r"\d", lines[j]) and not re.search(r"[A-Za-z]", lines[j]):
            # 进一步确认行中主要是数字（允许空格分隔）
            if re.match(r"^\s*[\d\s]+\s*$", lines[j]):
                numbers_idx = j
                break
    if numbers_idx is None:
        # 兼容性回退（按以前的策略：Direct 上方两行）
        numbers_idx = direct_idx - 2
        if numbers_idx < 0:
            raise ValueError("未能定位到原子数行（numbers line）")

    # 元素名行应在 numbers_idx 的上一行
    elements_idx = numbers_idx - 1
    if elements_idx < 0:
        raise ValueError("未能定位到元素名行")

    elem_line = lines[elements_idx]
    num_line = lines[numbers_idx]

    # 解析元素与对应的数量（数量行用正则提取整数）
    elements = elem_line.split()
    numbers = [int(x) for x in re.findall(r"\d+", num_line)]
    if len(elements) != len(numbers):
        raise ValueError(f"元素名数量与数目不匹配：elements={elements}, numbers={numbers}")

    # 构建原子名序列（按元素顺序展开）
    atom_names = []
    for el, cnt in zip(elements, numbers):
        atom_names.extend([el] * cnt)
    natoms = sum(numbers)

    # 读取坐标（从 Direct/Cartesian 下一行开始）
    coord_start = direct_idx + 1
    coords_with_mag = lines[coord_start: coord_start + natoms]
    if len(coords_with_mag) < natoms:
        raise ValueError("POSCAR 中的坐标行数小于期望的原子数")

    # 构建字典：元素 -> list of mags (每个为 list(3) 或 [])
    mags_dict = {el: [] for el in elements}
    for atom_el, line in zip(atom_names, coords_with_mag):
        parts = line.split()
        if len(parts) >= 6:
            # 转为 float，取绝对值，保留 1 位小数
            mags = [f"{float(x):.1f}" for x in parts[3:6]]
            mags_dict[atom_el].append(mags)
        else:
            mags_dict[atom_el].append([])

    return mags_dict

def get_magmom_abs(poscar_file) :
    """
    从 POSCAR（可为超胞或超胞）读取每个原子的磁矩并按元素分类返回字典。
    返回值示例：
      {
        "Co": [["0", "0", "2.6"], ["0", "0", "2.6"]],
        "F" : [["0", "0", "0"], ["0", "0", "0"], ...],
        ...
      }
    注：本函数取磁矩数值的绝对值。
    """
    with open(poscar_file, "r") as f:
        # 去掉空行（保持每行的内部空格）
        lines = [line.strip() for line in f if line.strip()]

    # 找到 Direct/Cartesian 行的位置（i 为该行索引）
    direct_idx = None
    for i, line in enumerate(lines):
        low = line.lower()
        if low.startswith("direc") or low.startswith("cartesian"):
            direct_idx = i
            break
    if direct_idx is None:
        raise ValueError("POSCAR 文件中找不到 Direct 或 Cartesian 行")

    # 从 Direct 行上方向上寻找 "numbers" 行（只含数字和空白的行）
    numbers_idx = None
    for j in range(direct_idx - 1, -1, -1):
        # 若该行至少含有一个数字，并且不含字母（排除元素名行）
        if re.search(r"\d", lines[j]) and not re.search(r"[A-Za-z]", lines[j]):
            # 进一步确认行中主要是数字（允许空格分隔）
            if re.match(r"^\s*[\d\s]+\s*$", lines[j]):
                numbers_idx = j
                break
    if numbers_idx is None:
        # 兼容性回退（按以前的策略：Direct 上方两行）
        numbers_idx = direct_idx - 2
        if numbers_idx < 0:
            raise ValueError("未能定位到原子数行（numbers line）")

    # 元素名行应在 numbers_idx 的上一行
    elements_idx = numbers_idx - 1
    if elements_idx < 0:
        raise ValueError("未能定位到元素名行")

    elem_line = lines[elements_idx]
    num_line = lines[numbers_idx]

    # 解析元素与对应的数量（数量行用正则提取整数）
    elements = elem_line.split()
    numbers = [int(x) for x in re.findall(r"\d+", num_line)]
    if len(elements) != len(numbers):
        raise ValueError(f"元素名数量与数目不匹配：elements={elements}, numbers={numbers}")

    # 构建原子名序列（按元素顺序展开）
    atom_names = []
    for el, cnt in zip(elements, numbers):
        atom_names.extend([el] * cnt)
    natoms = sum(numbers)

    # 读取坐标（从 Direct/Cartesian 下一行开始）
    coord_start = direct_idx + 1
    coords_with_mag = lines[coord_start: coord_start + natoms]
    if len(coords_with_mag) < natoms:
        raise ValueError("POSCAR 中的坐标行数小于期望的原子数")

    # 构建字典：元素 -> list of mags (每个为 list(3) 或 [])
    mags_dict = {el: [] for el in elements}
    for atom_el, line in zip(atom_names, coords_with_mag):
        parts = line.split()
        if len(parts) >= 6:
            # 转为 float，取绝对值，保留 1 位小数
            mags = [f"{abs(float(x)):.1f}" for x in parts[3:6]]
            mags_dict[atom_el].append(mags)
        else:
            mags_dict[atom_el].append([])

    return mags_dict



def write_magmom_nsoc(mags_dict, semi_path, poscar_file):
    """
    mags_dict: dict, {元素: [[mx,my,mz], ...]}，数值字符串或float都可
    semi_path: str, semiauto.py 的路径
    poscar_file: str, 原胞 POSCAR 文件路径
    """

    # Step 1: 读取原胞 POSCAR，获取元素和数量
    with open(poscar_file, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    elements = re.findall(r"\S+", lines[5])   # 第6行：元素名
    counts = [int(x) for x in re.findall(r"\S+", lines[6])]  # 第7行：数量

    # Step 2: 根据数量拼接磁矩字符串
    mags_flat = []
    for el, num in zip(elements, counts):
        if el not in mags_dict:
            raise ValueError(f"元素 {el} 在 mags_dict 中不存在！")

        # 取前 num 个磁矩
        sub_mags = mags_dict[el][:num]

        # 每个原子取三个方向中的最大值（绝对值），并保留1位小数
        for m in sub_mags:
            if m:  # 非空
                vals = [abs(float(x)) for x in m]
                mags_flat.append(f"{max(vals):.1f}")
            else:
                mags_flat.append("0.0")

    mag_str = " ".join(mags_flat)

    # Step 3: 修改 semi_path 文件中的 'MAGMOM' 行
    with open(semi_path, "r") as f:
        lines = f.readlines()

    new_lines = []
    pattern = re.compile(r"(^\s*'MAGMOM'\s*:\s*')(.+)('.*$)")
    for line in lines:
        m = pattern.match(line)
        if m:
            line = f"{m.group(1)}{mag_str}{m.group(3)}\n"
        new_lines.append(line)

    with open(semi_path, "w") as f:
        f.writelines(new_lines)


def write_magmom(mags_dict, semi_path, poscar_file):
    """
    mags_dict: dict, {元素: [[mx,my,mz], ...]}，数值字符串或float都可
    semi_path: str, semiauto.py 的路径
    poscar_file: str, 原胞 POSCAR 文件路径
    """

    # Step 1: 读取原胞 POSCAR，获取元素和数量
    with open(poscar_file, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    elements = re.findall(r"\S+", lines[5])   # 第6行：元素名
    counts = [int(x) for x in re.findall(r"\S+", lines[6])]  # 第7行：数量

    # Step 2: 根据数量拼接磁矩字符串
    mags_flat = []
    for el, num in zip(elements, counts):
        if el not in mags_dict:
            raise ValueError(f"元素 {el} 在 mags_dict 中不存在！")

        # 取前 num 个磁矩
        sub_mags = mags_dict[el][:num]

        # 每个原子写入三个方向的磁矩，绝对值并保留1位小数
        for m in sub_mags:
            if m:  # 非空
                vals = [abs(float(x)) for x in m]
                mags_flat.extend([f"{v:.1f}" for v in vals])
            else:
                mags_flat.extend(["0.0", "0.0", "0.0"])

    mag_str = " ".join(mags_flat)

    # Step 3: 修改 semi_path 文件中的 'MAGMOM' 行
    with open(semi_path, "r") as f:
        lines = f.readlines()

    new_lines = []
    pattern = re.compile(r"(^\s*'MAGMOM'\s*:\s*')(.+)('.*$)")
    for line in lines:
        m = pattern.match(line)
        if m:
            line = f"{m.group(1)}{mag_str}{m.group(3)}\n"
        new_lines.append(line)

    with open(semi_path, "w") as f:
        f.writelines(new_lines)


def write_nbands_lsorbit(socflag, semi_path, poscar_file, num_cpu=56):
    with open(semi_path, "r") as f:
        content = f.read()

    # 匹配 'LSORBIT': '任意内容'
    pattern = re.compile(r"('LSORBIT'\s*:\s*)'[^']*'")
    if socflag :
        content_new = pattern.sub(r"\1'.TRUE.'", content)
    else :
        content_new = pattern.sub(r"\1'.FALSE.'", content)

    
    #修改NBANDS
    with open(poscar_file, "r") as f:
        lines = [line.strip() for line in f if line.strip()]
    
    # 第 7 行是元素数目
    counts_line = lines[6]
    counts = [int(x) for x in re.findall(r"\d+", counts_line)]
    natoms = sum(counts)

    # Step 2: 计算 kn 和 NBANDS
    kn = (natoms * 9) // num_cpu
    new_nbands = num_cpu * (kn+1)

    # Step 3: 修改 semi_path 文件
    with open(semi_path, "r") as f:
        content = f.read()

    # 匹配 'NBANDS': 'xxxx'，只替换数字部分
    pattern = re.compile(r"('NBANDS'\s*:\s*')\d+(')")
    content_new = pattern.sub(lambda m: f"{m.group(1)}{new_nbands}{m.group(2)}", content_new)


    with open(semi_path, "w") as f:
        f.write(content_new)


if __name__ == "__main__":
    mags = get_magmom("./POSCAR")
    print(mags)
    write_magmom(mags, "./semiauto_dft_wannier.py", "./POSCAR_r")
    write_nbands_lsorbit(True, "./semiauto_dft_wannier.py", "./POSCAR_r", 56)