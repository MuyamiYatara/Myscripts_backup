import subprocess as sp
import re
import copy
import os
import time


#------------------functions------------------#
# All edit function can only edit the parameter the first time it occurs !!!!
def edit_fatband_generate(param_dict):
    # 读取文件内容
    with open(f'./fatband_generate.py', 'r') as f:
        lines = f.readlines()
        #print(lines)
    
    
    lines_edit = copy.deepcopy(lines)
    write_flag = False
    # 遍历字典中的每个参数
    for param_name, new_value in param_dict.items():
        write_flag = True
        # 遍历每一行，查找并修改多个参数
        for i, line in enumerate(lines):
            # 匹配当前参数的行
            
            match = re.search(rf"{param_name}", line)
            ifdelete = False
            if ((new_value=='d')or(new_value=='delete')):
                ifdelete = True
                write_flag = False
            if ifdelete :
                if match :
                    lines_edit[i] = " "
                    break
            else :
                if match :
                    lines_edit[i] = f"{param_name} = {new_value}\n"
                    write_flag = False
                    break        
    if write_flag :
        lines_edit.append(f"{param_name} = {new_value}\n")

    with open(f'./fatband_generate.py', 'w') as f:
        f.writelines(lines_edit)
    # 写回修改后的内容

    print("fatband_generate.py脚本参数修改完成。")

def edit_sbatch_script(param_dict):
    # 读取文件内容
    with open(f'./{sbatch_script_name}', 'r') as f:
        lines = f.readlines()
        #print(lines)
    
    
    lines_edit = copy.deepcopy(lines)
    write_flag = False
    # 遍历字典中的每个参数
    for param_name, new_value in param_dict.items():
        write_flag = True
        # 遍历每一行，查找并修改多个参数
        for i, line in enumerate(lines):
            # 匹配当前参数的行
            
            match = re.search(rf"{param_name}", line)
            ifdelete = False
            if ((new_value=='delete')):
                ifdelete = True
                write_flag = False
            
            if ifdelete :
                if match :
                    lines_edit[i] = " "
                    break
            else :
                if match :
                    lines_edit[i] = f"{param_name} {new_value}\n"
                    write_flag = False
                    break        
    if write_flag :
        lines_edit.append(f"{param_name} {new_value}\n")



            
    with open(f'./{sbatch_script_name}', 'w') as f:
        f.writelines(lines_edit)
    # 写回修改后的内容

    print("sbatch脚本参数修改完成。")

def edit_INCAR(param_dict):
    # 读取文件内容
    with open('./INCAR', 'r') as f:
        lines = f.readlines()
        #print(lines)
    
    
    lines_edit = copy.deepcopy(lines)
    write_flag = False
    # 遍历字典中的每个参数
    for param_name, new_value in param_dict.items():
        write_flag = True
        # 遍历每一行，查找并修改多个参数
        for i, line in enumerate(lines):
            # 匹配当前参数的行
            
            match = re.search(rf"{param_name}", line)
            ifdelete = False
            if ((new_value=='d')or(new_value=='delete')):
                ifdelete = True
                write_flag = False

            if ifdelete :
                if match :
                    lines_edit[i] = " "
                    break
            else :
                if match :
                    lines_edit[i] = f"{param_name} = {new_value}\n"
                    write_flag = False
                    break        
        if write_flag :
            lines_edit.append(f"{param_name} = {new_value}\n")


            
    with open('./INCAR', 'w') as f:
        f.writelines(lines_edit)
    # 写回修改后的内容

    print("INCAR参数修改完成。")

def get_ENCUT_value():
    with open(f"./POTCAR", 'r') as f:
        lines = f.readlines()
    enmax_values = []
    for i, line in enumerate(lines):        
        match = re.search(r"ENMAX\s*=\s*(\d+(\.\d+)?)", line)
        if match :
            enmax_values.append(match.group(1))
    enmax_values = [float(i) for i in enmax_values]
    encut = round(1.5*max(enmax_values)) 
    return str(encut)

def get_fermi_energy():
    # 正则表达式匹配 E-fermi 行，并提取数字部分
    pattern = r"E-fermi\s*:\s*([+-]?\d*\.\d+|\d+)"
    
    with open('OUTCAR', 'r') as file:
        for line in file:
            match = re.search(pattern, line)
            if match:
                fermi_energy = match.group(1)  # 提取匹配的数字部分
                return float(fermi_energy)  # 转换为浮点数并返回
    return None  # 如果没有找到匹配项

def modify_and_copy_file(input_file, output_file, ef):
    # 正则表达式模式，匹配 "set arrow from xxx, xxx to xxx, xxx nohead" 这种行
    arrow_pattern = r"set arrow from\s*(-?\d+\.\d+),\s*(-?\d+\.\d+)\s*to\s*(-?\d+\.\d+),\s*(-?\d+\.\d+)\s*nohead"
    xtics_pattern = r"set xtics \((.*?)\)"  # 匹配 set xtics 行

    with open(input_file, 'r') as file_a:
        lines = file_a.readlines()

    modified_lines = []  # 用来存放修改后的行

    # 遍历文件中的每一行
    for line in lines:
        # 查找 "set arrow" 行并进行修改
        match = re.match(arrow_pattern, line)
        if match:
            # 提取数字部分并替换
            start_x = match.group(1)
            start_y = match.group(2)
            end_x = match.group(3)
            end_y = match.group(4)
            
            # 替换为 lower_bound 和 upper_bound
            modified_line = f"set arrow from {start_x}, lower_bound to {end_x}, upper_bound nohead\n"
            modified_lines.append(modified_line)
        else:
            # 如果是 set xtics 行，直接添加到修改后的行列表
            match_xtics = re.match(xtics_pattern, line)
            if match_xtics:
                # 直接添加 xtics 行
                modified_lines.append(line)

    # 读取文件B的内容
    with open(output_file, 'r') as file_b:
        original_lines = file_b.readlines()

    # 计算倒数第三行的索引
    insert_index = -2

    # 在倒数第3行插入修改后的内容
    new_lines = original_lines[:insert_index] + modified_lines + original_lines[insert_index:]

    # 查找并修改fermi能
    pattern = r"(\$2 - )(\d+\.\d+)"

    fermi_str = "'" + str(ef) + "'"
    replaced_text = re.sub(pattern, r"\1" + fermi_str, new_lines[len(new_lines)-1])
    new_lines[len(new_lines)-1] = replaced_text
    # 将修改后的内容写回文件B
    with open(output_file, 'w') as file_b:
        file_b.writelines(new_lines)
    print("w9.gnu 修改完成")

def read_KPOINTS():
    with open("KPOINTS", 'r') as f:
        lines = f.readlines()

    # 筛选出含有坐标的行（排除空行和其他非数据行）
    kpoints_data = []
    for line in lines:
        kpoints_data = [line.strip() for line in lines[4:] if line.strip()]

    # 确保 kpoints_data 的长度是偶数
    if len(kpoints_data) % 2 != 0:
        raise ValueError("KPOINTS 文件的格式不完整，无法两两配对。")

    # 转换为目标格式
    formatted_lines = []
    for i in range(0, len(kpoints_data), 2):
        start = re.match(r"([\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+)\s+(\w+)", kpoints_data[i])
        end = re.match(r"([\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+)\s+(\w+)", kpoints_data[i+1])
        formatted_lines.append(f"{start.group(2)} {start.group(1)}     {end.group(2)} {end.group(1)}")


    return formatted_lines

def read_POSCAR():
    with open("POSCAR", 'r') as f:
        lines = f.readlines()
    
    lines = [i.strip() for i in lines]
    # 筛选出含有坐标的行（排除空行和其他非数据行）
    lattice_data = lines[2:5]
    atom_data = lines[8:]
    formatted_lines = []
    for i in range(0, len(atom_data)):
        #注意这里要求POSCAR的原子坐标后面必须接有原子名称，不然正则表达式会匹配不到
        start = re.match(r"([\d\.\-]+\s+[\d\.\-]+\s+[\d\.\-]+)\s+(\w+)", atom_data[i])
        formatted_lines.append(f"{start.group(2)} {start.group(1)}")
    
    return lattice_data, formatted_lines 

def edit_win(lattice_vec, atom_pos, kpath):
    with open("wannier90.win_backup", "r") as file:
        content = file.read()

    pattern_unit_cell_cart = re.compile(r"(begin unit_cell_cart)(.*?)(end unit_cell_cart)", re.DOTALL)
    pattern_atoms_frac = re.compile(r"(begin atoms_frac)(.*?)(end atoms_frac)", re.DOTALL)
    pattern_kpoint_path = re.compile(r"(begin kpoint_path)(.*?)(end kpoint_path)", re.DOTALL)
    lattice_vec_write = ""
    atom_pos_write = ""
    kpath_write = ""

    for i in lattice_vec :
        lattice_vec_write = lattice_vec_write + i + "\n"
    
    for i in atom_pos :
        atom_pos_write = atom_pos_write + i + "\n"
    
    for i in kpath :
        kpath_write = kpath_write + i + "\n"

    def replace_unit_cell_cart(match):
        return f"{match.group(1)}\n{lattice_vec_write}{match.group(3)}\n"

    def replace_atoms_frac(match):
        return f"{match.group(1)}\n{atom_pos_write}{match.group(3)}\n"
    
    def replace_kpath(match):
        return f"{match.group(1)}\n{kpath_write}{match.group(3)}\n"
    
    content = pattern_unit_cell_cart.sub(replace_unit_cell_cart, content)
    content = pattern_atoms_frac.sub(replace_atoms_frac, content)
    content = pattern_kpoint_path.sub(replace_kpath, content)

    with open("wannier90.win_backup", "w") as file:
        file.write(content)

    print("wannier90.win基本参数填写完成。")

def replace_hr_plot():
    """
    将文件中所有 'hr_plot' 替换为 'write_hr'。
    """
    file_path = "wannier90.win"
    with open(file_path, 'r') as f:
        content = f.read()

    new_content = content.replace('hr_plot', 'write_hr')
    with open(file_path, 'w') as f:
        f.write(new_content)


def submit_sbatch_script(script_path):
    """Submit a job with sbatch and return the job ID."""
    # 提交脚本并捕获输出
    result = sp.run(["sbatch", script_path], stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    
    if result.returncode != 0:
        print("submit unsuccessfully:", result.stderr)
        return None
    
    # 提取作业ID
    job_id = result.stdout.split()[-1]  # sbatch 输出通常以 "Submitted batch job <job_id>" 形式返回
    print(f"submit successfully,job ID: {job_id}")
    return job_id


def monitor_job(job_id):
    """Monitor the status of a job using squeue and display runtime."""
    start_time = time.time()  # 记录作业提交的开始时间
    
    while True:
        # 使用 squeue 查看作业状态
        result = sp.run(["squeue", "-j", job_id], stdout=sp.PIPE, stderr=sp.PIPE, text=True)
        
        
        # squeue 输出作业状态信息，若作业完成，squeue 不会显示该作业
        if job_id not in result.stdout:
            print(f"job {job_id} completed.")
            break
        
        # 计算已运行的时间
        elapsed_time = time.time() - start_time
        elapsed_minutes = elapsed_time // 60
        elapsed_seconds = int(elapsed_time % 60)

        status = result.stdout.split()[12]
        print(f"job {job_id} is running... running time: {int(elapsed_minutes)}分 {elapsed_seconds}秒  status: {status}")
        
        time.sleep(5)

#-#--------------------------------generate the input files--------------------------------#-#
#------------------self consistent------------------#
##### some input parameters #####
kmesh_accurary_level = "0.03"
sbatch_script_name = 'NANO.sh'



# input commands
nano = sp.Popen(["cp", f"/data/home/ycshen/Myscripts/dft-tools/{sbatch_script_name}", "./"])
vaspkit = sp.Popen(["vaspkit"], stdin=sp.PIPE)

# vaspkit task 102, Monkhorst-Pack Scheme
command_sc_inputs = f"102\n1\n{kmesh_accurary_level}\n".encode()

# copy the sbatch script and generate the basic input files
nano.communicate()
vaspkit.communicate(input=command_sc_inputs)

##### edit the sbatch script #####
sbatch_para = {
    '#SBATCH -N': 1,
    '#SBATCH -n': 56,
    '#SBATCH -A': 'hmt03',  
    '#SBATCH -p': 'regular,regular6430',
}
edit_sbatch_script(sbatch_para)


# get some parameters
encut = get_ENCUT_value()

##### edit the INCAR for SCF #####
SCF_para = {
    'ISTART': 'delete',
    'NSW': 'delete',
    'IBRION': 'delete',
    'ISIF': 'delete',
    'EDIFFG': 'delete',
    'ENCUT': encut,
    "LORBMOM": ".TRUE.",
    "LWAVE":".FALSE.",
    'ISPIN': '1',
    'LSORBIT': '.TRUE.',
    'MAGMOM': '0 0 0 0 0 0  2 2 2 2 2 2  0'
    # "NPAR":"32",
    'NBANDS':'112'

}
edit_INCAR(SCF_para)


#------------------band structure and fatband------------------#
# change the work directory
sp.run(['mkdir','-p','BD'])
os.chdir('./BD')

# input commands
#cd_BD = sp.Popen(["cd", './BD'])
cp_poscar = sp.Popen(["cp", f"../POSCAR", "./"])
cp_incar = sp.Popen(["cp", f"../INCAR", "./"])
cp_potcar = sp.Popen(["cp", f"../POTCAR", "./"])
cp_nano = sp.Popen(["cp", f"../{sbatch_script_name}", "./"])
fatband = sp.Popen(["cp", f"/data/home/ycshen/Myscripts/dft-tools/fatband_generate.py", "./"])

vaspkit = sp.Popen(["vaspkit"], stdin=sp.PIPE)

command_bd_inputs = f"303\n".encode()

# make bandstructure calculation directory and copy files 


#cd_BD.communicate()
fatband.communicate()
cp_incar.communicate()
cp_poscar.communicate()
cp_potcar.communicate()
cp_nano.communicate()

# edit INCAR 
BD_para = {
    'ISMEAR':'delete',
    'SIGMA':'delete',
    'NELM':'delete',
    'NELMIN':'delete',
    'EDIFF':'delete',
    'ICHARG':'11',
    'LORBIT':'11',


}
edit_INCAR(BD_para)

# generate KPOINTS for band calculation
vaspkit.communicate(command_bd_inputs)
sp.run(['cp', 'KPATH.in', 'KPOINTS'])

# modify the number of kpoints 
knum = '    70'
with open('KPOINTS', 'r') as f:
    lines = f.readlines()
if len(lines) >= 2:
    lines[1] = knum.rstrip('\n') + '\n'  # 替换第二行，确保只有一个换行符
with open('KPOINTS', 'w', encoding='utf-8') as f:
    f.writelines(lines)


# edit fatband_generate.py, without fermi energy

# 初始化列表
elements = []  # 第一个列表
orbitals = []  # 第二个列表

# 打开文件并读取第6行和第7行
with open('POSCAR', "r") as file:
    lines = file.readlines()
    element_names = lines[5].strip().split()  # 第6行：元素名称
    
# 生成第一个列表
for name in element_names:
    elements.extend([name] * 3 )  # 每个元素名称重复 4 次

# 生成第二个列表
total_elements = len(element_names)  # 总元素个数
orbitals = ["s", "p", "d"] * total_elements  # 循环 s, p, d


# 如果有稀土元素需要加入f轨道，则换成下面几行
#for name in element_names:
#    elements.extend([name] * 4 )  # 每个元素名称重复 4 次
#
#total_elements = len(element_names)  # 总元素个数
#orbitals = ["s", "p", "d", "f"] * total_elements  # 循环 s, p, d, f


fatband_para = {
    'modeofplot':"'o'",
    'targetelements':str(elements),
    'targetorbitals':str(orbitals),
    'Emin':"'-15'",
    'Emax':"'30'"
}


edit_fatband_generate(fatband_para)

# get the kpoints and poscar information and store for wannier90.win
klines = read_KPOINTS()
lattice, atoms = read_POSCAR()

# change back to the SCF directory
os.chdir('..')
#------------------wannier test------------------#
# change the work directory
sp.run(['mkdir','-p','WR'])
os.chdir('./WR')


# input commands
#cd_BD = sp.Popen(["cd", './BD'])
cp_poscar = sp.Popen(["cp", f"../POSCAR", "./"])
cp_incar = sp.Popen(["cp", f"../INCAR", "./"])
cp_potcar = sp.Popen(["cp", f"../POTCAR", "./"])
cp_kpoints = sp.Popen(["cp", f"../KPOINTS", "./"])
cp_nano = sp.Popen(["cp", f"../{sbatch_script_name}", "./"])
win = sp.Popen(["cp", f"/data/home/ycshen/Myscripts/dft-tools/wannier90.win_backup", "./"])
w90 = sp.Popen(["cp", f"/data/home/ycshen/Myscripts/dft-tools/W90.sh", "./"])
w9gnu = sp.Popen(["cp", f"/data/home/ycshen/Myscripts/dft-tools/w9.gnu", "./"])

win.communicate()
w90.communicate()
w9gnu.communicate()
cp_incar.communicate()
cp_poscar.communicate()
cp_potcar.communicate()
cp_kpoints.communicate()
cp_nano.communicate()

# edit INCAR 
WR_para = {
    'ISMEAR':'delete',
    'SIGMA':'delete',
    'NELM':'delete',
    'NELMIN':'delete',
    'EDIFF':'delete',
    'NPAR':'delete',
    'KPAR':'delete',
    'ICHARG':'11',
    'LWANNIER90':'.TRUE.',
    'ISYM':'-1'
}
edit_INCAR(WR_para)

# give a basic wannier90.win

edit_win(lattice, atoms, klines)

# change back to the SCF directory
os.chdir('..')


#-#---------------------------generate the input files---------------------------#-#

#-#---------------------------submit and watch---------------------------#-#
scf = 1
bd = 1
wrscf = 11
wr = 11
#efermi = get_fermi_energy()
#print(efermi)

if (scf == 1) :
    jobid = submit_sbatch_script(sbatch_script_name)
    monitor_job(jobid)
    print("SCF completed")


if (bd == 1) :
    efermi = get_fermi_energy()
    sp.run(["cp","CHGCAR","./BD"])
    os.chdir('./BD')
    jobid = submit_sbatch_script(sbatch_script_name)
    monitor_job(jobid)

    fatband_para={'fermienergy':"'"+str(efermi)+"'"}
    edit_fatband_generate(fatband_para)
    sp.run(["python", "fatband_generate.py"])
    sp.run(["python", "fatband.py"])

    vaspkit = sp.Popen(["vaspkit"], stdin=sp.PIPE)
    command_bd_inputs = f"211\n".encode()
    vaspkit.communicate(command_bd_inputs)
    print("BD completed")
    os.chdir('..')

#------------------wannier90.win prepared------------------#
if (wrscf == 1) :
    sp.run(["cp","CHGCAR","./WR"])
    os.chdir('./WR')
    jobid = submit_sbatch_script(sbatch_script_name)
    monitor_job(jobid)
    print("WR SCF completed")
    os.chdir('..')

if (wr == 1) :
    efermi = get_fermi_energy()
    os.chdir('./WR')
    sbatch_script_name = 'W90.sh'
    sbatch_para = {
        '#SBATCH -N': 1,
        '#SBATCH -n': 64,
        '#SBATCH -A': 'hmt03',  
        '#SBATCH -p': 'regular6430',
    }
    edit_sbatch_script(sbatch_para)
    jobid = submit_sbatch_script(sbatch_script_name)
    replace_hr_plot() #在这里用wannier90 v3.1并行运行wannier90,但前面是wannier90 v2.1的接口,因此要把hr_plot改成write_hr
    monitor_job(jobid)

    input_file = 'wannier90_band.gnu'
    output_file = 'w9.gnu'
    modify_and_copy_file(input_file, output_file, efermi)

    sp.run(["gnuplot", "w9.gnu"])
    sbatch_script_name = 'NANO.sh'

    os.chdir('..')




#-#---------------------------submit and watch---------------------------#-#
