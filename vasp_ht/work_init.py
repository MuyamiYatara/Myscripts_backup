import os
import re
import subprocess as sp
import prepare


# Step 1: 拷贝文件 
# 原始 POSCAR 文件所在目录
source_dir = "/data/home/ycshen/openMX_wannier/multi_test/test_POSCARS"
# 工作目录的根目录
work_root = "/data/home/ycshen/openMX_wannier/multi_test/work_soc"

# 若只用修改semiauto文件，则将if_restart改为True可节省时间
restart = False

# 控制计算类型的参数
fm_only = True
soc_cal = True

# 提交任务时用的总cpu个数
num_cpu = 56

# POSCAR文件名匹配规则
pattern = re.compile(r"^([\d.]+)_(.+)\.POSCAR$")

if restart :
    materials_dirs = [os.path.join(work_root, d) for d in os.listdir(work_root)
                  if os.path.isdir(os.path.join(work_root, d))]
    for n in materials_dirs :
        os.chdir(n)
        # Step 3: 提取磁矩，拷贝并修改 semiautowannier.py
        sp.run(["cp", "/data/home/ycshen/Myscripts/dft-tools/semiauto_dft_wannier.py", "./"], check=True)
        if fm_only :
            mags = prepare.get_magmom_abs("./MPOSCAR")    
        else :
            mags = prepare.get_magmom("./MPOSCAR") 
        
        if soc_cal :
            prepare.write_magmom(mags, "./semiauto_dft_wannier.py", "./POSCAR")    
        else :
            prepare.write_magmom_nsoc(mags, "./semiauto_dft_wannier.py", "./POSCAR")
        prepare.write_nbands_lsorbit(soc_cal, "./semiauto_dft_wannier.py", "./POSCAR", num_cpu)       

else :
    # 遍历目录
    for fname in os.listdir(source_dir):
        match = pattern.match(fname)
        if match:
            number, material = match.groups()
            full_path = os.path.join(source_dir, fname)

            # 工作目录路径 = work_root/material
            work_dir = os.path.join(work_root, material)
            os.makedirs(work_dir, exist_ok=True)

            # 拷贝文件到新目录，并命名为 POSCAR
            target_file = os.path.join(work_dir, "POSCAR")
            sp.run(["cp", full_path, target_file], check=True)


            os.chdir(work_dir)
            # Step 3: 调用 phonopy 生成原胞 POSCAR
            sp.run(["conda", "run", "-n", "PY_NEW", "phonopy", "--symmetry", "-c", "POSCAR"], check=True)
            sp.run(["cp", "POSCAR", "MPOSCAR"])
            sp.run(["cp", "PPOSCAR", "POSCAR"])

            # Step 4: 提取磁矩，拷贝并修改 semiautowannier.py
            sp.run(["cp", "/data/home/ycshen/Myscripts/dft-tools/semiauto_dft_wannier.py", work_dir], check=True)
            if fm_only :
                mags = prepare.get_magmom_abs("./MPOSCAR")    
            else :
                mags = prepare.get_magmom("./MPOSCAR") 
            
            if soc_cal :
                prepare.write_magmom(mags, "./semiauto_dft_wannier.py", "./POSCAR")    
            else :
                prepare.write_magmom_nsoc(mags, "./semiauto_dft_wannier.py", "./POSCAR")  
            prepare.write_nbands_lsorbit(soc_cal, "./semiauto_dft_wannier.py", "./POSCAR", num_cpu)
            













