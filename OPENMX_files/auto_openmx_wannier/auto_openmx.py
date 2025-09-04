import subprocess as sp
import re
import copy
import os
import time
import sys


def get_mpi_np() -> int:
    """优先用 Slurm 提供的进程数，兜底用本机核数。"""
    for var in ("SLURM_NTASKS", "PMI_SIZE", "OMPI_UNIVERSE_SIZE"):
        v = os.environ.get(var)
        if v and v.isdigit():
            return int(v)

    # 交互/单节点时可能只有每节点数
    v = os.environ.get("SLURM_NTASKS_PER_NODE")
    if v:
        # 可能类似 "64(x1)" 或 "32,32"
        try:
            first = v.split(",")[0]
            first = first.split("(")[0]
            return int(first)
        except Exception:
            pass

    # 兜底：用机器可见核数
    try:
        out = sp.check_output(["nproc"], text=True).strip()
        return int(out)
    except Exception:
        return 1

def run_openmx(input_dat: str, out_log: str = None):
    np = get_mpi_np()
    env = os.environ.copy()
    # 避免 OpenMP 抢核（OpenMX+MPI 常见做法）
    env.setdefault("OMP_NUM_THREADS", "1")

    cmd = ["mpirun", "-np", str(np), "/home/users/shenyc/openmx/openmx3.9/source/openmx", input_dat]
    if out_log:
        with open(out_log, "w") as f:
            sp.run(cmd, check=True, env=env, stdout=f, stderr=sp.STDOUT)
    else:
        sp.run(cmd, check=True, env=env)

name = sys.argv[1]
poscar_dir = sys.argv[2]

#make work directory
print(name)
match = re.match(r"(\d+(?:\.\d+)?)_", name)
if match :
    material_id = match.group(1)
else :
    print("POSCAR name error!")
    exit()

dir_name = "mat-" + material_id
#sp.run(['mkdir','-p', dir_name])
os.chdir(dir_name)

#------------------nsoc scf------------------#
sp.run(['cp', poscar_dir+name, './']) ######
sp.run(['cp', '/home/users/shenyc/Myscripts/OPENMX_files/auto_openmx_wannier/openmx.dat_nsoc', 'openmx.dat_nsoc'])
sp.run(['cp', '/home/users/shenyc/Myscripts/OPENMX_files/auto_openmx_wannier/poscar2openmx.py', 'poscar2openmx.py'])
sp.run(["conda", "run", "-n", "mp_api", "python", "poscar2openmx.py", name, "standard_basis", "openmx.dat_nsoc", "openmx.dat"])
sp.run(['cp', '/home/users/shenyc/Myscripts/OPENMX_files/auto_openmx_wannier/HM.sh', 'HM.sh'])

run_openmx("openmx.dat", out_log="openmx_step1_nsoc.out")
#------------------soc scf------------------#
sp.run(['cp', 'openmx.dat', './openmx.dat_step1_nsoc'])
sp.run(['cp', '/home/users/shenyc/Myscripts/OPENMX_files/auto_openmx_wannier/openmx.dat_soc', 'openmx.dat_soc'])
sp.run(["conda", "run", "-n", "mp_api", "python", "poscar2openmx.py", name, "standard_basis", "openmx.dat_soc", "openmx.dat"])

run_openmx("openmx.dat", out_log="openmx_step2_soc.out")


#------------------generate wt.in files and Hr files------------------#

sp.run(["/home/users/shenyc/openmx/openmx3.9/work/analysis_example", "openmx.scfout"])
sp.run(['cp', '/home/users/shenyc/Myscripts/OPENMX_files/auto_openmx_wannier/openmx2wt.py', 'openmx2wt.py'])
sp.run(["conda", "run", "-n", "mp_api", "python", "openmx2wt.py"])



