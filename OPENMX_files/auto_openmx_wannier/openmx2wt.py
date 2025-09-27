#!/home/jhli/anaconda3/envs/deeph/bin/python
# -*- coding: utf-8 -*
"""
Created on Sun Apr  9 18:00:59 2023

@author: ljhlx
"""

import numpy as np 
import re
import sys

def parse_atom_basis_string(s):
    # 以“-”分割字符串
    parts = s.split("-")
    if len(parts) < 2:
        return None
    
    # 获取后一部分
    suffix = parts[1]
    
    # 解析每个字符和数字对
    result = []
    i = 0
    while i < len(suffix):
        # 匹配字母
        char = suffix[i]
        # 匹配数字（如果存在）
        if i + 1 < len(suffix) and suffix[i + 1].isdigit():
            count = int(suffix[i + 1])
            result.extend([char] * count)
            i += 1
        else:
            result.append(char)
        i += 1
    
    # 进一步替换
    final_result = []
    for item in result:
        if item == "p":
            final_result.extend(["px", "py", "pz"])
        elif item == "d":
            final_result.extend(["dz2", "dx2-y2", "dxy", "dxz", "dyz"])
        elif item == "f":
            final_result.extend(["fz3", "fxz2", "fyz2", "fzx2", "fxyz", "fx3", "fy3x2"])
        else:
            final_result.append(item)
    
    return final_result

class openmx_input():
    def __init__(self, openmx_file="MBT.dat"):
        self.openmx_file = openmx_file 
        self.process_openmx_file()
        
    def process_openmx_file(self): # poscar must be sorted along z axis
        file = open(self.openmx_file)
        file_content = file.readlines()
        file.close()
        
        for ind_line, line in enumerate(file_content):
            if "Atoms.UnitVectors.Unit" in line:
                latt_unit = line.split()[1]
            if "<Atoms.UnitVectors" in line:
                latt_start_line = ind_line + 1
            if "Atoms.SpeciesAndCoordinates.Unit" in line:
                atom_info_unit = line.split()[1]
            if "<Atoms.SpeciesAndCoordinates" in line:
                atom_info_start_line = ind_line + 1 
            if "Atoms.Number" in line:
                atom_num_sum = int(line.split()[1])
            if "<Definition.of.Atomic.Species" in line:
                atom_spec_basis_line = ind_line + 1
            if "Species.Number" in line:
                atom_spec_num = int(line.strip().split()[-1])
        
        latt_list = []
        for line in file_content[latt_start_line : latt_start_line + 3]:
            if latt_unit.lower() == "ang":
                latt_list.append([float(x) for x in line.split()])
            elif latt_unit.lower() =="au":
                latt_list.append([float(x) * 0.529177249 for x in line.split()])
        
        atom_info_list = []
        for line in file_content[atom_info_start_line: atom_info_start_line + atom_num_sum]:
            atom_info = line.split()[1:5]
            atom_info_list.append(atom_info)
        
        atom_info_list.sort()
        atom_spec_dict = {}
            
        for ind_atom, atom_info in enumerate(atom_info_list):
            atom_spec = atom_info[0]
            #print(atom_info)
            if atom_spec not in atom_spec_dict:
                atom_spec_dict[atom_spec] = {}
                atom_spec_dict[atom_spec]["atom_num"] = 1
            else:
                atom_spec_dict[atom_spec]["atom_num"] += 1 
            
            if atom_info_unit.lower() == "au":
                atom_info_list[ind_atom] = [atom_spec, ] + [float(x) * 0.529177249 for x in atom_info[1:]]
            else:
                atom_info_list[ind_atom] = [atom_spec, ] + [float(x) for x in atom_info[1:]]
        
        for line in file_content[atom_spec_basis_line: atom_spec_basis_line + atom_spec_num]:
            atom_spec, atom_spec_basis_str = line.strip().split()[:2]
            atom_spec_basis_spdf_list = parse_atom_basis_string(atom_spec_basis_str)
            atom_spec_dict[atom_spec]["basis_list"] = atom_spec_basis_spdf_list
        
        self.latt_list = latt_list
        self.atom_spec_dict = atom_spec_dict
        self.atom_info_unit = atom_info_unit
        self.atom_info_list = atom_info_list
    
    
    def write_poscar_file(self, poscar_file="POSCAR"):
        file_out = open(poscar_file, "w")
        file_out.write("POSCAR converted from %s\n 1\n"%self.openmx_file)
        for latt in self.latt_list:
            file_out.write("%.8f  %.8f  %.8f\n"%(latt[0], latt[1], latt[2]))
        
        atom_spec_line = " "
        atom_num_line = " "
        for atom_spec, atom_spec_info in self.atom_spec_dict.items():
            atom_spec_line += atom_spec + " "
            atom_num_line += str(atom_spec_info["atom_num"]) + " "
        
        file_out.write(atom_spec_line + "\n")
        file_out.write(atom_num_line + "\n")
        if self.atom_info_unit.lower() == "ang":
            file_out.write("Cart \n")
        
        elif self.atom_info_unit.lower() == "frac":
            file_out.write("Direct \n")
        
        elif self.atom_info_unit.lower() == "au":
            file_out.write("Cart \n")
        
        for atom_info in self.atom_info_list:
            file_out.write("%.8f   %.8f   %.8f  %s\n"%(atom_info[1], atom_info[2], atom_info[3], atom_info[0]))
            
    def write_wt_file(self):
        file_out = open("wt.in", "w")
        
        # header file
        file_out.write("""&TB_FILE
Hrfile = 'H.dat'
Overlapfile = 'S.dat'
Package = 'OPENMX'
Is_Sparse_Hr = T
Is_Sparse = F
Orthogonal_Basis = F
/

""")

        # latt and atom pos 
        file_out.write("\nLATTICE\n")
        file_out.write("Angstrom\n")
        for latt in self.latt_list:
            file_out.write("%.8f  %.8f  %.8f\n"%(latt[0], latt[1], latt[2]))
        
        file_out.write("\n")
        file_out.write("ATOM_POSITIONS\n")
        file_out.write("%d   ! number of atoms for projectors\n"%len(self.atom_info_list))
        
        if self.atom_info_unit.lower() == "ang":
            file_out.write("Cart ! Direct or Cartisen coordinate\n")
        
        elif self.atom_info_unit.lower() == "frac":
            file_out.write("Direct ! Direct or Cartisen coordinate\n")
        
        for atom_info in self.atom_info_list:
            file_out.write("%s  %10.6f  %10.6f  %10.6f\n"%(atom_info[0], \
                           atom_info[1], atom_info[2], atom_info[3]))
        file_out.write("\n")
        
        # PROJECTOR
        file_out.write("PROJECTORS\n")
        for atom_info in self.atom_info_list:
            atom_spec = atom_info[0]
            file_out.write("%d "%len(self.atom_spec_dict[atom_spec]["basis_list"]))
        file_out.write(" ! number of projectors\n")
        
        for atom_info in self.atom_info_list:
            atom_spec = atom_info[0]
            atom_basis = " ".join(self.atom_spec_dict[atom_spec]["basis_list"])
            file_out.write("%s %s \n"%(atom_spec, atom_basis))
        file_out.write("\n")
        
        file_out.write("""
                       
!> bulk band structure calculation flag
&CONTROL
BulkBand_calc         = T
FindNodes_calc        = T
/
                       
&SYSTEM
NSLAB = 5              ! for thin film system
NumOccupied = 48        ! NumOccupied
SOC = 1                 ! soc
E_FERMI = -4.4195        ! e-fermi
/

&PARAMETERS
Nk1 = 9            ! number k points  odd number would be better
Nk2 = 9            ! number k points  odd number would be better
Nk3 = 9            ! number k points  odd number would be better
Gap_threshold = 0.0001
/

SURFACE            ! (001) surface
 1  0  0
 0  1  0
 
KPATH_BULK            ! k point path
4              ! number of k line only for bulk band
G 0.00000 0.00000 0.0000 Z 0.00000 0.00000 0.5000
Z 0.00000 0.00000 0.5000 F 0.50000 0.50000 0.0000
F 0.50000 0.50000 0.0000 G 0.00000 0.00000 0.0000
G 0.00000 0.00000 0.0000 L 0.50000 0.00000 0.0000

KPOINTS_3D
4
Direct
0.0  0.0  0.0
0.5  0.0  0.0
0.0  0.5  0.0
0.0  0.0  0.5

KCUBE_BULK
0 0 0
1 0 0
0 1 0
0 0 1

KPLANE_BULK
 0.00  0.00  0.00   ! Original point for 3D k plane
 1.00  0.00  0.00   ! The first vector to define 3d k space plane
 0.00  0.50  0.00   ! The second vector to define 3d k space plane

KPATH_SLAB
2        ! numker of k line for 2D case
K 0.33 0.67 G 0.0 0.0  ! k path for 2D case
G 0.0 0.0 M 0.5 0.5

KPLANE_SLAB
-0.1 -0.1      ! Original point for 2D k plane
 0.2  0.0      ! The first vector to define 2D k plane
 0.0  0.2      ! The second vector to define 2D k plane  for arc plots
""")
        
mbt_openmx_in = openmx_input(openmx_file="openmx.dat")
mbt_openmx_in.process_openmx_file()
mbt_openmx_in.write_poscar_file(poscar_file="POSCAR")
mbt_openmx_in.write_wt_file()
