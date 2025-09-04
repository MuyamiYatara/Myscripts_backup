#!/data/home/zy/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 10:54:07 2023

@author: ljhlx
"""

import os
import sys
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np
# OpenMX PAO basis set information
OPENMX_PAO = """
VPS 	Valence electrons 	Quick 	Standard 	Precise
H_PBE19 	1.0 	H5.0-s2 	H6.0-s2p1 	H7.0-s2p2d1
He_PBE19 	2.0 	He8.0-s1p1 	He8.0-s2p1 	He10.0-s2p2d1
Li_PBE19 	3.0 	Li8.0-s3p1 	Li8.0-s3p2 	Li8.0-s3p2d1
Be_PBE19 	2.0 	Be7.0-s2p1 	Be7.0-s2p2 	Be7.0-s3p2d1
B_PBE19 	3.0 	B7.0-s2p2 	B7.0-s2p2d1 	B7.0-s3p2d2
C_PBE19 	4.0 	C6.0-s2p2 	C6.0-s2p2d1 	C6.0-s3p2d2
N_PBE19 	5.0 	N6.0-s2p2 	N6.0-s2p2d1 	N6.0-s3p2d2
O_PBE19 	6.0 	O6.0-s2p2 	O6.0-s2p2d1 	O6.0-s3p2d2
F_PBE19 	7.0 	F6.0-s2p2 	F6.0-s2p2d1 	F6.0-s3p3d2f1
Ne_PBE19 	8.0 	Ne9.0-s2p2 	Ne9.0-s2p2d1 	Ne9.0-s3p2d2
Na_PBE19 	9.0 	Na9.0-s3p2 	Na9.0-s3p2d1 	Na9.0-s3p2d2
Mg_PBE19 	8.0 	Mg9.0-s2p2 	Mg9.0-s3p2d1 	Mg9.0-s3p2d2
Al_PBE19 	3.0 	Al7.0-s2p1d1 	Al7.0-s2p2d1 	Al7.0-s3p2d2
Si_PBE19 	4.0 	Si7.0-s2p1d1 	Si7.0-s2p2d1 	Si7.0-s3p3d2
P_PBE19 	5.0 	P7.0-s2p2d1 	P7.0-s2p2d1f1 	P7.0-s3p2d2f1
S_PBE19 	6.0 	S7.0-s2p2d1 	S7.0-s2p2d1f1 	S7.0-s3p2d2f1
Cl_PBE19 	7.0 	Cl7.0-s2p2d1 	Cl7.0-s2p2d1f1 	Cl7.0-s3p2d2f1
Ar_PBE19 	8.0 	Ar9.0-s2p2d1 	Ar9.0-s2p2d1f1 	Ar9.0-s3p2d2f1
K_PBE19 	9.0 	K10.0-s3p2 	K10.0-s3p2d1 	K10.0-s3p2d2
Ca_PBE19 	10.0 	Ca9.0-s3p2 	Ca9.0-s3p2d1 	Ca9.0-s3p2d2
Sc_PBE19 	11.0 	Sc9.0-s2p2d1 	Sc9.0-s3p2d1 	Sc9.0-s3p2d2
Ti_PBE19 	12.0 	Ti7.0-s2p2d1 	Ti7.0-s3p2d1 	Ti7.0-s3p2d2f1
V_PBE19 	13.0 	V6.0-s2p2d1 	V6.0-s3p2d1 	V6.0-s3p2d2f1
Cr_PBE19 	14.0 	Cr6.0-s2p2d1 	Cr6.0-s3p2d1 	Cr6.0-s3p2d2f1
Mn_PBE19 	15.0 	Mn6.0-s2p2d1 	Mn6.0-s3p2d1 	Mn6.0-s3p2d2f1
Fe_PBE19H 	16.0 	Fe5.5H-s2p2d1 	Fe5.5H-s3p2d1 	Fe5.5H-s3p2d2f1
Fe_PBE19S 	14.0 	Fe6.0S-s2p2d1 	Fe6.0S-s3p2d1 	Fe6.0S-s3p2d2f1
Co_PBE19H 	17.0 	Co6.0H-s2p2d1 	Co6.0H-s3p2d1 	Co6.0H-s3p2d2f1
Co_PBE19S 	15.0 	Co6.0S-s2p2d1 	Co6.0S-s3p2d1 	Co6.0S-s3p2d2f1
Ni_PBE19H 	18.0 	Ni6.0H-s2p2d1 	Ni6.0H-s3p2d1 	Ni6.0H-s3p2d2f1
Ni_PBE19S 	16.0 	Ni6.0S-s2p2d1 	Ni6.0S-s3p2d1 	Ni6.0S-s3p2d2f1
Cu_PBE19H 	19.0 	Cu6.0H-s2p2d1 	Cu6.0H-s3p2d1 	Cu6.0H-s3p2d2f1
Cu_PBE19S 	11.0 	Cu6.0S-s2p1d1 	Cu6.0S-s3p2d1 	Cu6.0S-s3p2d2f1
Zn_PBE19H 	20.0 	Zn6.0H-s2p2d1 	Zn6.0H-s3p2d1 	Zn6.0H-s3p2d2f1
Zn_PBE19S 	12.0 	Zn6.0S-s2p1d1 	Zn6.0S-s3p2d1 	Zn6.0S-s3p2d2f1
Ga_PBE19 	13.0 	Ga7.0-s2p2d1 	Ga7.0-s3p2d2 	Ga7.0-s3p2d2f1
Ge_PBE19 	4.0 	Ge7.0-s2p1d1 	Ge7.0-s3p2d2 	Ge7.0-s3p2d2f1
As_PBE19 	15.0 	As7.0-s3p2d1 	As7.0-s3p2d2 	As7.0-s3p2d2f1
Se_PBE19 	6.0 	Se7.0-s3p2d1 	Se7.0-s3p2d2 	Se7.0-s3p2d2f1
Br_PBE19 	7.0 	Br7.0-s3p2d1 	Br7.0-s3p2d2 	Br7.0-s3p2d2f1
Kr_PBE19 	8.0 	Kr10.0-s2p2d1 	Kr10.0-s3p2d2 	Kr10.0-s3p2d2f1
Rb_PBE19 	9.0 	Rb11.0-s2p2d1 	Rb11.0-s3p2d2 	Rb11.0-s3p2d2f1
Sr_PBE19 	10.0 	Sr10.0-s2p2d1 	Sr10.0-s3p2d2 	Sr10.0-s3p3d2f1
Y_PBE19 	11.0 	Y10.0-s3p2d1 	Y10.0-s3p2d2 	Y10.0-s3p3d2f1
Zr_PBE19 	12.0 	Zr7.0-s3p2d1 	Zr7.0-s3p2d2 	Zr7.0-s3p2d2f1
Nb_PBE19 	13.0 	Nb7.0-s3p2d1 	Nb7.0-s3p2d2 	Nb7.0-s3p2d2f1
Mo_PBE19 	14.0 	Mo7.0-s3p2d1 	Mo7.0-s3p2d2 	Mo7.0-s3p2d2f1
Tc_PBE19 	15.0 	Tc7.0-s3p2d1 	Tc7.0-s3p2d2 	Tc7.0-s3p2d2f1
Ru_PBE19 	14.0 	Ru7.0-s3p2d1 	Ru7.0-s3p2d2 	Ru7.0-s3p2d2f1
Rh_PBE19 	15.0 	Rh7.0-s3p2d1 	Rh7.0-s3p2d2 	Rh7.0-s3p2d2f1
Pd_PBE19 	16.0 	Pd7.0-s3p2d1 	Pd7.0-s3p2d2 	Pd7.0-s3p2d2f1
Ag_PBE19 	17.0 	Ag7.0-s3p2d1 	Ag7.0-s3p2d2 	Ag7.0-s3p2d2f1
Cd_PBE19 	12.0 	Cd7.0-s3p2d1 	Cd7.0-s3p2d2 	Cd7.0-s3p2d2f1
In_PBE19 	13.0 	In7.0-s3p2d1 	In7.0-s3p2d2 	In7.0-s3p2d2f1
Sn_PBE19 	14.0 	Sn7.0-s3p2d1 	Sn7.0-s3p2d2 	Sn7.0-s3p2d2f1
Sb_PBE19 	15.0 	Sb7.0-s3p2d1 	Sb7.0-s3p2d2 	Sb7.0-s3p2d2f1
Te_PBE19 	16.0 	Te7.0-s3p2d2 	Te7.0-s3p2d2f1 	Te7.0-s3p3d2f1
I_PBE19 	7.0 	I7.0-s3p2d2 	I7.0-s3p2d2f1 	I7.0-s3p3d2f1
Xe_PBE19 	8.0 	Xe11.0-s3p2d1 	Xe11.0-s3p2d2 	Xe11.0-s3p2d2f1
Cs_PBE19 	9.0 	Cs12.0-s3p2d1 	Cs12.0-s3p2d2 	Cs12.0-s3p2d2f1
Ba_PBE19 	10.0 	Ba10.0-s3p2d1 	Ba10.0-s3p2d2 	Ba10.0-s3p2d2f1
La_PBE19 	11.0 	La8.0-s3p2d1f1 	La8.0-s3p2d2f1 	La8.0-s3p3d2f1
Ce_PBE19 	12.0 	Ce8.0-s3p2d1f1 	Ce8.0-s3p2d2f1 	Ce8.0-s3p3d2f1
Pr_PBE19 	13.0 	Pr8.0-s3p2d1f1 	Pr8.0-s3p2d2f1 	Pr8.0-s3p3d2f1
Nd_PBE19 	14.0 	Nd8.0-s3p2d1f1 	Nd8.0-s3p2d2f1 	Nd8.0-s3p3d2f1
Pm_PBE19 	15.0 	Pm8.0-s3p2d1f1 	Pm8.0-s3p2d2f1 	Pm8.0-s3p3d2f1
Sm_PBE19 	16.0 	Sm8.0-s3p2d1f1 	Sm8.0-s3p2d2f1 	Sm8.0-s3p3d2f1
Dy_PBE19 	20.0 	Dy8.0-s3p2d1f1 	Dy8.0-s3p2d2f1 	Dy8.0-s3p3d2f1
Ho_PBE19 	21.0 	Ho8.0-s3p2d1f1 	Ho8.0-s3p2d2f1 	Ho8.0-s3p3d2f1
Lu_PBE19 	11.0 	Lu8.0-s3p2d2 	Lu8.0-s3p2d2f1 	Lu8.0-s3p3d2f1
Hf_PBE19 	12.0 	Hf9.0-s3p2d2 	Hf9.0-s3p2d2f1 	Hf9.0-s3p3d2f1
Ta_PBE19 	13.0 	Ta7.0-s3p2d2 	Ta7.0-s3p2d2f1 	Ta7.0-s3p3d2f1
W_PBE19 	12.0 	W7.0-s3p2d2 	W7.0-s3p2d2f1 	W7.0-s3p3d2f1
Re_PBE19 	15.0 	Re7.0-s3p2d2 	Re7.0-s3p2d2f1 	Re7.0-s3p3d2f1
Os_PBE19 	14.0 	Os7.0-s3p2d2 	Os7.0-s3p2d2f1 	Os7.0-s3p3d2f1
Ir_PBE19 	15.0 	Ir7.0-s3p2d2 	Ir7.0-s3p2d2f1 	Ir7.0-s3p3d2f1
Pt_PBE19 	16.0 	Pt7.0-s3p2d2 	Pt7.0-s3p2d2f1 	Pt7.0-s3p3d2f1
Au_PBE19 	17.0 	Au7.0-s3p2d2 	Au7.0-s3p2d2f1 	Au7.0-s3p3d2f1
Hg_PBE19 	18.0 	Hg8.0-s3p2d2 	Hg8.0-s3p2d2f1 	Hg8.0-s3p3d2f1
Tl_PBE19 	19.0 	Tl8.0-s3p2d2 	Tl8.0-s3p2d2f1 	Tl8.0-s3p3d2f1
Pb_PBE19 	14.0 	Pb8.0-s3p2d2 	Pb8.0-s3p2d2f1 	Pb8.0-s3p3d2f1
Bi_PBE19 	15.0 	Bi8.0-s3p2d2 	Bi8.0-s3p2d2f1 	Bi8.0-s3p3d2f1
"""

def parse_openmx_pao_dict():
    """Parse OpenMX PAO basis set data into a dictionary."""
    openmx_pao_file = OPENMX_PAO.split("\n")
    openmx_pao_dict = {}

    for line in openmx_pao_file[2:]:
        if not line.strip():
            continue
        tmp = line.split()
        element_string = tmp[0].split("_")[0]
        vps_name = tmp[0]
        num_ele = float(tmp[1])
        quick_basis = tmp[2]
        standard_basis = tmp[3]
        precise_basis = tmp[4]

        # Store the information in the dictionary
        if element_string not in openmx_pao_dict:
            openmx_pao_dict[element_string] = {
                "num_ele": num_ele,
                "quick_basis": quick_basis,
                "standard_basis": standard_basis,
                "precise_basis": precise_basis,
                "vps_name" : vps_name
            }

    return openmx_pao_dict

def calculate_kmesh(lattice, dq_grid=0.05):
        """根据晶格参数计算k网格"""
        # 计算晶格矩阵
        avec = np.array(lattice)
        
        # 计算倒易晶格矩阵
        bvec = 2 * np.pi * np.linalg.inv(avec).T
        
        # 计算原胞体积
        v_uc = abs(np.linalg.det(avec))
        
        # 计算晶格常数的范数（每个晶格矢量的长度）
        norm = np.zeros(3)
        for i in range(3):
            norm[i] = np.sqrt(np.dot(avec[i], avec[i]))
        
        # 计算常数 c
        c = (norm[0] * norm[1] * norm[2] / v_uc)**(1/3.0)
        
        # 计算每个方向的 k-点数
        nq = np.zeros(3, dtype=int)
        for i in range(3):
            nq[i] = round(c / norm[i] / dq_grid)
            if nq[i] == 0:
                nq[i] = 1  # 确保至少有一个 k-点
            elif nq[i] > 24:
                nq[i] = 24  # 限制最大k点数，避免计算过重
        
        return nq[0], nq[1], nq[2]


openmx_pao_dict = parse_openmx_pao_dict()

class OpenMXInput:
    """Class to handle conversion of POSCAR to OpenMX input format."""
    
    def __init__(self, poscar_file="POSCAR"):
        """Initialize with POSCAR file and extract element list."""
        self.structure = Structure.from_file(poscar_file)
        self.element_list = [x.symbol for x in self.structure.composition.elements]

    def append_latt_info(self, in_file="openmx.dat_ori", out_file="openmx.dat", basis_set="standard_basis"):
        """Append lattice and atomic species information to OpenMX input file."""
        
        # Validate the basis set
        basis_set_mapping = {
            "s": "standard_basis",
            "q": "quick_basis",
            "p": "precise_basis"
        }
        
        basis_set = basis_set_mapping.get(basis_set[0].lower())
        if not basis_set:
            raise ValueError("Basis set must start with 's', 'q', or 'p'")

        # Copy input file to output file
        cp_command = "cp" if os.name == "posix" else "copy"
        os.system(f"{cp_command} {in_file} {out_file}")

        # Open output file to append information
        with open(out_file, "a") as openmx_in_file:

            # KPOINTS 
            openmx_in_file.write("\nscf.Kgrid                     ")
            k1, k2, k3 = calculate_kmesh(self.structure.lattice.matrix, 0.03)
            openmx_in_file.write(f"{k1} {k2} {k3}\n")
            
            # Crystal lattice
            openmx_in_file.write("\nAtoms.UnitVectors.Unit  Ang\n<Atoms.UnitVectors\n")
            for latt in self.structure.lattice.matrix:
                openmx_in_file.write(f"   {latt[0]:16.10f}   {latt[1]:16.10f}   {latt[2]:16.10f}\n")
            openmx_in_file.write("Atoms.UnitVectors>\n\n")
            
            # Atom species and definitions
            openmx_in_file.write(f"Species.Number       {len(self.element_list)}\n")
            openmx_in_file.write("<Definition.of.Atomic.Species\n")
            for element in self.element_list:
                vps = openmx_pao_dict[element]["vps_name"]
                basis_type = openmx_pao_dict[element][basis_set]
                openmx_in_file.write(f"{element}   {basis_type}   {vps}\n")
            openmx_in_file.write("Definition.of.Atomic.Species>\n\n")
            
            # Atomic species and coordinates
            openmx_in_file.write(f"Atoms.Number        {len(self.structure)}\n")
            openmx_in_file.write("Atoms.SpeciesAndCoordinates.Unit   Ang # Ang|AU\n<Atoms.SpeciesAndCoordinates\n")
            for ind_atom, atom in enumerate(self.structure):
                atom_spec = atom.species_string
                atom_line = f"{ind_atom + 1:>6d} {atom_spec:>6s} "
                atom_cart_coords = atom.coords
                atom_line += f"{atom_cart_coords[0]:17.10f} {atom_cart_coords[1]:17.10f} {atom_cart_coords[2]:17.10f} "
                val_ele_num = openmx_pao_dict[atom_spec]["num_ele"]
                atom_line += f"{val_ele_num/2:8.4f} {val_ele_num/2:8.4f} 0.0 0.0\n"
                openmx_in_file.write(atom_line)
            openmx_in_file.write("Atoms.SpeciesAndCoordinates>\n")

# Main execution
if len(sys.argv) != 5:
    print("Usage: python poscar2openmx.py POSCAR standard_basis openmx.dat_ori openmx.dat")
    sys.exit(1)

poscar_file = sys.argv[1]
basis_set = sys.argv[2]
in_file = sys.argv[3]
out_file = sys.argv[4]

mbt_film = OpenMXInput(poscar_file=poscar_file)
mbt_film.append_latt_info(in_file=in_file, basis_set=basis_set, out_file=out_file)

# Usage example: python poscar2openmx.py POSCAR standard_basis openmx.dat_ori openmx.dat
