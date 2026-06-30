import numpy as np
import os, glob
from ase.io import read, write
from collections import Counter
from ase.data import atomic_masses, atomic_numbers

#data_pattern = "*.vasp"
#p_path = '/home/users/nus/zhongpc/scratch/zzh/lmp-mace/MACE-matpes-r2scan-omat-ft.model-lammps.pt'
#files = sorted(glob.glob(data_pattern))
#T = 100
#L = 0.3

def input_NVT(pos, p_path, T, fs):
    ff = pos.split('.')[0]
    T_dir = f"{T:.1f}"+'_dir'
    ff2 = ff+'/NVT/'+T_dir
    os.system("mkdir " + ff+'/NVT/')
    os.system("mkdir " + ff2)

    atoms = read(pos)
    symbols = atoms.get_chemical_symbols()
    count = Counter(symbols)

    elem = []
    num = []
    mass = []
    for e in count:
        elem.append(e)
        num.append(count[e])
        Z = atomic_numbers[e]
        mass.append(atomic_masses[Z])
    N = sum(num)


    file_hes = open(ff2+'/in.lammps','w')
    file_hes.write('# I_L'+'\n'+'\n')
    file_hes.write("units           metal"+'\n')
    file_hes.write("boundary        p p p"+'\n')
    file_hes.write("atom_style      atomic"+'\n')
    file_hes.write("atom_modify     map yes"+'\n')
    file_hes.write("newton          on"+'\n'+'\n')
    file_hes.write("neighbor        2.0 bin"+'\n')
    file_hes.write("neigh_modify    every 10 delay 0 check no"+'\n'+'\n')
    file_hes.write("read_data       SUPERCELL.lmp"+'\n')
    for i in range(len(elem)):
        file_hes.write("mass            "+str(i+1)+f" {mass[i]:.3f}"+'\n')
    file_hes.write('\n'+'\n')
    file_hes.write("pair_style      mace no_domain_decomposition"+'\n')
    file_hes.write("pair_coeff      * * "+p_path+' ')
    for i in range(len(elem)):
        file_hes.write(elem[i]+" ")
    file_hes.write('\n'+'\n')
    file_hes.write(f"velocity        all create {T:.1f} 8257 dist gaussian mom yes"+'\n'+'\n')
    file_hes.write(f"fix             1 all nvt temp {T:.1f} {T:.1f} 0.1"+'\n')
    file_hes.write("timestep        0.001"+'\n')
    file_hes.write("thermo_style    custom step pe ke etotal temp press vol"+'\n')
    file_hes.write("thermo          100"+'\n')
    file_hes.write("dump            1 all custom 100 SUPERCELL.dump id type xu yu zu vx vy vz"+'\n'+'\n')
    file_hes.write(f"run             {fs}"+'\n')
    file_hes.write("write_data        final.data"+'\n')
    file_hes.close()


def input_NPT(pos, p_path, T_dir, T, fs):
    ff = pos.split('.')[0]
    T_dir = T_dir
    ff2 = ff+'/NPT/'+T_dir
    os.system("mkdir " + ff+'/NPT/')
    os.system("mkdir " + ff2)

    atoms = read(pos)
    symbols = atoms.get_chemical_symbols()
    count = Counter(symbols)

    elem = []
    num = []
    mass = []
    for e in count:
        elem.append(e)
        num.append(count[e])
        Z = atomic_numbers[e]
        mass.append(atomic_masses[Z])
    N = sum(num)


    file_hes = open(ff2+'/in.lammps','w')
    file_hes.write('# I_L'+'\n'+'\n')
    file_hes.write("units           metal"+'\n')
    file_hes.write("boundary        p p p"+'\n')
    file_hes.write("atom_style      atomic"+'\n')
    file_hes.write("atom_modify     map yes"+'\n')
    file_hes.write("newton          on"+'\n'+'\n')
    file_hes.write("neighbor        2.0 bin"+'\n')
    file_hes.write("neigh_modify    every 10 delay 0 check no"+'\n'+'\n')
    file_hes.write("read_data       SUPERCELL.lmp"+'\n')
    for i in range(len(elem)):
        file_hes.write("mass            "+str(i+1)+f" {mass[i]:.3f}"+'\n')
    file_hes.write('\n'+'\n')
    file_hes.write("pair_style      mace no_domain_decomposition"+'\n')
    file_hes.write("pair_coeff      * * "+p_path+' ')
    for i in range(len(elem)):
        file_hes.write(elem[i]+" ")
    file_hes.write('\n'+'\n')
    file_hes.write(f"velocity        all create {T:.1f} 8257 dist gaussian mom yes"+'\n'+'\n')
    file_hes.write(f"fix             npt_run all npt   temp {T:.1f} {T:.1f} 0.1   iso 1 1 1.0"+'\n')
    file_hes.write("timestep        0.001"+'\n')
    file_hes.write("thermo_style    custom step pe ke etotal temp press vol"+'\n')
    file_hes.write("thermo          100"+'\n')
    file_hes.write("dump            1 all custom 10 SUPERCELL.dump id type xu yu zu vx vy vz"+'\n'+'\n')
    file_hes.write(f"run             {fs}"+'\n')
    file_hes.write("write_data        final.data"+'\n')
    file_hes.close()
