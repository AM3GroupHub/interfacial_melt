import numpy as np
import os, glob
from ase.io import read, write
from ase.geometry import cell_to_cellpar
#from pymatgen.io.ase import AseAtomsAdaptor
#from pymatgen.io.lammps.data import LammpsData

#data_pattern = "*.vasp"
#files = sorted(glob.glob(data_pattern))

def read_SC(pos,cif):
    ff = pos.split('.')[0]
    atoms = read(cif)
    cellpar = cell_to_cellpar(atoms.cell)
    a, b, c, alpha, beta, gamma = cellpar
    comment = f"# CELL(abcABC): {a:.6f} {b:.6f} {c:.6f} {alpha:.3f} {beta:.3f} {gamma:.3f}"+" positions{angstrom} cell{angstrom}"
    os.system('mkdir '+ff)
    write(ff+"/SUPERCELL.xyz", atoms, format="xyz", comment=comment)
    write(ff+"/SUPERCELL.lmp-ase", atoms, format="lammps-data")
    write(ff+"/POSCAR", atoms)

    file3 = open(ff+"/SUPERCELL.lmp-ase",'r')
    s = file3.readlines()
    for i in range(len(s)):
        s[i] = str(s[i]).replace('\n','')
    file3.close()

    file3 = open(ff+"/SUPERCELL.lmp",'w')
    for i in range(len(s)):
        if s[i] == 'Atoms # atomic':
            atomic = i
    for i in range(atomic+2):
        file3.write(s[i]+'\n')

    Natom = int(s[2].split()[0])
    aya = 1
    file3.write(s[atomic+2].split()[0]+'  '+str(aya)+'  '+s[atomic+2].split()[2]+'  '+s[atomic+2].split()[3]+'  '+s[atomic+2].split()[4]+'  '+'\n')
    for i in range(Natom-1):
        if int(s[atomic+2+i].split()[1]) != int(s[atomic+2+i+1].split()[1]):
            aya = aya +1
        file3.write(s[atomic+2+i+1].split()[0]+'  '+str(aya)+'  '+s[atomic+2+i+1].split()[2]+'  '+s[atomic+2+i+1].split()[3]+'  '+s[atomic+2+i+1].split()[4]+'  '+'\n')
    file3.close()



def SC(pos):
    ff = pos.split('.')[0]
    os.system('mkdir SC && cp ' +pos +' SC/POSCAR')
    os.system('cd SC && phonopy --symmetry --tolerance=0.1')
    os.system('cp SC/BPOSCAR SC/POSCAR')

    init_pos = read("SC/POSCAR")
    N_init = len(init_pos)
    if N_init < 5:
        os.system("cd SC && phonopy -d --dim='5 5 5'")
    elif N_init >= 5 and N_init < 8:
        os.system("cd SC && phonopy -d --dim='4 4 4'")
    elif N_init >= 8 and N_init < 19:
        os.system("cd SC && phonopy -d --dim='3 3 3'")
    elif N_init >= 19 and N_init < 81:
        os.system("cd SC && phonopy -d --dim='2 2 2'")
    elif N_init >= 81:
        os.system("cp SC/POSCAR SC/SPOSCAR")

    os.system('mkdir '+ff+' && cp SC/SPOSCAR '+ff+'/POSCAR && rm -rf SC')

    atoms = read(ff+"/POSCAR")
    cellpar = cell_to_cellpar(atoms.cell)
    a, b, c, alpha, beta, gamma = cellpar
    comment = f"# CELL(abcABC): {a:.6f} {b:.6f} {c:.6f} {alpha:.3f} {beta:.3f} {gamma:.3f}"+" positions{angstrom} cell{angstrom}"
    write(ff+"/SUPERCELL.xyz", atoms, format="xyz", comment=comment)
    write(ff+"/SUPERCELL.lmp-ase", atoms, format="lammps-data")
#    adaptor = AseAtomsAdaptor()
#    structure = adaptor.get_structure(atoms)
#    data = LammpsData.from_structure(structure, atom_style="atomic")
#    data.write_file(ff+"/SUPERCELL.lmp-bak")

    file3 = open(ff+"/SUPERCELL.lmp-ase",'r')
    s = file3.readlines()
    for i in range(len(s)):
        s[i] = str(s[i]).replace('\n','')
    file3.close()

    file3 = open(ff+"/SUPERCELL.lmp",'w')
    for i in range(len(s)):
        if s[i] == 'Atoms # atomic':
            atomic = i
    for i in range(atomic+2):
        file3.write(s[i]+'\n')
    
    Natom = int(s[2].split()[0])
    aya = 1
    file3.write(s[atomic+2].split()[0]+'  '+str(aya)+'  '+s[atomic+2].split()[2]+'  '+s[atomic+2].split()[3]+'  '+s[atomic+2].split()[4]+'  '+'\n')
    for i in range(Natom-1):
        if int(s[atomic+2+i].split()[1]) != int(s[atomic+2+i+1].split()[1]):
            aya = aya +1
        file3.write(s[atomic+2+i+1].split()[0]+'  '+str(aya)+'  '+s[atomic+2+i+1].split()[2]+'  '+s[atomic+2+i+1].split()[3]+'  '+s[atomic+2+i+1].split()[4]+'  '+'\n')
    file3.close()



def SC_d(pos,X,Y,Z):
    ff = pos.split('.')[0]
    os.system('mkdir SC && cp ' +pos +' SC/POSCAR')
    os.system('cd SC && phonopy --symmetry --tolerance=0.1')
    os.system('cp SC/BPOSCAR SC/POSCAR')

    init_pos = read("SC/POSCAR")
    N_init = len(init_pos)
    os.system(f"cd SC && phonopy -d --dim='{X} {Y} {Z}'")

    os.system('mkdir '+ff+' && cp SC/SPOSCAR '+ff+'/POSCAR && rm -rf SC')

    atoms = read(ff+"/POSCAR")
    cellpar = cell_to_cellpar(atoms.cell)
    a, b, c, alpha, beta, gamma = cellpar
    comment = f"# CELL(abcABC): {a:.6f} {b:.6f} {c:.6f} {alpha:.3f} {beta:.3f} {gamma:.3f}"+" positions{angstrom} cell{angstrom}"
    write(ff+"/SUPERCELL.xyz", atoms, format="xyz", comment=comment)
#    write(ff+"/SUPERCELL.lmp", atoms, format="lammps-data")

    write(ff+"/SUPERCELL.lmp-ase", atoms, format="lammps-data")
#    adaptor = AseAtomsAdaptor()
#    structure = adaptor.get_structure(atoms)
#    data = LammpsData.from_structure(structure, atom_style="atomic")
#    data.write_file(ff+"/SUPERCELL.lmp-bak")

    file3 = open(ff+"/SUPERCELL.lmp-ase",'r')
    s = file3.readlines()
    for i in range(len(s)):
        s[i] = str(s[i]).replace('\n','')
    file3.close()

    file3 = open(ff+"/SUPERCELL.lmp",'w')
    for i in range(len(s)):
        if s[i] == 'Atoms # atomic':
            atomic = i
    for i in range(atomic+2):
        file3.write(s[i]+'\n')

    Natom = int(s[2].split()[0])
    aya = 1
    file3.write(s[atomic+2].split()[0]+'  '+str(aya)+'  '+s[atomic+2].split()[2]+'  '+s[atomic+2].split()[3]+'  '+s[atomic+2].split()[4]+'  '+'\n')
    for i in range(Natom-1):
        if int(s[atomic+2+i].split()[1]) != int(s[atomic+2+i+1].split()[1]):
            aya = aya +1
        file3.write(s[atomic+2+i+1].split()[0]+'  '+str(aya)+'  '+s[atomic+2+i+1].split()[2]+'  '+s[atomic+2+i+1].split()[3]+'  '+s[atomic+2+i+1].split()[4]+'  '+'\n')
    file3.close()
