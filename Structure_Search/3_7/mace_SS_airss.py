import numpy as np
import os
import datetime
import random
from ase.io import read, write
from ase.optimize import BFGS
from mace.calculators import MACECalculator
from ase.constraints import ExpCellFilter
from ase import units
#from ase.spacegroup import get_spacegroup
from collections import Counter
import spglib

N = 1000
infile = 'BFe.cell'
calc = MACECalculator(model_paths="/nfs/home/svu/zihan.zhang/software/mace_model/MACE-matpes-r2scan-omat-ft.model", device="cuda") 

seed = random.randint(0, 100000)
start_t = datetime.datetime.now()
path = os.getcwd()
for i in range(N):

    pos = str(seed)+'_'+str(i)

#    os.system('buildcell < '+ infile +' > ' + pos +'.cell')
    os.system('timeout 60s buildcell < '+ infile +' > ' + pos +'.cell')
    if os.path.exists(pos +'.cell'):
        os.system('cabal cell cif < ' + pos + '.cell > ' + pos + '.cif')
        os.system('rm ' + pos + '.cell')

        atoms = read(pos + '.cif')
        atoms.set_pbc([True, True, True])
        atoms.set_calculator(calc)
        ecf = ExpCellFilter(atoms,mask=[1, 1, 1, 0, 0, 0])
        dyn = BFGS(ecf, trajectory=pos + ".traj", logfile=pos + ".log", maxstep=0.05)
        dyn.run(fmax=0.05, steps=1000)

        if dyn.nsteps < 999:
            write(pos + "_optimized.cif", atoms) 
            os.system('rm ' + pos + '.traj')

            fileRes = open(pos+'.res','w')
 
            stress = atoms.get_stress(voigt=True)
            p = -np.mean(stress[:3]) * 160.21766208
            v = atoms.get_volume()
            e = atoms.get_potential_energy()
            TOT_Atoms = len(atoms)
#            SYS = get_spacegroup(atoms, symprec=1e-2).symbol
            cell = (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers())
            try:
                dataset = spglib.get_symmetry_dataset(cell, symprec=1e-2)
                sg = dataset["international"]
            except Exception:
                sg = "Unknown"
            SYS = sg

            a, b, c, alpha, beta, gamma = atoms.cell.cellpar()
            symbols = atoms.get_chemical_symbols()
            counts = Counter(symbols)
            n_elements = len(counts)

            fileRes.write('TITL '+pos+'  '+str(p)+'  '+str(v)+'  '+str(e)+' 0 0 '+str(TOT_Atoms)+' ('+SYS+')  n - 1'+'\n')
            fileRes.write('REM'+'\n')
            fileRes.write('REM Run started in '+ path +'\n')
            fileRes.write('REM'+'\n')
            fileRes.write('CELL 1.54180    '+str(round(a,6))+'    '+str(round(b,6))+'    '+str(round(c,6))+'    '+str(round(alpha,6))+'    '+str(round(beta,6))+'    '+str(round(gamma,6))+'    '+'\n')
            fileRes.write('LATT -1'+'\n')
            fileRes.write('SFAC ')
            for elem, num in counts.items():
                fileRes.write(elem+'  ')
#            for j in range(len(Elements)):
#                fileRes.write(str(Elements[j])+'  ')
            fileRes.write('\n')
            for j, atom in enumerate(atoms):
                x, y, z = atoms.get_scaled_positions()[j]
                fileRes.write(str(atom.symbol)+'  ')
                fileRes.write(str(1)+'  ')
                fileRes.write(str(x)+'  ')
                fileRes.write(str(y)+'  ')
                fileRes.write(str(z)+'  ')
                fileRes.write('1.0  ')
                fileRes.write('\n')
  
            fileRes.write('END'+'\n')
            fileRes.close()


end_t = datetime.datetime.now()
elapsed_sec = (end_t - start_t).total_seconds()
