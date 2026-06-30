import numpy as np
import os, glob

import gen_sc, gen_input

data_pattern = "*.vasp"
p_path = '/scratch/users/nus/zhongpc/zzh/entropy/phase_diagram/Fe-B-final/mace-mh-1_finetuned_FeB_mixed.model-lammps.pt'
run_lmp = 'mpirun /home/users/nus/zhongpc/softwares/mace-lmp-ipi/lammps/build-ampere/lmp -k on g 1 -sf kk -pk kokkos newton on neigh half -in in.lammps > log'
T1 = []
Time1 = [10000]

T2 = [1800,1500,1200,900]
Time2 = [150000]

T3 = [2800]
Time3 = []

files = sorted(glob.glob(data_pattern))

for pos in files:
    ff = pos.split('.')[0]
#####################################supercell############################
#    cif = pos.split('.')[0]+'.cif'
    gen_sc.SC_d(pos,2,5,2)
#####################################Melting##############################
    for i in T1:
        for j in Time1:
            L_T_dir = f"{i:.1f}"+'_dir'
            ff2 = ff+'/NVT/'+L_T_dir
            gen_input.input_NVT(pos, p_path, i, j)
            os.system("cp "+ff+'/SUPERCELL.lmp ' +ff2)
            os.system('cd '+ff2+' && '+run_lmp)
            os.system("cp "+ff2+'/final.data ' +ff)
#####################################Quenching############################
    for i in T2:
        for j in Time2:
            L_T_dir = f"{i:.1f}"+'_dir'
            ff2 = ff+'/NPT/'+L_T_dir
            gen_input.input_NPT(pos, p_path, L_T_dir, i, j)
            os.system("cp "+ff+'/final.data ' +ff2+'/SUPERCELL.lmp ')
            os.system('cd '+ff2+' && '+run_lmp)
            os.system("cp "+ff2+'/final.data ' +ff)
#####################################Relaxtion############################
    for i in T3:
        for j in Time3:
            L_T_dir = f"{i:.1f}"+'_R_dir'
            ff2 = ff+'/NPT/'+L_T_dir
            gen_input.input_NPT(pos, p_path, L_T_dir, i, j)
            os.system("cp "+ff+'/final.data ' +ff2+'/SUPERCELL.lmp ')
            os.system('cd '+ff2+' && '+run_lmp)
            os.system("cp "+ff2+'/final.data ' +ff)
