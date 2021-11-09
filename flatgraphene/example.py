# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 13:15:05 2021

@author: danpa
"""
from tqdm import tqdm
import ase.io
import os
import numpy as np
import sys
import flatgraphene as fg


if __name__=="__main__":
    npoints=2
    guess_t_=np.linspace(10,30,npoints)
    #real_t=np.zeros(npoints)
    a=2.529
    sep=3.35
    a_nn=a/np.sqrt(3)
    base_path=os.getcwd()
    for t in tqdm(guess_t_):
        p_found, q_found, theta_comp = fg.twist.find_p_q(t,a_tol=4)
        test_sheet = fg.twist.make_graphene(cell_type='hex',n_layer=2,
                                   p=p_found,q=q_found,lat_con=0.0,a_nn=a_nn,sym=["B","Ti"],
                                   mass=[12.01,12.02],sep=sep,h_vac=3)
        
        real_t=theta_comp
        name=str(np.round(real_t,decimals=2)).split(".")
        name="_".join(name)
        folder=base_path+"io_t"+name
        if not os.path.exists(folder):
            os.mkdir(folder)
        
        
        #write_lammps_data(folder+"/lammps_data"+name+".data",test_sheet)
        ase.io.write(folder+"/lammps_data"+name+".data",test_sheet,format="lammps-data",atom_style="full")
    
    folder=base_path+"io_AB"
    atoms=fg.shift.make_graphene(stacking=['A','B'],cell_type='hex',n_layer=2,
		        n_1=5,n_2=5,lat_con=0.0,a_nn=a_nn,sep=sep,sym=['B','Ti'],h_vac=3)
    ase.io.write(folder+"/lammps_dataAB.data",atoms,format="lammps-data",atom_style="full")
    
    folder=base_path+"io_AA"
    atoms=fg.shift.make_graphene(stacking=['A','A'],cell_type='hex',n_layer=2,
		        n_1=5,n_2=5,lat_con=0.0,a_nn=a_nn,sep=sep,sym=['B','Ti'],h_vac=3)
    ase.io.write(folder+"/lammps_dataAA.data",atoms,format="lammps-data",atom_style="full")