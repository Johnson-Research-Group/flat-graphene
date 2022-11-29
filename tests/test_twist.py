
import ase
from ase.visualize import view
import flatgraphene as fg

#test number of atoms in a cell
n_21_78 = 28 #number of atoms in hexagonal cell of 21.78 degree twist angle
p_found, q_found, theta_comp = fg.twist.find_p_q(21.79)
atoms = fg.twist.make_graphene(cell_type='hex',n_layer=2,
                               p=p_found,q=q_found,lat_con=0.0,a_nn=1.5,
                               sep=3.35,h_vac=3)
n_fg = (atoms.get_positions()).shape[0] #number of atoms in returned structure
if (n_fg != n_21_78):
    print('ERROR: 21.78 degree hex cell has wrong number of atoms')

#small bilayer system with simple inputs
p_found, q_found, theta_comp = fg.twist.find_p_q(21.79)
atoms = fg.twist.make_graphene(cell_type='hex',n_layer=2,
                               p=p_found,q=q_found,lat_con=0.0,a_nn=1.5,
                               sep=3.35,h_vac=3)
ase.visualize.view(atoms)

#large, arbitrary 5 layer system with per-layer inputs
p_found, q_found, theta_comp = fg.twist.find_p_q(9.43)
atoms = fg.twist.make_graphene(cell_type='hex',n_layer=5,
                               p=p_found,q=q_found,lat_con=0.0,a_nn=1.5,
                               sep=[1.0,1.5,1.5,3.0,1.0],sym=['C','Kr','C','Kr','C'],
                               mass=[1,2,3,4,5],mol_id=[1,101,2,202,3],h_vac=3)
ase.visualize.view(atoms)
