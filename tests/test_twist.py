
import ase
from ase.visualize import view
from .context import flatgraphene as fg

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
                               sep=[1.0,1.5,1.5,3.0,1.0],mass=12.02,h_vac=3)
ase.visualize.view(atoms)
