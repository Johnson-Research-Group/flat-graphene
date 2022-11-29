
import numpy as np
import ase
from ase.visualize import view
import flatgraphene as fg

#single unit cell (hex/triclinic) with distance test
test_dist = 1.44
atoms = fg.shift.make_graphene(stacking='A',cell_type='hex',n_1=1,n_2=1,
                               lat_con=0.0,n_layer=1,sep=1.0,a_nn=test_dist)
pos = atoms.get_positions()
dist = np.linalg.norm(pos[0] - pos[1])
if (not np.isclose(dist,test_dist)):
    print('ERROR: single unit cell test failed')


#standard hexagonal AB bilayer with simple inputs
atoms = fg.shift.make_graphene(stacking=['A','B'],cell_type='hex',n_1=3,n_2=3,
	    		       lat_con=0.0,n_layer=2,sep=1.0,a_nn=1.5)
ase.visualize.view(atoms)

#arbitrary 5 layer system with per-layer inputs, and using every feature
atoms = fg.shift.make_graphene(stacking=['A','B','A','SP','SP'],cell_type='rect',n_1=3,n_2=3,
                               lat_con=0.0,n_layer=5,sep=[1.0,1.5,1.5,3.0,1.0],a_nn=1.5,
                               sym=['C','B','N','O','Kr'],mass=[12.01,11,10,9,8],
                               mol_id=[1,2,3,4,5],h_vac=6)
ase.visualize.view(atoms)
