
import ase
from ase.visualize import view
from .context import flatgraphene as fg


#standard hexagonal AB bilayer with simple inputs
atoms=fg.shift.make_graphene(stacking=['A','B'],cell_type='hex',n_1=3,n_2=3,
			     lat_con=0.0,n_layer=2,sep=1.0,a_nn=1.5)
ase.visualize.view(atoms)

#arbitrary 5 layer system with per-layer inputs
atoms=fg.shift.make_graphene(stacking=['A','B','A','SP','SP'],cell_type='rect',n_1=3,n_2=3,
                             lat_con=0.0,n_layer=5,sep=[1.0,1.5,1.5,3.0,1.0],a_nn=1.5,
                             sym='C',mass=12,h_vac=6)
ase.visualize.view(atoms)
