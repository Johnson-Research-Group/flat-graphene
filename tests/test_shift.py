
import ase
from ase.visualize import view
from .context import flatgraphene as fg

#standard AB bilayer with simple inputs
atoms=fg.shift.make_graphene(alignment='AB',cell_type='rect',n_layer=2,
                                      n_1=3,n_2=3,lat_con=0.0,a_nn=1.5,sep=1.0)
ase.visualize.view(atoms)


#arbitrary 5 layer system with per-layer inputs
atoms=fg.shift.make_graphene(alignment=['AB','AA','SP','SP'],cell_type='rect',n_layer=5,
                                      n_1=3,n_2=3,lat_con=0.0,a_nn=1.5,sep=[1.0,1.5,1.5,3.0],
                                      sym='C',mass=12)
ase.visualize.view(atoms)
