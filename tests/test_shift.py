
import ase
from ase.visualize import view
from .context import flatgraphene as fg

#note the inputs are all given with their names for clarity, but this is not necessary
#the nearest neighbor distance (in-plane) a_nn is optional, and overrides the lat_con variable
#  meaning the value of lat_con is unused
atoms=fg.shift.make_graphene(alignment='AB',cell_type='rect',n_layer=3,
                                      n_1=3,n_2=3,lat_con=0.0,a_nn=1.5,sep=1.0,
                                      sym='C',mass=12)
ase.visualize.view(atoms)
