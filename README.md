
## flat-graphene

Python module to generate flat (uncorrugated) graphene configurations with one or more layers, as well as twisted bilayer graphene.

### Documentation

Non-twisted graphene may be created using the `make_graphene` function, which returns an ASE atoms object. Parameters include cell type (rectangular or hexagonal), alignment (AA,AB,SP), number of unit cells (in-plane), interlayer spacing, and lattice constant.

```python
def make_graphene(alignment,cell_type,n_layer,n_1,n_2,lat_con,a_nn=None,sep=None):
    """
    Creates and returns ASE atoms object with specified graphene's geometry
    ---Input---
    alignment: specification of alignment for layers above first,
               ('AA','AB','SP') relative to first layer, list of strings
               or a single string
               *NOTE*: single string inputs result in the input string
                       alternated with 'AA'
    cell_type: unit cell type, 'rect' or 'hex', string
    n_layer: number of graphene layers (1 or 2), integer
    n_1: number of unit cells in x direction, integer
    n_2: number of unit cells in y direction, integer
    lat_con: in-plane lattice constant, float [Angstroms]
    a_nn: optional argument to specify distance between
          nearest neighbors and override lat_con, float [Angstroms]
    sep: interlayer separation for n_layer>1, n_layer-1 list of separations
         (relative to layer below) or float (uniform separations)
    ---Output---
    atoms: graphene stack, ASE atoms object
    """
```

Twisted graphene has not been officially implemented yet, though functions for it exist in `twist_bilayer_graphene.py`.

### Examples

This example creates an AB trilayer graphene system. Here AB means the odd layers (including bottom) have no in-plane shift relative to the bottom layer while even layers are shifted by a single nearest neighbor distance.
```python
import ase
from ase.visualize import view
#note the inputs are all given with their names for clarity, but this is not necessary
#the nearest neighbor distance (in-plane) a_nn is optional, and overrides the lat_con variable
#  meaning the value of lat_con is unused
atoms=make_graphene(alignment='AB',cell_type='rect',n_layer=3,
		    n_1=3,n_2=3,lat_con=0.0,a_nn=1.5,sep=1.0,
		    sym='C',mass=12)
ase.visualize.view(atoms)
```

This example gives the same result as the above, but specifies the relevant properties per layer using lists (instead of with a single value to be assumed for all layers). Note that while there are 3 layers, the alignment and interlayer separations lists are only `n_layer-1=3-1=2` elements long, since alignment is relative to the bottom layer and interlayer separation is relative to the layer below. However, the atomic symbol and mass lists are `n_layer` long, since it makes sense to specify these quantities for every layer. See the documentation int the section above for more information on this.

```python
import ase
from ase.visualize import view
#the comments from the above example apply here as well
atoms=make_graphene(alignment=['AB','AA'],cell_type='rect',n_layer=3,
		    n_1=3,n_2=3,lat_con=0.0,a_nn=1.5,sep=[1.0,1.0],
		    sym=['C','C','C'],mass=[12,12,12])
ase.visualize.view(atoms)
```