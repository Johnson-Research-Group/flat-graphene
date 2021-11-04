
Python module to generate flat (uncorrugated) graphene configurations with one or more layers, as well as twisted bilayer graphene.

## Documentation

Non-twisted (shfited) graphene may be created using the `shift.make_graphene` function, which returns an ASE atoms object. Parameters include cell type (rectangular or hexagonal), alignment (A,B,SP), number of unit cells (in-plane), interlayer spacing, and lattice constant.

```python
#in flatgraphene/shift.py
def make_graphene(stacking,cell_type,n_1,n_2,lat_con,n_layer,sep,
                  a_nn=None,sym='C',mass=12.01,h_vac=None):
    """
    Generates untwisted, uncorrugated graphene and returns ASE atoms object
    with specified graphene's geometry
    ---Input---
    stacking: specification of stacking for layers above first,
              ('A','B','SP') relative to first layer,
              single string, or list of strings, or numpy array of shape
			  (n_layer-1,2)
              *NOTE*: single string inputs result in the input string
                      alternated with 'A' for n_layers
    cell_type: unit cell type, 'rect' or 'hex', string
    n_1: number of unit cells in x direction, integer
    n_2: number of unit cells in y direction, integer
    lat_con: in-plane lattice constant, float [Angstroms]
    n_layer: number of graphene layers, integer
    sep: interlayer separation(s) for n_layer>1, n_layer list of separations
         (relative to layer below) or float (uniform separations)
         last element specifies distance between top layer and top surface of box
    a_nn: optional argument to specify distance between
          nearest neighbors and override lat_con, float [Angstroms]
    sym: optional atomic symbol(s), list of length n_layer containing
         characters/strings or single character/string if same symbol
         for every layer
    mass: optional mass, list of length n_layer containing numeric values
          or single numerical value if every layer has the same mass
    h_vac: height of the vacuum layer above and below outer layers, float [Angstroms]
    ---Output---
    atoms: graphene stack, ASE atoms object
    """
```

Twisted multi-layer graphene may be similarly created via `twist.make_graphene`. Parameters are similar to the shifted case, but note that twist angle **is not** an input parameters, rather the `p,q` specification from "Electronic structure of turbostratic graphene" by Shallcross is preferred. One can use the function `twist.find_p_q(theta)` to find the `p,q` corresponding the desired twist angle, then use these as the input to the graphene generator.

```python
#in flatgraphene/twist.py
def make_graphene(cell_type,p,q,lat_con,n_layer,sep,a_nn=None,sym='C',
                  mass=12.01,h_vac=None):
    """
    Generates twisted, uncorrugated graphene and returns ASE atoms object
    with specified graphene's geometry
    NOTE: This function does not allow the input of twist angle, rather
          the user should determine the proper (p,q) via the provided function
          find_q_p(theta), then use the computed (p,q) here.
    ---Input---
    cell_type: unit cell type, 'rect' or 'hex', string
    p : p value from "Electronic structure of turbostratic graphene" by Shallcross et al, integer 
    q : q value from "Electronic structure of turbostratic graphene" by Shallcross et al, integer 
    lat_con: in-plane lattice constant, float [Angstroms]
    n_layer: number of graphene layers, integer
    sep: interlayer separation(s) for n_layer>1, n_layer list of separations
         (relative to layer below) or float (uniform separations)
         last element specifies distance between top layer and top surface of box
    a_nn: optional argument to specify distance between
          nearest neighbors and override lat_con, float [Angstroms]
    sym: optional atomic symbol(s), list of length n_layer containing
         characters/strings or single character/string if same symbol
         for every layer
    mass: optional mass, list of length n_layer containing numeric values
          or single numerical value if every layer has the same mass
    h_vac: height of the vacuum layer above and below outer layers, float [Angstroms]
    ---Output---
    atoms: graphene stack, ASE atoms object
    """
```

## Examples

### Shift
This example creates an ABA trilayer graphene system using. Any stacking designation (A,B,SP) is defined relative to the bottom layer in the stack. The bottom layer is assumed A, and when one string input is used for stacking, that input is alternated with A stacked layers.
```python
import ase
from ase.visualize import view
import flatgraphene as fg
#note the inputs are all given with variable name for clarity,
#  but this is not necessary
#the nearest neighbor distance (in-plane) a_nn is optional, and
#  overrides the lat_con variable  meaning the value of lat_con is unused
atoms=fg.shift.make_graphene(stacking='B',cell_type='rect',n_1=3,n_2=3,
                             lat_con=0.0,n_layer=3,sep=1.0,a_nn=1.5,
                             sym='C',mass=12,h_vac=3.0)
ase.visualize.view(atoms)
```

This example gives the same result as the above, but specifies the relevant properties per layer using lists (instead of with a single value to be assumed for all layers).
Note that while there are 3 layers, the alignment and lists are only `n_layer-1=3-1=2` elements long, since alignment is relative to the bottom layer. Interlayer separation is relative to the layer below, and the last element is used to define spacing between top layer and top of supercell box.
Atomic symbol and mass lists are `n_layer` long, since it makes sense to specify these quantities for every layer.
See the documentation in the section above for more information on this.

```python
import ase
from ase.visualize import view
import flatgraphene as fg
#the comments from the above example apply here as well
atoms=fg.shift.make_graphene(stacking=['AB','AA'],cell_type='rect',n_1=3,n_2=3,
                             lat_con=0.0,n_layer=3,sep=[1.0,1.0,1.0],a_nn=1.5,
                             sym=['C','C','C'],mass=[12,12,12],h_vac=3.0)
ase.visualize.view(atoms)
```

### Twist
This example creates a 21.79 degree twisted system by first computing the proper `p,q`, then using these as inputs to `make_graphene`. All of the properties from the shifted case which can be set here also allow the same variety of input formats (scalar, list, etc.).
```python
import ase
from ase.visualize import view
import flatgraphene as fg
p_found, q_found, theta_comp = fg.twist.find_p_q(21.79)
atoms=fg.twist.make_graphene(cell_type='hex',n_layer=2,
                             p=p_found,q=q_found,lat_con=0.0,a_nn=1.5,
                             sep=3.35,h_vac=3)
ase.visualize.view(atoms)
```
