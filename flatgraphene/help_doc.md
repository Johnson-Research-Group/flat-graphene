
View with markdown syntax highlighting at https://github.com/Johnson-Research-Group/flat-graphene
at location `flatgraphene/help_doc.md`.

## Documentation
Non-twisted (shfited) graphene may be created using the `shift.make_graphene` function, which returns an ASE atoms object.
Parameters include cell type (rectangular or hexagonal), alignment (A,B,C,SP), number of unit cells (in-plane), interlayer spacing, and lattice constant.

```python
#in flatgraphene/shift.py
def make_graphene(stacking,cell_type,n_1,n_2,lat_con,n_layer,sep,a_nn=None,sym='C',
                  mass=12.01,mol_id=None,h_vac=None):
    """
    Generates untwisted, uncorrugated graphene and returns ASE atoms object
    with specified graphene's geometry
    ---Input---
    stacking: specification of stacking for layers above first,
              ('A','B','C','SP') relative to origin of cell, 
              single string, or list of strings, or numpy array of shape
	      (n_layer,2)
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
    mol_id: optional molecular IDs for each layer, integer or list of integers
            of length n_layer
    h_vac: height of the vacuum layer above and below outer layers, float [Angstroms]
    ---Output---
    atoms: ASE atoms object
    """
```

Twisted multi-layer graphene may be similarly created via `twist.make_graphene`.
Parameters are similar to the shifted case, but note that twist angle **is not** an input parameters, rather the `p,q` specification from "Electronic structure of turbostratic graphene" by Shallcross is preferred.
Use the function `twist.find_p_q(theta)` (documentation below) to find the `p,q` corresponding the desired twist angle, then use computed `p,q` as the input to the graphene generator.

```python
#in flatgraphene/twist.py
def make_graphene(cell_type,p,q,lat_con,n_layer,sep,a_nn=None,sym='C',
                  mass=12.01,mol_id=None,h_vac=None):
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
    mol_id: optional molecular IDs for each layer, integer or list of integers
            of length n_layer
    h_vac: height of the vacuum layer above and below outer layers, float [Angstroms]
    ---Output---
    atoms: graphene stack, ASE atoms object
    """
```

Compute the integers `p,q` corresponding to input twist angle using `twist.find_p_q`.
If no match to input angle is found with default parameters, consider increasing `q_max` or decreasing `a_tol`.


```python
def find_p_q(theta_deg,q_max=100,a_tol=1e-2):
    """
    Computes the p_q that generate a twist of theta degrees
    ---Inputs---
    theta_deg : {float}
        desired twist angle, [angular degrees]
    q_max : {integer}
        q >= p > 0, so q_max controls how many pairs are checked
    a_tol : {float}
        acceptable absolute difference between computed angle and desired angle
    ---Outputs---
    p : {integer}
        p value from "Electronic structure of turbostratic graphene" by Shallcross et al
    q : {integer}
        q value from "Electronic structure of turbostratic graphene" by Shallcross et al
    """
```
