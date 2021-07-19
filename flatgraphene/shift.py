
import copy
import numpy as np
import ase
from ase.lattice.orthorhombic import SimpleOrthorhombicFactory
from ase.visualize import view


"""
TO DO:
  -change cell height based on number of layers
  -don't forget periodic boundary conditions
  -make hexagonal unit cell generators
"""

class GrapheneFactoryRectangular(SimpleOrthorhombicFactory):
    """
     y                     horiz_right
     ^                     |
     |__>x               O---O
             diag_up- > /     \
                       O       O
     (0,0,z_box_len/2)-^
    """
    #construct a single rectangular unit cell
    horiz_right=np.array([1.,0.,0.])
    diag_up=np.array([np.cos(np.pi/3), np.sin(np.pi/3), 0])
    bravais_basis_real=np.array([ [0.,0.,0.],   #atomic basis for real system with nearest neighbor distance=1
                                    diag_up,
                                    diag_up+horiz_right,
                                    2*horiz_right])  #standard atomic basis at z=0
    bravais_basis_square=copy.deepcopy(bravais_basis_real)
    bravais_basis_square[:,0]*=1/3 #scale x values to fit in square with side length 1
    bravais_basis_square[:,1]*=1/(2*np.sin(np.pi/3)) #scale y values

    #set properties of class
    bravais_basis=bravais_basis_square
    element_basis=(0,0,0,0)

#class GrapheneFactoryHexagonal(SimpleHexagonalFactory)
######two atom basis

def make_layer(alignment,cell_type,n_1,n_2,lat_con,z_val,sym,mass):
    """
    Creates and returns a single layer of graphene which has been shifted
    according to its alignment (AA, AB, SP) and z_value
    ---Input---
    alignment: specification of alignment for layers above first,
               ('AA','AB','SP') relative to first layer, list of strings
               or a single string
               *NOTE*: single string inputs result in the input string
                       alternated with 'AA'
    cell_type: unit cell type, 'rect' or 'hex', string
    n_1: number of unit cells in x direction, integer
    n_2: number of unit cells in y direction, integer
    lat_con: in-plane lattice constant, float [Angstroms]
    z_val: z coordinates of atoms in layer, float [Angstroms]
    ---Output---
    atoms: a single graphene layer shifted appropriately, ASE object
    """

    a_nn=lat_con/(2*np.sin(np.pi/3)) #compute nearest neighbor distance
    vert_shift=np.array([0.,0.,z_val])
    if (cell_type=='rect'):
        fact=GrapheneFactoryRectangular()
        if (alignment=='AA'):
            horz_shift=np.array([0.,0.,0.])
        elif (alignment=='AB'):
            horz_shift=np.array([a_nn,0.,0.])
        elif (alignment=='SP'):
            horz_shift=np.array([0.,lat_con/2,0.])
        latticeconstant_a_nn_1=np.array([3,2*np.sin(np.pi/3),1]) #lattice constants to scale unit cell from 1x1x1 such that afterwards the nearest neighbor distances are all 1
        scale_x_y_z=np.array([a_nn,a_nn,1])
        latticeconstant_scaled=tuple(scale_x_y_z*latticeconstant_a_nn_1) #scale up a_nn=1 unit cell by a_nn computed from lattice constant
        atoms=fact(directions=[[1,0,0],[0,1,0],[0,0,1]],
                   size=(n_1,n_2,1),
                   latticeconstant=latticeconstant_scaled,
                   symbol=sym)
        n_atoms_layer=atoms.get_masses().shape[0] #extract number of atoms in layer
        atoms.set_masses(n_atoms_layer*np.ones(n_atoms_layer)) #set masses to mass (rather than ASE default by symbol)
        atoms.translate(horz_shift+vert_shift) #shift layer according to alignment and z_val
    elif (cell_type=='hex'):
        print('ERROR: hexagonal unit cells not yet implemented')

    return atoms


def make_graphene(alignment,cell_type,n_layer,n_1,n_2,lat_con,a_nn=None,sep=None,sym='C',mass=12.01,h_vac=6):
    """
    Generates untwisted, uncorrugated graphene and returns ASE atoms object
    with specified graphene's geometry
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
    sep: interlayer separation(s) for n_layer>1, n_layer-1 list of separations
         (relative to layer below) or float (uniform separations)
    sym: optional atomic symbol(s), list of length n_layer containing
         characters/strings or single character/string if same symbol
         for every layer
    mass: optional mass, list of length n_layer containing numeric values
          or single numerical value if every layer has the same mass
    h_vac: height of the vacuum layer, float [Angstroms]
    ---Output---
    atoms: graphene stack, ASE atoms object
    """

    n_layer=int(n_layer) #clean input

    #optionally use a_nn to override lat_con
    if (a_nn):
        lat_con=2*a_nn*np.sin(np.pi/3)

    #check errors in alignment (make into list if necessary)
    if (isinstance(alignment,list)):
        if (len(alignment) != (n_layer-1) ):
            print('ERROR: specifying alignment as list requires list of length n_layer-1')
            return
        else:
            alignment=['AA']+alignment #add "hidden" 'AA' alignment for bottom layer
    elif (isinstance(alignment,str)):
          if (alignment not in ['AA','AB','SP']):
            print('ERROR: alignment string not \'AA\',\'AB\',\'SP\'')
            return
          else: #make string input into n_layer length string
              alignment_string=alignment
              alignment=[0]*(n_layer) #includes specification of bottom layer
              for i_layer in range(len(alignment)):
                  if (i_layer%2 == 0):
                      alignment[i_layer]='AA'
                  else:
                      alignment[i_layer]=alignment_string
    else:
        print('ERROR: alignment input must be list (of strings), or string')

    #check errors in cell_type
    if (cell_type not in ['rect','hex']):
        print('ERROR: only rectangular and hexagonal unit cells supported')
        return

    #check errors in sep (turn into list if necessary)
    if (n_layer == 1):
        z_abs=np.array([0.0]) #monolayer z height
    elif (n_layer > 1):
        if (not sep):
            print('ERROR: multilayer systems require optional input sep')
            return
        elif (isinstance(sep,list)):
            if (len(sep) != (n_layer-1)):
                print('ERROR: specifying sep as list requires list length n_layer-1')
                return
            else:
                sep_input=np.array([0.0]+sep,dtype=float) #snuck in leading 0.0 for first layer
                z_abs=np.empty(sep_input.shape[0],dtype=float)
                for i_sep in range(z_abs.shape[0]):
                    #turn offsets from (relative to layer below) to (relative to
                    #  bottom layer)
                    z_abs[i_sep]=np.sum(sep_input[0:i_sep+1])
        elif (isinstance(sep,(float,int))):
            sep_input=np.array([0.0]+[sep]*(n_layer-1),dtype=float) #sneak in 0.0
            z_abs=sep*np.arange(n_layer) #[0, sep, 2*sep, ...]
    else:
        print('ERROR: n_layer must be a positive integer')

    #check errors in sym
    if (isinstance(sym,list)):
        if (len(sym) != n_layer):
            print('ERROR: specifying sym as list requires length n_layer')
            return
    elif (isinstance(sym,str)):
        sym=[sym]*n_layer #convert to list of length n_layer
    else:
        print('ERROR: optional sym inputs must be list of characters/strings or character/string')

    #check errors in mass
    if (isinstance(mass,list)):
        if (len(mass) != n_layer):
            print('ERROR: specifying mass as list requires length n_layer')
            return
    elif (isinstance(mass,(float,int))):
        mass=mass*np.ones(n_layer)
    else:
        print('ERROR: optional mass inputs must be list or numeric')

            
    #create specified geometry
    atoms=make_layer(alignment[0],cell_type,n_1,n_2,lat_con,z_abs[0],sym[0],mass[0]) #make monolayer as atoms
    #add layers on top on at a time
    for i_layer in range(1,n_layer):
        #add new atoms to object
        atoms+=make_layer(alignment[i_layer],cell_type,n_1,n_2,lat_con,z_abs[i_layer],sym[i_layer],mass[i_layer]) 
        #adjust z-height of simulation cell
        cur_cell=atoms.get_cell()
        cur_cell[2]=(z_abs[i_layer]+sep_input[i_layer])*np.eye(3)[:,2] #set z-height as vector, set buffer above to previous interlayer separation (fine in most cases)
        atoms.set_cell(cur_cell)

    #add vacuum layer
    ase.build.add_vacuum(atoms,h_vac)

    return atoms

        

if (__name__=="__main__"):
    """
    atoms=make_graphene(alignment='AB',cell_type='rect',n_layer=3,
		        n_1=3,n_2=3,lat_con=0.0,a_nn=1.5,sep=1.0,sym='C',mass=12)
    atoms=make_graphene(alignment=['AB','AA'],cell_type='rect',n_layer=3,
		        n_1=3,n_2=3,lat_con=0.0,a_nn=1.5,sep=[1.0,1.0],
                        sym=['C','C','C'],mass=[12,12,12])
    """
    atoms=make_graphene(alignment='SP',cell_type='rect',n_layer=2,
		        n_1=3,n_2=3,lat_con=0.0,a_nn=1.5,sep=2.0)
    ase.visualize.view(atoms)
