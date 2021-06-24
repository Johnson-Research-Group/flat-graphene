
import copy
import numpy as np
import ase
from ase.lattice.orthorhombic import SimpleOrthorhombicFactory
from ase.visualize import view


"""
TO DO:
  -rethink how to assemble multilayers (can't handle AA with differrent atoms
   in different layers)
    -solution: only have AAGrapheneFactoryRectangular() and then assemble more complicated Atoms objects by joining AA monolayers
    -this allows n-layer specification with optional property list for each layer (atom type, mass, etc.)
  -don't forget periodic boundary conditions
  -make hexagonal unit cell generators
"""

#how to deal with multiple layers?

class AAGrapheneFactoryRectangular(SimpleOrthorhombicFactory):
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

class AAGrapheneFactoryRectangular(SimpleOrthorhombicFactory):
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

def make_graphene(alignment,cell_type,n_layer,n_x,n_y,lat_con,a_nn=None,sep=None,disp=False):
    """
    Writes specified graphene's geometry to file.
    ---Input---
    alignment: 'AA' or 'AB', string
    cell_type: unit cell type, 'rect' or 'hex', string
    n_layer: number of graphene layers (1 or 2), integer
    n_x: number of unit cells in x direction, integer
    n_y: number of unit cells in y direction, integer
    lat_con: in-plane lattice constant, float [Angstroms]
    a_nn: optional argument to specify distance between
          nearest neighbors and override lat_con, float [Angstroms]
    sep: used to determine interlayer separation when n_layer=2
    ---Output---
    NONE: geometry is written to specified file is written
    """
    #catch input errors
    if (alignment not in ['AA','AB']):
        print('ERROR: alignment not \'AA\'or \'AB\'')
        return
    if (n_layer not in [1,2]):
        print('ERROR: number of layers not equal to 1 or 2')
        return
    if (cell_type not in ['rect','hex']):
        print('ERROR: only rectangular and hexagonal unit cells supported')
        return
    
    #optionally use a_nn to override lat_con
    if (a_nn):
        lat_con=2*a_nn*np.sin(np.pi/3)

    #alias classes
    AA=AAGrapheneFactoryRectangular()

    #create specified geometry
    if (n_layer==1):
        if (cell_type=='rect'):
            latticeconstant_a_nn_1=np.array([3,2*np.sin(np.pi/3),1]) #lattice constants to scale unit cell from 1x1x1 such that afterwards the nearest neighbor distances are all 1
            scale_x_y_z=np.array([lat_con/(2*np.sin(np.pi/3)),lat_con/(2*np.sin(np.pi/3)),1])
            latticeconstant_scaled=tuple(scale_x_y_z*latticeconstant_a_nn_1) #scale up a_nn=1 unit cell by a_nn computed from lattice constant
            atoms=AA(directions=[[1,0,0],[0,1,0],[0,0,1]],
                    size=(n_x,n_y,1),
                    symbol='C',
                    latticeconstant=latticeconstant_scaled)
        elif (cell_type=='hex'):
                pass
    elif (n_layer==2):
        if (not sep):
            print('ERROR: bilayer requires defined interlayer separation (sep)')
            return
        if (cell_type=='rect'):
            if (alignment == 'AA'):

                latticeconstant_a_nn_1=np.array([3,2*np.sin(np.pi/3),1]) #lattice constants to scale unit cell from 1x1x1 such that afterwards the nearest neighbor distances are all 1
                scale_x_y_z=np.array([lat_con/(2*np.sin(np.pi/3)),lat_con/(2*np.sin(np.pi/3)),sep])
                latticeconstant_scaled=tuple(scale_x_y_z*latticeconstant_a_nn_1) #scale up a_nn=1 unit cell by a_nn computed from lattice constant
                atoms=AA(directions=[[1,0,0],[0,1,0],[0,0,1]],
                        size=(n_x,n_y,2),
                        symbol='C',
                        latticeconstant=latticeconstant_scaled)
            elif (alignment == 'AB'):
                pass
        elif (cell_type=='hex'):
            if (alignment == 'AA'):
                pass
            elif (alignment == 'AB'):
                pass
            pass

    return atoms

        

if (__name__=="__main__"):
    #def make_graphene(alignment,cell_type,n_layer,n_x,n_y,lat_con,a_nn=None,sep=None,disp=False):
    atoms=make_graphene(alignment='AA',cell_type='rect',n_layer=2,
                        n_x=10,n_y=10,lat_con=0.0,a_nn=1.0,sep=1.0,disp=True)
    ase.visualize.view(atoms)
