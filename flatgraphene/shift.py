
import copy
import numpy as np
import ase
from ase.lattice.orthorhombic import SimpleOrthorhombicFactory
from ase.visualize import view



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

def make_layer(stacking,cell_type,n_1,n_2,lat_con,z_val,sym,mass):
    """
    Creates and returns a single layer of graphene which has been shifted
    according to its stacking (AA, AB, SP) and z_value
    ---Input---
    stacking: specification of stacking for layers above first,
              ('AA','AB','SP') relative to first layer,
              numpy array of shape (n_layer-1,2), list of strings
              or single string
    cell_type: unit cell type, 'rect' or 'hex', string
    n_1: number of unit cells in x direction, integer
    n_2: number of unit cells in y direction, integer
    lat_con: in-plane lattice constant, float [Angstroms]
    z_val: z coordinates of atoms in layer, float [Angstroms]
    ---Output---
    layer: a single graphene layer shifted appropriately, ASE object
    """

    a_nn=lat_con/(2*np.sin(np.pi/3)) #compute nearest neighbor distance
    vert_shift=np.array([0.,0.,z_val])
    horz_shift=np.zeros(3)
    if (cell_type=='rect'):
        fact=GrapheneFactoryRectangular()
        if (type(stacking).__module__ == 'numpy'):
            horz_shift[0:2]=stacking #insert x, y components
        elif (stacking=='AA'):
            horz_shift=np.array([0.,0.,0.])
        elif (stacking=='AB'):
            horz_shift=np.array([a_nn,0.,0.])
        elif (stacking=='SP'):
            horz_shift=np.array([0.,lat_con/2,0.])
        latticeconstant_a_nn_1=np.array([3,2*np.sin(np.pi/3),1]) #lattice constants to scale unit cell from 1x1x1 such that afterwards the nearest neighbor distances are all 1
        scale_x_y_z=np.array([a_nn,a_nn,1])
        latticeconstant_scaled=tuple(scale_x_y_z*latticeconstant_a_nn_1) #scale up a_nn=1 unit cell by a_nn computed from lattice constant
        layer=fact(directions=[[1,0,0],[0,1,0],[0,0,1]],
                   size=(n_1,n_2,1),
                   latticeconstant=latticeconstant_scaled,
                   symbol=sym)
        layer.translate(horz_shift+vert_shift) #shift layer according to stacking and z_val
        n_atoms_layer=layer.get_masses().shape[0] #extract number of atoms in layer
        layer.set_masses(mass*np.ones(n_atoms_layer)) #set masses to mass (rather than ASE default by symbol)
    elif (cell_type=='hex'):
        print('ERROR: hexagonal unit cells not yet implemented')

    return layer


def make_graphene(stacking,cell_type,n_layer,n_1,n_2,lat_con,a_nn=None,sep=None,sym='C',mass=12.01,h_vac=None):
    """
    Generates untwisted, uncorrugated graphene and returns ASE atoms object
    with specified graphene's geometry
    ---Input---
    stacking: specification of stacking for layers above first,
              ('AA','AB','SP') relative to first layer,
              numpy array of shape (n_layer-1,2), list of strings
              or single string
              *NOTE*: single string inputs result in the input string
                      alternated with 'AA' for n_layers
    cell_type: unit cell type, 'rect' or 'hex', string
    n_layer: number of graphene layers, integer
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
    h_vac: height of the vacuum layer above and below outer layers, float [Angstroms]
    ---Output---
    atoms: graphene stack, ASE atoms object
    """

    n_layer=int(n_layer) #clean input

    #optionally use a_nn to override lat_con
    if (a_nn):
        lat_con=2*a_nn*np.sin(np.pi/3)

    #check errors in stacking (make into list if necessary)
    if (type(stacking).__module__ == 'numpy'):
        if (stacking.shape != (n_layer-1,2)):
            print('ERROR: specifying stacking as numpy array requires shape (n_layer-1,2)')
            return
        else:
            stacking = np.vstack((np.zeros(2),stacking)) #add "hidden" AA stacking for bottom layer
    elif (isinstance(stacking,list)):
        if (len(stacking) != (n_layer-1) ):
            print('ERROR: specifying stacking as list of strings requires list of length n_layer-1')
            return
        else:
            stacking=['AA']+stacking #add "hidden" 'AA' stacking for bottom layer
    elif (isinstance(stacking,str)):
          if (stacking not in ['AA','AB','SP']):
            print('ERROR: stacking string not in {\'AA\',\'AB\',\'SP\'}')
            return
          else: #make string input into n_layer length string
              stacking_string=stacking
              stacking=[0]*(n_layer) #includes specification of bottom layer
              for i_layer in range(len(stacking)):
                  if (i_layer%2 == 0):
                      stacking[i_layer]='AA'
                  else:
                      stacking[i_layer]=stacking_string
    else:
        print('ERROR: stacking input must be string, list of strings, or numpy \
array with shape (n_layers-1,2)')

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
                for i_sep in range(z_abs.shape[0]): #compute offsets relateive to bottom layer
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
    atoms=make_layer(stacking[0],cell_type,n_1,n_2,lat_con,z_abs[0],sym[0],mass[0]) #make monolayer as atoms
    #add layers on top on at a time
    for i_layer in range(1,n_layer):
        #add new atoms to object
        atoms+=make_layer(stacking[i_layer],cell_type,n_1,n_2,lat_con,z_abs[i_layer],sym[i_layer],mass[i_layer]) 
        #adjust z-height of simulation cell
        cur_cell=atoms.get_cell()
        cur_cell[2]=(z_abs[i_layer]+sep_input[i_layer])*np.eye(3)[:,2] #set z-height as vector, set buffer above to previous interlayer separation (fine in most cases)
        atoms.set_cell(cur_cell)

    #add vacuum layer of h_vac around outermost layers
    if h_vac: 
        #TURN OFF Z PERIODICITY WHEN PERIODICITY IS ADDRESSED
        cur_cell=atoms.get_cell()
        cur_cell[2]=z_abs[i_layer]*np.eye(3)[:,2] #not z-periodic, remove assumed periodic space
        atoms.set_cell(cur_cell)
        z_last=z_abs[-1] #height of last layer
        tot_vac=2*h_vac #total thickness of vacuum layer
        ase.build.add_vacuum(atoms,tot_vac) #all vacuum layer added to top
        coords=atoms.get_positions()
        n_atoms=coords.shape[0] #get number of atoms
        coords[:,2]+=h_vac*np.ones(n_atoms) #shift atoms up in vacuum so it's symmetric about center
        atoms.set_positions(coords) #rewrite "vacuum-centered coordinates"

    return atoms

        

if (__name__=="__main__"):
    """
    atoms=make_graphene(stacking='AB',cell_type='rect',n_layer=3,
		        n_1=3,n_2=3,lat_con=0.0,a_nn=1.5,sep=1.0,sym='C',mass=12)
    atoms=make_graphene(stacking=['AB','AA'],cell_type='rect',n_layer=3,
		        n_1=3,n_2=3,lat_con=0.0,a_nn=1.5,sep=[1.0,1.0],
                        sym=['C','C','C'],mass=[12,12,12])
    """

    atoms=make_graphene(stacking='SP',cell_type='rect',n_layer=2,
		        n_1=3,n_2=3,lat_con=0.0,a_nn=1.5,sep=2.0,h_vac=3.0)
    ase.visualize.view(atoms)
