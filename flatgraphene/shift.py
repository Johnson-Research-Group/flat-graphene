
import copy
import numpy as np
import ase
from ase.lattice.orthorhombic import SimpleOrthorhombicFactory
from ase.lattice.triclinic import TriclinicFactory
from ase.lattice.hexagonal import HexagonalFactory
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


class GrapheneFactoryTriclinic(TriclinicFactory):
    """
    
    y              /--------/
    ^             /        /
    |__>x        /  o     /
                o--------/
          (0,0)-^
     
    """
    #construct basis for triclinic unit cell in terms of triclinic basis vectors
    bravais_basis_unit=np.array([ [0.,0.,0.],  #atomic basis for real system with lattice vector lengths=1
                                  (1/3)*np.array([1.,1.,0.])])

    #set properties of class
    bravais_basis=bravais_basis_unit
    element_basis=(0,0)


def make_layer(stacking,cell_type,n_1,n_2,lat_con,z_val,sym,mass):
    """
    Creates and returns a single layer of graphene which has been shifted
    according to its stacking and z_value
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
    sym: optional atomic symbol(s), list of length n_layer containing
         characters/strings or single character/string if same symbol
         for every layer
    mass: optional mass, list of length n_layer containing numeric values
          or single numerical value if every layer has the same mass
    ---Output---
    layer: a single graphene layer shifted appropriately, ASE object
    """
    vert_shift=np.array([0.,0.,z_val])
    horz_shift=np.zeros(3)
    if (cell_type=='rect'):
        fact=GrapheneFactoryRectangular()
        a_nn=lat_con/(2*np.sin(np.pi/3)) #compute nearest neighbor distance
        if (type(stacking).__module__ == 'numpy'):
            horz_shift[0:2]=stacking #insert x, y components
        elif (stacking=='AA'):
            horz_shift=np.array([0.,0.,0.])
        elif (stacking=='AB'):
            horz_shift=np.array([a_nn,0.,0.])
        elif (stacking=='SP'):
            horz_shift=np.array([0.,lat_con/2,0.])
        latticeconstant_scaled = tuple(np.array([3*a_nn,lat_con,1.])) #lengths of 3 lattice vectors
        layer=fact(directions=[[1,0,0],[0,1,0],[0,0,1]],
                   size=(n_1,n_2,1),
                   latticeconstant=latticeconstant_scaled,
                   symbol=sym)
        layer.translate(horz_shift+vert_shift) #shift layer according to stacking and z_val
        n_atoms_layer=layer.get_masses().shape[0] #extract number of atoms in layer
        layer.set_masses(mass*np.ones(n_atoms_layer)) #set masses to mass (rather than ASE default by symbol)
    elif (cell_type=='hex'):
        fact = GrapheneFactoryTriclinic()
        a_nn = (np.sqrt(2*(1+np.cos(np.pi/3)))/3)*lat_con #compute nearest neighbor distance
        if (type(stacking).__module__ == 'numpy'):
            horz_shift[0:2]=stacking #insert x, y components
        elif (stacking=='AA'):
            horz_shift=np.array([0.,0.,0.])
        elif (stacking=='AB'):
            horz_shift=np.array([0.,a_nn,0.])
        elif (stacking=='SP'):
            horz_shift=np.array([lat_con/2,0.,0.])
        triclinic_side_lengths = np.array([lat_con,lat_con,1.])
        triclinic_angles = np.array([90,90,60]) #[alpha, beta, gamma]
        triclinic_cell_info = tuple(np.hstack((triclinic_side_lengths,triclinic_angles)))
        layer=fact(size=(n_1,n_2,1),
                   latticeconstant=triclinic_cell_info,
                   symbol=sym)
        layer.translate(horz_shift+vert_shift) #shift layer according to stacking and z_val
        n_atoms_layer=layer.get_masses().shape[0] #extract number of atoms in layer
        layer.set_masses(mass*np.ones(n_atoms_layer)) #set masses to mass (rather than ASE default by symbol)
    return layer


def make_graphene(stacking,cell_type,n_1,n_2,lat_con,n_layer,sep,a_nn=None,sym='C',mass=12.01,h_vac=None):
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

    #check errors in cell_type
    if (cell_type not in ['rect','hex']):
        print('ERROR: only rectangular and hexagonal unit cells supported')
        return

    #optionally use a_nn to override lat_con
    if ((a_nn) and (cell_type == 'rect')):
        lat_con=2*a_nn*np.sin(np.pi/3)
    elif ((a_nn) and (cell_type == 'hex')):
        lat_con = (3/np.sqrt(2*(1+np.cos(np.pi/3))))*a_nn

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

    #check n_layer
    if (not (n_layer - int(n_layer) == 0.0 )):
        print('ERROR: n_layer must be a positive integer')
        return
    else:
        n_layer = int(n_layer) #clean input

    #check errors in sep (turn into list if necessary)
    if (not sep):
        print('ERROR: sep required even for monolayer (to specify z height of box)')
        return
    elif (isinstance(sep,list)):
        if (len(sep) != n_layer):
            print('ERROR: specifying sep as list requires list length n_layer')
            return
        else:
            sep_input=np.array([0.0]+sep,dtype=float) #sneak in leading 0.0 for first layer
            z_abs=np.empty(sep_input.shape[0],dtype=float)
            for i_sep in range(z_abs.shape[0]): #compute offsets relative to bottom layer
                z_abs[i_sep]=np.sum(sep_input[0:i_sep+1])
    elif (isinstance(sep,(float,int))):
        sep_input=np.array([0.0]+[sep]*(n_layer-1),dtype=float) #sneak in 0.0
        z_abs=sep*np.arange(n_layer+1) #[0, sep, 2*sep, ...]
    print(z_abs)

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

    #create specified geometry layer by layer
    #atoms=make_layer(stacking[0],cell_type,n_1,n_2,lat_con,z_abs[0],sym[0],mass[0]) #make monolayer as atoms
    #add layers on top on at a time
    for i_layer in range(0,n_layer):
        if (i_layer == 0): #create new atoms object
            atoms = make_layer(stacking[i_layer],cell_type,n_1,n_2,lat_con,z_abs[i_layer],sym[i_layer],mass[i_layer]) 
        else: #add atoms to object
            atoms += make_layer(stacking[i_layer],cell_type,n_1,n_2,lat_con,z_abs[i_layer],sym[i_layer],mass[i_layer]) 
        #adjust z-height of simulation cell
        cur_cell=atoms.get_cell()
        cur_cell[2]=z_abs[i_layer+1]*np.eye(3)[:,2] #set z-height as vector, set buffer above to previous interlayer separation (fine in most cases)
        atoms.set_cell(cur_cell)

    #add vacuum layer of h_vac around outermost layers
    if (h_vac): 
        #TURN OFF Z PERIODICITY WHEN PERIODICITY IS ADDRESSED
        cur_cell=atoms.get_cell()
        cur_cell[2]=z_abs[i_layer]*np.eye(3)[:,2] #not z-periodic, remove assumed periodic space
        cur_cell[2] += 2*h_vac*np.eye(3)[:,2] #add full thickness of vacuum on top
        atoms.set_cell(cur_cell)
        coords=atoms.get_positions()
        n_atoms=coords.shape[0] #get number of atoms
        coords[:,2]+=h_vac*np.ones(n_atoms) #shift atoms up so vacuum symmetric about center
        atoms.set_positions(coords) #rewrite "vacuum-centered" coordinates

    return atoms

        

if (__name__=="__main__"):
    #example to modify when working on module
    atoms=make_graphene(stacking='AA',cell_type='hex',n_layer=2,
		        n_1=1,n_2=1,lat_con=0.0,a_nn=1.5,sep=2.0)
    print(atoms.get_positions())
    ase.visualize.view(atoms)













