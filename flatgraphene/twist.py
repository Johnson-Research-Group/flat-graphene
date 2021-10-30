
import copy
import numpy as np
import ase
import shift
from ase.lattice.orthorhombic import SimpleOrthorhombicFactory
from ase.lattice.triclinic import TriclinicFactory
from ase.lattice.hexagonal import HexagonalFactory
from ase.visualize import view

def rot_mat_z(theta):
    """
    Generates a 3D rotation matrix which results in rotation of
    theta about the z axis
    """
    rot_mat = np.array([[np.cos(theta), -np.sin(theta), 0],
                        [np.sin(theta), np.cos(theta), 0],
                        [0, 0, 1]])
    return rot_mat

def cut_out_of_bounds(atoms_object):
    temp_atoms_object = copy.deepcopy(atoms_object)
    temp_atoms_object.wrap()
    original_positions = atoms_object.get_positions()
    wrapped_positions = temp_atoms_object.get_positions()
    in_bounds_positions = [original_positions[i] for i in range(original_positions.shape[0])
                           if (original_positions[i] == wrapped_positions[i])] #only unwrapped atoms were in bounds
    atoms_object.set_positions(in_bounds_positions)
    return atoms_object


def make_graphene(cell_type,p,q,lat_con,n_layer,sep,a_nn=None,sym='C',mass=12.01,h_vac=None):
    """
    Generates twisted, uncorrugated graphene and returns ASE atoms object
    with specified graphene's geometry
    ---Input---
    cell_type: unit cell type, 'rect' or 'hex', string
    q : XXXX, integer diophantine sort of thing from Shallcross
    q : XXXX, integer ''
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
    #purify inputs
    p = int(p)
    q = int(q)

    #optionally use a_nn to override lat_con
    if ((a_nn) and (cell_type == 'rect')):
        lat_con=2*a_nn*np.sin(np.pi/3)
    elif ((a_nn) and (cell_type == 'hex')):
        lat_con = (3/np.sqrt(2*(1+np.cos(np.pi/3))))*a_nn

    #compute integer vectors m, n
    n_full = np.array([p+3*q,-2*p]) + np.array([2*p,-p+3*q])
    m_full = np.array([-p+3*q,2*p]) + np.array([-2*p,p+3*q])
    n = n_full/np.gcd.reduce(n_full) #reduce to smallest integer vector of same direction
    m = m_full/np.gcd.reduce(m_full)
    n = n.astype(int) #convert back into integer vectors
    m = m.astype(int)

    #compute twist angle
    theta = np.arccos((3*np.power(q,2) - np.power(p,2))/(3*np.power(q,2) + np.power(p,2))) #defined by p, q

    #compute commensuration vectors
    delta = int(3/np.gcd(p,3))
    gamma = np.gcd(int(3*q+p),int(3*q-p))
    if (delta == 1):
        t_1 = (1/gamma)*np.array([[p + 3*q],
                                [-2*p]])
        t_2 = (1/gamma)*np.array([[2*p],
                                [-p + 3*q]])
    elif (delta == 3):
        t_1 = (1/gamma)*np.array([[-p - q],
                                [2*q]])
        t_2 = (1/gamma)*np.array([[2*q],
                                [-p + q]])

    #---just for testing befor for loop over layers goes in---#
    z_temp = 0.0
    sym_temp = 'C'
    mass_temp = 12.01
    #---just for testing befor for loop over layers goes in---#

    """
    NOTES ON CURRENT APPROACH:
    - I can't just generate a raw cell then rotate, wrap, and rotate, wrap because the lattice vectors of the untwisted system and the less twisted system are not the same length
    - at least I can't without first generating an initially larger system, then cutting atoms out
    - generating initially large system seems okay for now, but need to make sure that end of first lattice vector of less twisted system aligns perfectly with an atoms as in Shallcross paper
    - you'll want to use n vector to determine dyanmically how much overage to create in the non-twisted cell
    ALTERNATIVELY:
    - could I find a way to systematically/analytically generate the first, less twisted cell, then proceed to use the rotate and wrap strategy?
    """

    #generate sheet will first box vector along x direction
    #def make_layer(stacking,cell_type,n_1,n_2,lat_con,z_val,sym,mass):
    magic_size = 5 #magic variable, get rid of by setting dynamically large cell via check of n vector
    nontwisted_sheet = shift.make_layer('A','hex',magic_size,magic_size,lat_con,z_temp,sym_temp,mass_temp) #sheet with first lattice vector in x direction
    nontwisted_cell = nontwisted_sheet.get_cell() #extract box/cell vectors

    #get angle between sheet with first lattice vector along x, and sheet with first lattice vector
    lat_vec_mat_nontwisted = nontwisted_cell[:2,:2].T/magic_size #matrix of 2D lattice vectors from box vectors, one per column
    n_loc = lat_vec_mat_nontwisted@n #convert from number of lattice hops (vector n) to actual spatial position

    theta_nontwisted_less = np.arccos(np.dot(n_loc,np.eye(2)[0])/np.linalg.norm(n_loc)) #angle between unrotated sheet and lesser rotated sheet

    #contruct less twisted sheet
    max_n = np.max(n) #ensures cell is always "square" (equal repititions along both lattice vectors)
    less_twisted_sheet = copy.deepcopy(nontwisted_sheet)
    less_twisted_cell = np.eye(3)
    less_twisted_cell[:2,:2] = lat_vec_mat_nontwisted@(max_n*np.eye(2))
    #less_twisted_cell = (rot_mat_z(theta_nontwisted_less)@nontwisted_cell.T).T #rotate nontwisted cell vectors
    less_twisted_sheet.set_cell(less_twisted_cell)
    #less_twisted_sheet.wrap() #wrap atoms back across new box vectors
                              

    #normal_sheet = generate one normal layer of graphene as atoms object, n_max=m_max=max(n vector)
    #twisted_1 = copy(normal_sheet)
    #rotate twisted_1's lattice vectors by angle between n vector and [1,0]
    #wrap atoms of twisted_1 since lattice vectors have changed
    #twisted_2 = copy(twisted_1)
    #rotate twisted_2's lattice vectors by theta
    #wrap atoms of twisted_2 since lattice vectors have changed


    print('theta (degrees):',theta*(180/np.pi))
    print('n:',n)
    print('m:',m)
    print('theta_normal_1 (degrees):',theta_nontwisted_less*(180/np.pi))
    ase.visualize.view(nontwisted_sheet)
    ase.visualize.view(less_twisted_sheet)

if (__name__=="__main__"):
    #example to modify when working on module
    make_graphene(cell_type='hex',n_layer=1,
		  p=1,q=3,lat_con=0.0,a_nn=1.5,
                  sep=1)
