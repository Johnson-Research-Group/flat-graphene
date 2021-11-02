
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


def remove_atoms_outside_cell(atoms_object):
    """
    Deletes all atoms outside of the box/cell of atoms_object
    ---Inputs---
    atoms_object : {ASE atoms object}
    ---Outputs---
    NONE, atoms_object is edited directly
    """
    original_positions = atoms_object.get_positions()
    temp_atoms_object = copy.deepcopy(atoms_object)
    temp_atoms_object.wrap()
    wrapped_positions = temp_atoms_object.get_positions()
    out_of_bounds_indices = [atoms_object[i].index for i in range(original_positions.shape[0])
                           if (not np.allclose(original_positions[i],wrapped_positions[i]))]
    del atoms_object[out_of_bounds_indices] #remove out of bounds atoms by known index


def rotate_to_standard(atoms_object):
    """
    Rotates a cell/box such that it's first lattice vector is along [1,0,0]
    ---Inputs---
    atoms_object : {ASE atoms object}
    ---Outputs---
    NONE, atoms_object is edited directly
    """
    #cell vectors
    unrotated_cell_vecs = atoms_object.get_cell() #cell vectors are rows
    first_lat_vec = unrotated_cell_vecs[0]
    if (np.isclose(first_lat_vec[2],0)):
        theta_back = np.arccos(np.dot(first_lat_vec,np.eye(3)[0])/np.linalg.norm(first_lat_vec))
        print('theta_back:',theta_back)
    else:
        print('ERROR: expected unrotated cell vector with only nonzero x, y components')
    rot_mat_back = rot_mat_z(theta_back).T #matrix which rotates vectors CLOCKWISE by theta_back
    rotated_cell_vecs = (rot_mat_back@unrotated_cell_vecs.T).T
    atoms_object.set_cell(rotated_cell_vecs)

    #positions
    unrotated_positions = atoms_object.get_positions()
    rotated_positions = (rot_mat_back@unrotated_positions.T).T
    atoms_object.set_positions(rotated_positions)
    

def make_unique_twisted_layers(p,q,lat_con):
    """
    Generates and returns the two unique layers in twisted bilayer graphene
    ---Inputs---
    p : {integer}
        p value from "Electronic structure of turbostratic graphene" by Shallcross et al
    q : {integer}
        q value from "Electronic structure of turbostratic graphene" by Shallcross et al
    lat_con : {float}
        lattice constant of a single layer
    ---Outputs--
    less_twisted_layer : {ASE atoms object}
        atoms object corresponding to layer which is less twisted compared to untwisted cell
    more_twisted_layer : {ASE atoms object}
        atoms object corresponding to layer which is more twisted compared to untwisted cell
    theta : {float}
        twist angle
    """
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

    #DETERMINE SIZE OF OVERSIZED UNTWISTED SHEET
    unit_sheet = shift.make_layer('A','hex',1,1,lat_con,1.0,'C',12.01) #sheet with first lattice vector in x direction
    unit_cell = unit_sheet.get_cell() #extract box/cell vectors

    lat_vec_mat_nontwisted = unit_cell[:2,:2].T #matrix of 2D lattice vectors from box vectors, one per column
    n_loc = lat_vec_mat_nontwisted@n #convert from number of lattice hops (vector n) to actual spatial position, also first box vector of less twisted cell
    m_loc = lat_vec_mat_nontwisted@m #convert from number of lattice hops (vector m) to actual spatial position, also first box vector of more twisted cell
    n_loc_10 = np.linalg.norm(n_loc)*np.eye(2)[0] #vector of length n_loc in direction [1,0]
    rot_mat_60 = rot_mat_z(np.pi/3)[:2,:2] #matrix to generate 60 degree rotation about z-axis
    x_range = np.ceil((n_loc_10 + rot_mat_60@n_loc_10)[0]) #size of necessary oversized system in x direction (a bit overkill, but it works), [A]
    y_range = (m_loc + rot_mat_60@m_loc)[1] #size of necessary oversized system in y direction, [A]
    n_x = np.ceil(x_range/lat_vec_mat_nontwisted[0,0]) + 1
    n_y = np.ceil(y_range/lat_vec_mat_nontwisted[1,1]) + 1 
    n_x = int(n_x) #should be able to do int cast with ceil, but oh well
    n_y = int(n_y)

    #GENERATE OVERSIZED NONTWISTED SHEET (with first box vector along [1,0,0])
    nontwisted_sheet = shift.make_layer('A','hex',n_x,n_y,lat_con,0.0,'C',12.01) #sheet with first lattice vector in x direction; symbol, mass, etc. will be overwritten later and are irrelevant
    theta_nontwisted_less = np.arccos(np.dot(n_loc,np.eye(2)[0])/np.linalg.norm(n_loc)) #angle between unrotated sheet with lattice vector along [1,0,0] and lesser rotated sheet with first lattice vector along n_loc

    #CONTRUCT LESS TWISTED SHEET
    less_twisted_sheet = copy.deepcopy(nontwisted_sheet)
    #shift atoms so that cell vectors do not start at left edge of oversized sheet
    unshifted_positions = less_twisted_sheet.get_positions()
    num_atoms_unshifted = unshifted_positions.shape[0]
    shift_vec = np.zeros(3)
    x_shift_multiplier = np.ceil(((rot_mat_60@n_loc)[1]/np.tan(np.pi/3)-(rot_mat_60@n_loc)[0])/lat_vec_mat_nontwisted[0,0])
    x_shift_multiplier = int(x_shift_multiplier)
    shift_vec[:2] = x_shift_multiplier*lat_vec_mat_nontwisted[:,0]
    shift_array = np.tile(shift_vec,(num_atoms_unshifted,1))
    shifted_positions = unshifted_positions - shift_array
    less_twisted_sheet.set_positions(shifted_positions)
    #set proper lattice vectors
    less_twisted_cell = np.eye(3)
    less_twisted_cell[:2,:2] = np.array([n_loc,rot_mat_60@n_loc])
    less_twisted_sheet.set_cell(less_twisted_cell)
    remove_atoms_outside_cell(less_twisted_sheet) #trim extra atoms outside cell

    #CONSTRUCT MORE TWISTED SHEET
    more_twisted_sheet = copy.deepcopy(nontwisted_sheet)
    #shift atoms so that cell vectors do not start at left edge
    unshifted_positions = more_twisted_sheet.get_positions()
    num_atoms_unshifted = unshifted_positions.shape[0]
    x_shift_multiplier = np.ceil((rot_mat_60@m_loc)[1]/(np.tan(np.pi/3)*lat_vec_mat_nontwisted[0,0]))
    x_shift_multiplier = np.ceil(((rot_mat_60@m_loc)[1]/np.tan(np.pi/3)-(rot_mat_60@m_loc)[0])/lat_vec_mat_nontwisted[0,0])
    x_shift_multiplier = int(x_shift_multiplier)
    shift_vec = np.zeros(3)
    shift_vec[:2] = x_shift_multiplier*lat_vec_mat_nontwisted[:,0]
    shift_array = np.tile(shift_vec,(num_atoms_unshifted,1))
    shifted_positions = unshifted_positions - shift_array
    more_twisted_sheet.set_positions(shifted_positions)
    #set proper lattice vectors
    m_loc_rot_60 = rot_mat_z(np.pi/3)[:2,:2]@m_loc #second box vector of more twisted cell
    more_twisted_cell = np.eye(3)
    more_twisted_cell[:2,:2] = np.array([m_loc,m_loc_rot_60])
    more_twisted_sheet.set_cell(more_twisted_cell)
    remove_atoms_outside_cell(more_twisted_sheet) #trim extra atoms outside cell

    n_less_twisted = (less_twisted_sheet.get_positions()).shape[0]
    n_more_twisted = (more_twisted_sheet.get_positions()).shape[0]
    if (n_less_twisted != n_more_twisted):
        print('ERROR: number of atoms in two different layers is not the same, please open an issue on GitHub')
        return

    #ROTATE BOTH SHEETS TO ALIGN SUCH THAT FIRST LATTICE VECTOR ALONG [1,0,0]
    rotate_to_standard(less_twisted_sheet)
    rotate_to_standard(more_twisted_sheet)

    return less_twisted_sheet, more_twisted_sheet, theta


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
    #make two layers
    #assemble those layers with proper distances, etc.
    #purify inputs
    p = int(p)
    q = int(q)

    #optionally use a_nn to override lat_con
    if ((a_nn) and (cell_type == 'rect')):
        lat_con=2*a_nn*np.sin(np.pi/3)
    elif ((a_nn) and (cell_type == 'hex')):
        lat_con = (3/np.sqrt(2*(1+np.cos(np.pi/3))))*a_nn

    less_twisted_sheet, more_twisted_sheet, theta = make_unique_twisted_layers(p,q,lat_con)

    print('theta (degrees):',theta*(180/np.pi))



if (__name__=="__main__"):
    #example(s) to modify when working on module
    """
    #27.78 degrees
    make_graphene(cell_type='hex',n_layer=1,
		  p=1,q=3,lat_con=0.0,a_nn=1.5,
                  sep=1)
    """
    #??? degrees
    make_graphene(cell_type='hex',n_layer=1,
		  p=1,q=53,lat_con=0.0,a_nn=1.5,
                  sep=1)
