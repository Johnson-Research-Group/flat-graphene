
import copy
import numpy as np
import ase
from flatgraphene import shift
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


def find_p_q(theta_deg,q_max=100,a_tol=1e-2):
    """
    Computes the p_q that generate a twist of theta radians
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
    for p in range(1,q_max+1):
        for q in range(p,q_max+1):
            theta_comp_rad = np.arccos((3*np.power(q,2) - np.power(p,2))/(3*np.power(q,2) + np.power(p,2)))
            theta_comp_deg = theta_comp_rad*(180/np.pi)
            if (np.isclose(theta_comp_deg,theta_deg,atol=a_tol)):
                return p, q, theta_comp_deg
    print('ERROR: was not able to find (p,q) that give desired twist angle, consider increasing optional argument q_max')
    return

    
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
        twist angle, radians
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

    #DETERMINE SIZE OF OVERSIZED UNTWISTED LAYER
    unit_layer = shift.make_layer('A','hex',1,1,lat_con,1.0,'C',12.01) #layer with first lattice vector in x direction
    unit_cell = unit_layer.get_cell() #extract box/cell vectors

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

    #GENERATE OVERSIZED NONTWISTED LAYER (with first box vector along [1,0,0])
    nontwisted_layer = shift.make_layer('A','hex',n_x,n_y,lat_con,0.0,'C',12.01) #layer with first lattice vector in x direction; symbol, mass, etc. will be overwritten later and are irrelevant
    theta_nontwisted_less = np.arccos(np.dot(n_loc,np.eye(2)[0])/np.linalg.norm(n_loc)) #angle between unrotated layer with lattice vector along [1,0,0] and lesser rotated layer with first lattice vector along n_loc

    #CONTRUCT LESS TWISTED LAYER
    less_twisted_layer = copy.deepcopy(nontwisted_layer)
    #shift atoms so that cell vectors do not start at left edge of oversized layer
    unshifted_positions = less_twisted_layer.get_positions()
    num_atoms_unshifted = unshifted_positions.shape[0]
    shift_vec = np.zeros(3)
    x_shift_multiplier = np.ceil(((rot_mat_60@n_loc)[1]/np.tan(np.pi/3)-(rot_mat_60@n_loc)[0])/lat_vec_mat_nontwisted[0,0])
    x_shift_multiplier = int(x_shift_multiplier)
    shift_vec[:2] = x_shift_multiplier*lat_vec_mat_nontwisted[:,0]
    shift_array = np.tile(shift_vec,(num_atoms_unshifted,1))
    shifted_positions = unshifted_positions - shift_array
    less_twisted_layer.set_positions(shifted_positions)
    #set proper lattice vectors
    less_twisted_cell = np.eye(3)
    less_twisted_cell[:2,:2] = np.array([n_loc,rot_mat_60@n_loc])
    less_twisted_layer.set_cell(less_twisted_cell)
    remove_atoms_outside_cell(less_twisted_layer) #trim extra atoms outside cell

    #CONSTRUCT MORE TWISTED LAYER
    more_twisted_layer = copy.deepcopy(nontwisted_layer)
    #shift atoms so that cell vectors do not start at left edge
    unshifted_positions = more_twisted_layer.get_positions()
    num_atoms_unshifted = unshifted_positions.shape[0]
    x_shift_multiplier = np.ceil((rot_mat_60@m_loc)[1]/(np.tan(np.pi/3)*lat_vec_mat_nontwisted[0,0]))
    x_shift_multiplier = np.ceil(((rot_mat_60@m_loc)[1]/np.tan(np.pi/3)-(rot_mat_60@m_loc)[0])/lat_vec_mat_nontwisted[0,0])
    x_shift_multiplier = int(x_shift_multiplier)
    shift_vec = np.zeros(3)
    shift_vec[:2] = x_shift_multiplier*lat_vec_mat_nontwisted[:,0]
    shift_array = np.tile(shift_vec,(num_atoms_unshifted,1))
    shifted_positions = unshifted_positions - shift_array
    more_twisted_layer.set_positions(shifted_positions)
    #set proper lattice vectors
    m_loc_rot_60 = rot_mat_z(np.pi/3)[:2,:2]@m_loc #second box vector of more twisted cell
    more_twisted_cell = np.eye(3)
    more_twisted_cell[:2,:2] = np.array([m_loc,m_loc_rot_60])
    more_twisted_layer.set_cell(more_twisted_cell)
    remove_atoms_outside_cell(more_twisted_layer) #trim extra atoms outside cell

    n_less_twisted = (less_twisted_layer.get_positions()).shape[0]
    n_more_twisted = (more_twisted_layer.get_positions()).shape[0]
    if (n_less_twisted != n_more_twisted):
        print('ERROR: number of atoms in two different layers is not the same, please open an issue on GitHub')
        return

    #ROTATE BOTH LAYERS TO ALIGN SUCH THAT FIRST LATTICE VECTOR ALONG [1,0,0]
    rotate_to_standard(less_twisted_layer)
    rotate_to_standard(more_twisted_layer)

    return less_twisted_layer, more_twisted_layer, theta


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
    #purify inputs
    p = int(p)
    q = int(q)

    #optionally use a_nn to override lat_con
    if ((a_nn) and (cell_type == 'rect')):
        lat_con=2*a_nn*np.sin(np.pi/3)
    elif ((a_nn) and (cell_type == 'hex')):
        lat_con = (3/np.sqrt(2*(1+np.cos(np.pi/3))))*a_nn

    #check n_layer
    if (not (n_layer - int(n_layer) == 0.0 )):
        print('ERROR: n_layer must be a positive integer')
        return
    else:
        n_layer = int(n_layer) #clean input

    #check errors in sep (turn into list if necessary)
    if (not sep):
        print('ERROR: parameter sep required (even for monolayer, to specify z height of box)')
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

    #make two layers
    less_twisted_layer, more_twisted_layer, theta = make_unique_twisted_layers(p,q,lat_con)
    n_atoms_layer = (less_twisted_layer.get_positions()).shape[0] #number of atoms in single layer

    #create specified geometry layer by layer (assumes stacking less-more-less-...)
    #add layers on top one at a time
    for i_layer in range(0,n_layer):
        #make copy of correct layer  based on layer number
        if ((i_layer+1)%2 == 0): #odd layer number
            cur_layer = copy.deepcopy(less_twisted_layer)
        else: #even layer number (1...n_layer)
            cur_layer = copy.deepcopy(more_twisted_layer)

        #overwrite dummy properties of layer set in make_unique_twisted
        cur_positions = cur_layer.get_positions()
        cur_positions[:,2] = z_abs[i_layer]*np.ones(n_atoms_layer) #set proper z value for layer
        cur_layer.set_positions(cur_positions)
        
        cur_layer.symbols[:] = sym[i_layer] #set symbol
        cur_layer.set_masses(mass[i_layer]*np.ones(n_atoms_layer)) #set masses (to overwrite default symbol mass)

        #create stack or add layer to stack
        if (i_layer == 0): #create new atoms object
            atoms = copy.deepcopy(cur_layer)
        else:
            atoms += cur_layer

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
    #example(s) to modify when working on module
    #p_found, q_found, theta_comp = find_p_q(21.79)
    p_found, q_found, theta_comp = find_p_q(9.43)
    print('generating system with a {:.2f} degree twist'.format(theta_comp))
    test_sheet = make_graphene(cell_type='hex',n_layer=2,
                               p=p_found,q=q_found,lat_con=0.0,a_nn=1.5,
                               sep=3.4,sym=['C','Ne'],mass=[12.01,12.02],h_vac=3)
    ase.visualize.view(test_sheet)
                        
