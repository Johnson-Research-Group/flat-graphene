
import copy
import numpy as np
#import matplotlib.pyplot as plt

"""
TO DO:
  -create "file writers" for other formats (POSCAR, etc.; only LAMMPS format currently)
"""

def make_xyz_honeycomb_rectangular(n_cell_x,n_cell_y,a_nn_input,z_box_len,lat_con=None):
    """
    Generates a honecomb lattice and writes to file in specified format.
    ---Inputs---
    n_cell_x: number of unit cells in x direction, integer
    n_cell_y: number of unit cells in y direction, integer
    a_nn_input: distance between nearest neighbors, float [Angstroms]
    z_box_len: magnitude of unit cell in z direction, z_box_len/2 gives the
               z-coordinate of the lattice, float [Angstroms]
    lat_con: lattice constant of honeycomb (overrides a_nn_input,
                      if present), float
    ---Outputs---
    super_cell_xyz: coordinates of all atoms in supercell, numpy array, shape(supercell=(n_atoms, 3)
    box_vector: 3x3 array with each row as a box vector (x y z), numpy array
    """

    rad_per_deg=np.pi/180
    if (lat_con): #use lattice_constant to compute nearest neighbor distance if present
        a_nn=lat_con/(2*np.sin(60*rad_per_deg))
    else:
        a_nn=a_nn_input #nearest neighbor distance

    """
     y                     horiz_right
     ^                     |
     |__>x               O---O
             diag_up- > /     \
                       O       O
     (0,0,z_box_len/2)-^
    """
    #construct a single rectangular unit cell
    z_height=np.array([0.,0.,z_box_len/2])
    horiz_right=a_nn*np.array([1.,0.,0.])
    diag_up=a_nn*np.array([np.cos(60*rad_per_deg),np.sin(60*rad_per_deg),0.])
    atomic_basis_xyz=np.array([ [0.,0.,0.],   #atomic basis vectors
                                diag_up,
                                diag_up+horiz_right,
                                2*horiz_right])  #standard atomic basis at z=0
    atomic_basis_xyz[:,2]+=(z_box_len/2)*np.ones(atomic_basis_xyz.shape[0]) #shift standard atomic basis to proper z height

    #lattice vectors
    a1=3*horiz_right #horizontal lattice vector
    a2=2*diag_up[1]*np.array([0,1,0])  #vertical lattice vector
    n_atom_in_cell=atomic_basis_xyz.shape[0] #number of atoms in primary unit cell

    #use unit cell to create "supercell" (many copies of unit cell)
    n_atom=int(n_atom_in_cell)*int(n_cell_x)*int(n_cell_y) #total number of atoms in full supercell
    super_cell_xyz=np.empty((n_atom,3))
    n_atom_placed=0
    for i_cell_x in range(n_cell_x):
        for i_cell_y in range(n_cell_y):
            for i_atom_in_cell in range(n_atom_in_cell):
                super_cell_xyz[n_atom_placed,:]=atomic_basis_xyz[i_atom_in_cell,:] + i_cell_x*a1 + i_cell_y*a2 #shift atomic basis vectors by multiple of lattice vector
                n_atom_placed+=1 #increase row counter

    box_vector_x=n_cell_x*a1
    box_vector_y=n_cell_y*a2
    box_vector_z=z_box_len*np.array([0.,0.,1.])
    box_vectors=np.array([box_vector_x,
                          box_vector_y,
                          box_vector_z])
    return super_cell_xyz, box_vectors


def add_mass_lammps_info(xyz,mass,molecule_tag,atom_type,q):
    """
    Prepends mass information and "LAMMPS info" (molecule tag, atom type,
    charge) columns to left of xyz data and returns this new object
    ---Inputs---
    *NOTE: all input properties will be assigned to every atom in xyz*
    xyz: coordinates of all atoms, numpy array, shape(xyz=(n_atoms, 3)
    mass: mass of atoms in xyz, float
    molecule_tag: tag specifying molecule to which atom belongs, integer
    atom_type: number associated with atom type, integer
    q: charge of atom, float
    ---Outputs---
    mass_lammps_info_xyz: mass and "LAMMPS data columns prepended to
                          the left of xyz
    """
    #mass [LAMMPS meta data] [coordinate data] ---equivalent to-->
    #mass mol_tag atom_type q x y z
    mass_lammps_info=np.array([mass,molecule_tag,atom_type,q])
    mass_lammps_info_columns=np.vstack(tuple([mass_lammps_info]*xyz.shape[0]))
    mass_lammps_info_xyz=np.hstack((mass_lammps_info_columns,xyz))
    return mass_lammps_info_xyz
    

def mlixyz_to_lammps_data(file_name,mass_lammps_info_xyz,box_vectors):
    """
    mass, LAMMPS info, xyz -> mlixyz
    Takes numpy array with columns [mass "LAMMPS info" x y z] and writes
    a LAMMPS data file.
    ---Inputs---
    file_name: LAMMPS data file to write, string
    mass_lammps_info_xyz: mass and "LAMMPS info" columns prepended to the left of xyz
    box_vectors: 3x3 array with each row as a box vector (x y z), numpy array
    ---Outputs---
    NONE, the information is written to file defined by file_name
    """

    masses=mass_lammps_info_xyz[:,0] #peel off mass column
    lammps_info_xyz=mass_lammps_info_xyz[:,1:] 
    user_atom_type_list,indices_u=np.unique(lammps_info_xyz[:,1],return_index=True) #get unique atom types from lammps_info_xyz
    mass_list=masses[indices_u] #use unique atom type indices to extract unique masses
    
    n_atoms=lammps_info_xyz.shape[0] #extract total number of atoms
    n_atom_types=user_atom_type_list.shape[0]

    #remap user defined atom types onto 1...N=n_atom_types
    atom_type_list_1_N=range(1,n_atom_types+1)
    for i_atom in range(n_atoms):
        user_atom_type=lammps_info_xyz[i_atom,1]
        lammps_info_xyz[i_atom,1]=np.where(user_atom_type_list==user_atom_type)[0]+1 #get index of user atom type in sorted unique list and write this (+1) as the atom type of that atom

    file_object=open(file_name,'w+')
    #header and cumulative data
    file_object.write('LAMMPS Description\n \n')
    file_object.write(str(n_atoms) + ' atoms\n \n')
    file_object.write(str(n_atom_types) + ' atom types\n \n')

    #box size
    dim_labels=['x','y','z']
    for i_dim,dim_label in enumerate(dim_labels):
        cur_box_vec=box_vectors[i_dim,:]
        min_max=[str(np.min(cur_box_vec)),str(np.max(cur_box_vec))]
        extent_labels=[dim_label+'lo',dim_label+'hi\n']
        file_object.write(' '.join(min_max)+' '+' '.join(extent_labels))

    #atom type - mass mapping
    file_object.write('\nMasses\n\n')
    for i_type,cur_type in enumerate(atom_type_list_1_N):
        cur_mass=mass_list[i_type]
        file_object.write(str(int(cur_type)) + ' ' + str(cur_mass)+'\n')

    #information for each atom
    file_object.write('\nAtoms\n\n')
    for i_row,cur_row in enumerate(lammps_info_xyz):
        data_string=[str(int(i_row+1))] #global ID
        data_string+=[str(int(elem)) for elem in cur_row[0:3]] #mol_tag type q
        data_string+=[f"{elem:.5f}" for elem in cur_row[3:]] #x y z
        file_object.write(' '.join(data_string)+'\n')
        
    file_object.close()
    print('\nWrote LAMMPS data file:',file_name)
    print('Number of atoms:',n_atoms,'\n')


def return_to_box(xyz,origin=None):
    """
    PROBABLY WON'T END UP NEEDING TO USE THIS FUNCTION
    Return atoms outside of box vectors to domain inside that specified by
    the box vectors.
    ---Inputs---
    xyz: coordinates of all atoms, numpy array, shape(xyz=(n_atoms, 3)
    origin: optional value to set origin (where tails of box vectors meet), 3 element 1D numpy array
    ---Outputs---
    meta_xyz: metadata columns prepended to the left of xyz
    """
    #loop over atoms
    #for each atom while loop on each coordinate until it is within bounds
    pass


def monolayer_lammps_data(file_name,n_cell_x,n_cell_y,a_nn,z_box_len,lat_con=None):
    """
    Generates a LAMMPS data file containing a monolayer graphene lattice
    of the specified size and geometry.
    ---Inputs---
    n_cell_x: number of unit cells in x direction, integer
    n_cell_y: number of unit cells in y direction, integer
    a_nn: distance between nearest neighbors, float [Angstroms]
    z_box_len: magnitude of unit cell in z direction, z_box_len/2 gives the
               z-coordinate of the lattice, float [Angstroms]
    lat_con: lattice constant of honeycomb (overrides a_nn,
                      if present), float
    ---Outputs---
    NONE, the information is written to file defined by file_name
    """
    #set contstants
    mass=12.01 #mass of carbon
    molecule_tag=1
    atom_type=1
    q=0 #charge
    #create geometry of lattice, positions and box vectors only
    xyz,box_vector_array=make_xyz_honeycomb_rectangular(n_cell_x,n_cell_y,a_nn,z_box_len,lat_con)
    #add indentity information to atoms in lattice (mass, charge, etc.)
    mass_lammps_info_xyz=add_mass_lammps_info(xyz,mass,molecule_tag,atom_type,q)
    #use above arrays to write LAMMPS data file
    mlixyz_to_lammps_data(file_name,mass_lammps_info_xyz,box_vector_array)


def bilayer_aa_lammps_data(file_name,n_cell_x,n_cell_y,a_nn,interlayer_sep,z_box_len,lat_con=None,mass_list=None):
    """
    Generates a LAMMPS data file containing a bilayer AA graphene lattice
    of the specified size and geometry.
    ---Inputs---
    n_cell_x: number of unit cells in x direction, integer
    n_cell_y: number of unit cells in y direction, integer
    a_nn: distance between nearest neighbors, float [Angstroms]
    interlayer_sep: interlayer separation of graphene sheets, float [Angstroms]
    z_box_len: magnitude of unit cell in z direction, z_box_len/2 gives the
               z-coordinate of the lattice, float [Angstroms]
    lat_con: lattice constant of honeycomb (overrides a_nn,
                      if present), float
    mass_lower_upper: overrides having the same mass for both layers with
                      list [m_lower,m_upper] (only use when file needs to
                      be run with LATTE, as LATTE only identifies atoms
                      by mass), 2 element list
    ---Outputs---
    NONE, the information is written to file defined by file_name
    """
    #optionally assign different masses for layers
    if (mass_list):
        mass_lower=mass_list[0]
        mass_upper=mass_list[1]
    else:
        mass_lower=12.01
        mass_upper=mass_lower
            
    #set constants (per layer, as necessary)
    molecule_tag_lower=1
    molecule_tag_upper=2
    atom_type_lower=1
    atom_type_upper=2
    q=0
    
    #create geometry of lattice, positions and box vectors only
    #lower layer
    xyz_lower,box_vector_array=make_xyz_honeycomb_rectangular(n_cell_x,n_cell_y,a_nn,z_box_len,lat_con)
    interlayer_sep_over_2_shift=np.vstack(tuple([np.array([0.,0.,interlayer_sep/2])]*xyz_lower.shape[0])) #create array of same shape as xyz_lower to shift all atoms up by interlayer_sep/2
    xyz_lower-=interlayer_sep_over_2_shift #shift atoms in lower layer from center of box down to final height
    #upper layer
    xyz_upper=copy.deepcopy(xyz_lower)
    xyz_upper+=2*interlayer_sep_over_2_shift #shift xyz_upper from lower layer height to final upper layer height

    #add mass and "LAMMPS info" to xyz
    mlixyz_lower=add_mass_lammps_info(xyz_lower,mass_lower,molecule_tag_lower,atom_type_lower,q)
    mlixyz_upper=add_mass_lammps_info(xyz_upper,mass_upper,molecule_tag_upper,atom_type_upper,q)

    #append upper layer to lower layer
    mlixyz=np.vstack((mlixyz_lower,mlixyz_upper))

    #write LAMMPS data file
    mlixyz_to_lammps_data(file_name,mlixyz,box_vector_array)


def bilayer_ab_lammps_data(file_name,n_cell_x,n_cell_y,a_nn,interlayer_sep,z_box_len,lat_con=None,mass_list=None):
    """
    Generates a LAMMPS data file containing a bilayer AB graphene lattice
    of the specified size and geometry.
    ---Inputs---
    n_cell_x: number of unit cells in x direction, integer
    n_cell_y: number of unit cells in y direction, integer
    a_nn: distance between nearest neighbors, float [Angstroms]
    interlayer_sep: interlayer separation of graphene sheets, float [Angstroms]
    z_box_len: magnitude of unit cell in z direction, z_box_len/2 gives the
               z-coordinate of the lattice, float [Angstroms]
    lat_con: lattice constant of honeycomb (overrides a_nn,
                      if present), float
    mass_lower_upper: overrides having the same mass for both layers with
                      list [m_lower,m_upper] (only use when file needs to
                      be run with LATTE, as LATTE only identifies atoms
                      by mass), 2 element list
    ---Outputs---
    NONE, the information is written to file defined by file_name
    """
    #optionally assign different masses for layers
    if (mass_list):
        mass_lower=mass_list[0]
        mass_upper=mass_list[1]
    else:
        mass_lower=12.01
        mass_upper=mass_lower
            
    #set constants (per layer, as necessary)
    molecule_tag_lower=1
    molecule_tag_upper=2
    atom_type_lower=1
    atom_type_upper=2
    q=0
    
    #create geometry of lattice, positions and box vectors only
    #lower layer
    xyz_lower,box_vector_array=make_xyz_honeycomb_rectangular(n_cell_x,n_cell_y,a_nn,z_box_len,lat_con)
    interlayer_sep_over_2_shift=np.vstack(tuple([np.array([0.,0.,interlayer_sep/2])]*xyz_lower.shape[0])) #create array of same shape as xyz_lower to shift all atoms up by interlayer_sep/2
    a_nn_used=np.linalg.norm(xyz_lower[1]-xyz_lower[0]) #extract the nearest neighbor distance from lower layer, since one can't tell this distance from inputs to this function (possible override by lat_con INSIDE make_xyz_honeycomb_rectangular) 
    ab_horizontal_shift=np.vstack(tuple([np.array([a_nn_used,0.,0.])]*xyz_lower.shape[0])) #create array of same shape as xyz_lower to shift all atoms up by interlayer_sep/2
    xyz_lower-=interlayer_sep_over_2_shift #shift atoms in lower layer from center of box down to final height
    #upper layer
    xyz_upper=copy.deepcopy(xyz_lower)
    xyz_upper+=2*interlayer_sep_over_2_shift+ab_horizontal_shift #shift xyz_upper from lower layer height to final upper layer height and shift horizontally to make AB from AA

    #add mass and "LAMMPS info" to xyz
    mlixyz_lower=add_mass_lammps_info(xyz_lower,mass_lower,molecule_tag_lower,atom_type_lower,q)
    mlixyz_upper=add_mass_lammps_info(xyz_upper,mass_upper,molecule_tag_upper,atom_type_upper,q)

    #append upper layer to lower layer
    mlixyz=np.vstack((mlixyz_lower,mlixyz_upper))

    #write LAMMPS data file
    mlixyz_to_lammps_data(file_name,mlixyz,box_vector_array)



if (__name__=="__main__"):
    #def monolayer_graphene_lammps_data(file_name,n_cell_x,n_cell_y,a_nn,z_box_len,lat_con=None):
    monolayer_lammps_data('test_mono.data',2,1,0.,10,lat_con=1.0)

    #def bilayer_aa_lammps_data(file_name,n_cell_x,n_cell_y,a_nn,interlayer_sep,z_box_len,lat_con=None,mass_list=None):
    bilayer_aa_lammps_data('test_aa.data',2,2,1.0,1,10,mass_list=[12.0101,12.01])

    #def bilayer_ab_lammps_data(file_name,n_cell_x,n_cell_y,a_nn,interlayer_sep,z_box_len,lat_con=None,mass_list=None):
    bilayer_ab_lammps_data('test_ab.data',2,2,1.0,1,10,mass_list=[12.0101,12.01])
