


"""
Written by Tawfiqur Rakib. Awaiting integration into flat_graphene as a function.
"""

# CSL Anges= 21.7867893, 13.17355112, 9.43000791, 6.0089832, 7.34099302, 5.08584781,4.408455008,3.89023817,3.4...,2.875894634,2.0046278307, 1.4701297259, 1.1606527199, 1.08454904924, 1.0501208798, 0.98743029788, 0.1000899404333
# 0.5050308328


import numpy as np

def POSCAR_writer(filename, atom_num, a1, a2, b2, c3, xyz):
    f = open(filename, "w")
    f.write('written by TR \n')
    f.write('%.9f \n'%(1.0))
    f.write('%.9f %.9f %.9f \n'%(a1, 0.0, 0.0))
    f.write('%.9f %.9f %.9f \n'%(a2,b2, 0.0))
    f.write('%.9f %.9f %.9f \n'%(0.0, 0.0, c3))
    f.write('Type1 \n')
    f.write('%d\n'%(atom_num))
    f.write('Cartesian \n')

    for i in range (atom_num):
        f.write('%1.8f %1.8f %1.8f \n'%(xyz[i,0],xyz[i,1],xyz[i,2]))

    
    f.close()

def RotMatrix(theta_rad):

	return np.array(([np.cos(theta_rad), -np.sin(theta_rad),0],[np.sin(theta_rad), np.cos(theta_rad), 0],[0,0,1]))



if __name__=="__main__":
    t = 21.7867893

    filename = "structure_lammps.txt"
    
    a=2.456
    
    p=1  
    
    
    theta= (np.pi/180)*t
    
    
    
        
    
    latvec1a = a*np.array([1/np.sqrt(2),1/np.sqrt(2)])
    
    latvec1b = RotMatrix(np.pi/3)[:2,:2].dot(latvec1a)
    
    lv1 = np.array([latvec1a,latvec1b])
    
    
    
    latvec2a = p*RotMatrix(theta)[:2,:2].dot(latvec1a)
    
    latvec2b = RotMatrix(np.pi/3)[:2,:2].dot(latvec2a)
    
    lv2 = np.array([latvec2a,latvec2b])
    
    
    
    	
    
    M = lv2.dot(np.linalg.inv(lv1))
    
    
    
    P = np.linalg.inv(np.identity(len(M))-M).dot(M)
    
    
    
    lvm = P.dot(lv1)
    
    
    
    m=lvm.dot(np.linalg.inv(lv1))
    
            
    
    print('m = %f \n' %-m[0,0]) 
    
    print('n = %f \n' %m[1,1]) 
    
    
    
    
    
    ss = 1
    
    lvmp=np.array([[lvm[0,0],lvm[0,1]],[(2*lvm[1,0])-lvm[0,0],2*(lvm[1,1])-lvm[0,1]]])
    
    xmod=np.sqrt(np.square(lvm[0,0]) + np.square(lvm[0,1]) )
    
    
    
    xaxis=[1,0]
    
    lvmpx=np.array([lvm[0,0],lvm[0,1]])
    
    
    
    dotprod=lvmpx.dot(xaxis)/xmod
    
    phi=np.arccos(lvmpx.dot(xaxis)/xmod)
    
    phi0=180*phi/np.pi
    
    
    
    np.linalg.norm(lvm,axis=1)
    
    
    
    mag = np.linalg.norm(lvm,axis=1).mean()
    
    
    
    lvmp=np.array([[lvm[0,0],lvm[0,1]],[(2*lvm[1,0])-lvm[0,0],2*(lvm[1,1])-lvm[0,1]]])
    
    rotlvmp1= RotMatrix(-phi)[:2,:2].dot(lvmp[0])
    
    rotlvmp2= RotMatrix(-phi)[:2,:2].dot(lvmp[1]) 
    
    rotlvmp=[rotlvmp1,rotlvmp2] 
    
    
    
    #expand lattice of lower layer
    
    latt_x = int(1000)
    
    latt_y = int(1000)
    
    tot_mat = int(400000)
    
    rot_num = latt_x*latt_y*2
    
    half_x = int(latt_x/2)
    
    half_y = int(latt_y/2)
    
    lattice_points1= np.ndarray(shape=(latt_x,latt_y,3))
    
    
    
    for j in range(-half_x,half_y):
    
        for i in range(-half_x,half_y):
    
                lattice_points1[i][j][0] = i*latvec1a[0] + j*latvec1b[0]
    
                lattice_points1[i][j][1] = i*latvec1a[1] + j*latvec1b[1]
    
                lattice_points1[i][j][2] = 0
    
    
    
    #Add atom positions to lattice
    
    
    
    
    
    rb0=[0,0]
    
    #rb1=[0,0]
    
    rb1 = rb0 + (latvec1a + latvec1b)/3
    
    
    
    atom_positions1=np.ndarray(shape=(latt_x,latt_y,2,3))
    
    for j in range(-half_x,half_y):
    
        for i in range(-half_x,half_y):
    
            atom_positions1[i][j][0][0]=lattice_points1[i][j][0] +rb0[0]
    
            atom_positions1[i][j][0][1]=lattice_points1[i][j][1] +rb0[1]
    
            atom_positions1[i][j][0][2]=0
    
            atom_positions1[i][j][1][0]=lattice_points1[i][j][0] +rb1[0]
    
            atom_positions1[i][j][1][1]=lattice_points1[i][j][1] +rb1[1]
    
            atom_positions1[i][j][1][2]=0
    
       
    
            
    
    #expand lattice of over layer
    
    lattice_points2= np.ndarray(shape=(latt_x,latt_y,3))
    
    
    
    for j in range(-half_x,half_y):
    
        for i in range(-half_x,half_y):
    
                lattice_points2[i][j][0] = i*latvec2a[0] + j*latvec2b[0]
    
                lattice_points2[i][j][1] = i*latvec2a[1] + j*latvec2b[1]
    
                lattice_points2[i][j][2] = 0
    
    
    
    #Add atom positions to lattice
    
    
    
    
    
    ro0=[0,0]
    
    #ro1=[0,0]
    
    ro1 = ro0 + (latvec2a + latvec2b)/3
    
    
    
    atom_positions2=np.ndarray(shape=(latt_x,latt_y,2,3))
    
    for j in range(-half_x,half_y):
    
        for i in range(-half_x,half_y):
    
            atom_positions2[i][j][0][0]=lattice_points2[i][j][0] +ro0[0]
    
            atom_positions2[i][j][0][1]=lattice_points2[i][j][1] +ro0[1]
    
            atom_positions2[i][j][0][2]=2.5
    
            atom_positions2[i][j][1][0]=lattice_points2[i][j][0] +ro1[0]
    
            atom_positions2[i][j][1][1]=lattice_points2[i][j][1] +ro1[1]
    
            atom_positions2[i][j][1][2]=2.5
    
    
    
    
    
    
    
    
    
    coor=[atom_positions1 + atom_positions2]
    
    #coordinates1=[[0]]*800
    
    #for i in range(800):
    
    #    coordinates1[i]=[0]*3
    
    coordinates1=np.ndarray(shape=(rot_num,3))    
    
    b=0
    
    while b < rot_num :
    
        for j in range(-half_x,half_y):
    
            for i in range(-half_x,half_y):
    
                coordinates1[b][0]=atom_positions1[i][j][0][0]
    
                coordinates1[b][1]=atom_positions1[i][j][0][1]
    
                coordinates1[b][2]=atom_positions1[i][j][0][2]
    
                coordinates1[b+1][0]=atom_positions1[i][j][1][0]
    
                coordinates1[b+1][1]=atom_positions1[i][j][1][1]
    
                coordinates1[b+1][2]=atom_positions1[i][j][1][2]
    
                b+=2
    
    
    
    rotcoordinates1=np.ndarray(shape=(rot_num ,3))
    
    
    
    for b in range(0,rot_num):
    
        rotcoordinates1[b,:2]=RotMatrix(-phi)[:2,:2].dot(coordinates1[b,:2])
    
        
    
    #supercell of rotated moire pattern upper
    
    fincoord1=np.ndarray(shape=(tot_mat,3))
    
    fincoord1.fill(0)
    
    
    
    v=0
    
    for c in range(0,rot_num):
    
        if (rotcoordinates1[c,0]>-0.0001 and rotcoordinates1[c,0]<ss*rotlvmp1[0] and rotcoordinates1[c,1]>-0.0001 and rotcoordinates1[c,1]<ss*rotlvmp2[1]):
    
            fincoord1[v,0]=rotcoordinates1[c,0]
    
            fincoord1[v,1]=rotcoordinates1[c,1]
    
            fincoord1[v,2]=rotcoordinates1[c,2]
    
            v+=1
    
           
    
    
    
    coordinates2=np.ndarray(shape=(rot_num,3))    
    
    b=0
    
    while b < rot_num :
    
        for j in range(-half_x,half_y):
    
            for i in range(-half_x,half_y):
    
                coordinates2[b][0]=atom_positions2[i][j][0][0]
    
                coordinates2[b][1]=atom_positions2[i][j][0][1]
    
                coordinates2[b][2]=atom_positions2[i][j][0][2]
    
                coordinates2[b+1][0]=atom_positions2[i][j][1][0]
    
                coordinates2[b+1][1]=atom_positions2[i][j][1][1]
    
                coordinates2[b+1][2]=atom_positions2[i][j][1][2]
    
                b+=2
    
    
    
    #rotated moire pattern lower
    
    rotcoordinates2=np.ndarray(shape=(rot_num ,3))
    
    
    
    for b in range(0,rot_num):
    
        rotcoordinates2[b,:2]=RotMatrix(-phi)[:2,:2].dot(coordinates2[b,:2])
    
    
    
    #supercell of rotated moire pattern upper
    
    fincoord2=np.ndarray(shape=(tot_mat,3))
    
    fincoord2.fill(0)
    
    v=0
    
    for a in range(0,rot_num):
    
        if (rotcoordinates2[a,0]>-0.0001 and rotcoordinates2[a,0]<ss*rotlvmp1[0] and rotcoordinates2[a,1]>-0.0001 and rotcoordinates2[a,1]<ss*rotlvmp2[1]):
    
            fincoord2[v,0]=rotcoordinates2[a,0]
    
            fincoord2[v,1]=rotcoordinates2[a,1]
    
            fincoord2[v,2]=3.4
    
            v+=1
    
            
    
                    
    
    
    
    ntype1=1
    
    ntype2=2    
    
    xlim=[0,rotlvmp1[0]]
    
    ylim=[0,rotlvmp2[1]]
    
    zlim=[0,18.4]
    
    
    p = 0
    for i in range(len(fincoord1)-1):
        if (fincoord1[i,0] != 0 or fincoord1[i+1,0] != 0):
            p = p+1
    
    for i in range(len(fincoord2)-1):
        if (fincoord2[i,0] != 0 or fincoord2[i+1,0] != 0):
            p=p+1
      
    
    natom = int(p)
    
    #write to files  
    
    
    
    f=open(filename,'w+')
    
    print("writing")
    
    f.write('(written by TR)\n\n')
    
    f.write('%d\tatoms\n'%(p))
    
    f.write('2\tatom types\n')
    
    f.write('%.9g\t%.9g\txlo xhi\n'%(xlim[0]*ss,xlim[1]*ss))
    
    f.write('%.9g\t%.9g\tylo yhi\n'%(ylim[0]*ss,ylim[1]*ss))
    
    f.write('%.9g\t%.9g\tzlo zhi\n'%(zlim[0],zlim[1]))
    
    
    
    f.write('\n\n')
    
    f.write('Masses\n\n')
    
    
    
    f.write('1\t12\n')
    
    f.write('2\t12\n\n')
    
    f.write('Atoms\n\n')
    
    p = 1
    
    for i in range(int(natom/2)):
    
    
    
        f.write('%d %d %d %d %1.8f %1.8f %1.8f %d %d %d\n'%(p,ntype1,ntype1,0,fincoord1[i,0],fincoord1[i,1],7.5,0,0,0))
    
        p = p+1
    
    for i in range(int(natom/2)):
    
    
        f.write('%d %d %d %d %1.8f %1.8f %1.8f %d %d %d\n'%(p,ntype2,ntype2,0,fincoord2[i,0],fincoord2[i,1],10.9,0,0,0))
    
        p=p+1
    
    f.close()
    
    
