import ase
import flatgraphene as fg

#example to modify when working on module
"""
atoms = fg.shift.make_graphene(stacking=['A','B','C'],cell_type='hex',
                               n_layer=3,n_1=3,n_2=3,lat_con=0.0,a_nn=1.5,
                               sep=[3.0,2,4],sym=['O','F','N'],mass=[0,1,2])
ase.visualize.view(atoms)
"""

#p_found, q_found, theta_comp = find_p_q(21.79)
p_found, q_found, theta_comp = fg.twist.find_p_q(9.43)
print('generating system with a {:.2f} degree twist'.format(theta_comp))
test_sheet = fg.twist.make_graphene(cell_type='hex',n_layer=2,
                                    p=p_found,q=q_found,lat_con=0.0,a_nn=1.5,
                                    sep=3.4,sym=['C','Ne'],mass=[12.01,12.02],
                                    mol_id=0,h_vac=3)
ase.visualize.view(test_sheet)
