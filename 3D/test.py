#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 11:41:46 2024

@author: benjaminlear
"""
from plotlyMol3D import draw_3D_mol

#%% smiles test
moldraw = draw_3D_rep(smiles = "CCNCOCSC", 
            mode = "ball+stick", 
            ambient = 0.1
            ) 

#%% xyztest

moldraw = draw_3D_rep(xyzfile = "/Users/benjaminlear/Downloads/test.xyz", 
            mode = "ball+stick", 
            ambient = 0.1
            ) 



#%% Cubetest

moldraw = draw_3D_rep(
            cubefile = "/Users/benjaminlear/Documents/GitHub/plotlyMol/3D/anto_occ_1-min2.cube", 
            molfile = "/Users/benjaminlear/Documents/GitHub/plotlyMol/3D/cube.mol",
            mode = "ball+stick", 
            ambient = 0.1,
            cubedraw = "orbitals",
            orbital_opacity = 0.25,
            orbital_colors = ["darkorange", "darkblue"],
            ) 
#%% multidraw test...

moldraw = draw_3D_rep(smiles = "CCNCOCSC", xyzfile = "/Users/benjaminlear/Downloads/test.xyz", 
            mode = "ball+stick", 
            ambient = 0.1
            ) 


#%%
m = Chem.MolFromMolFile('/Users/benjaminlear/Documents/GitHub/plotlyMol/3D/cube.mol')


#%%

testblock = 1
testblock = cubefile_to_xyzblock("/Users/benjaminlear/Documents/GitHub/plotlyMol/3D/anto_occ_1-min2.cube")


    


# def process_xyz_coords(xyz, charge = 0): # should be a block of coords, expected by rdkit
#     raw_mol = Chem.MolFromXYZBlock(ind)
#     from rdkit.Chem import rdDetermineBonds <-- hangs on import
#     conn_mol = Chem.Mol(raw_mol)
    
#     rdDetermineBonds.DetermineConnectivity(conn_mol)
#     rdDetermineBonds.DetermineBondOrders(conn_mol, charge=charge)
    
#     atoms = conn_mol.GetAtoms()
#     bonds = conn_mol.GetBonds()
    
#     return atoms, bonds


#raw_mol = Chem.MolFromXYZFile('/Users/benjaminlear/Downloads/test.xyz')

# testblock = f'''9
# 	Energy:       4.3719968
# C         -5.60141        3.84674       -0.07080
# C         -4.10301        3.83876       -0.08074
# H         -3.72298        4.87682        0.04590
# H         -3.75288        3.44948       -1.06129
# N         -3.59011        2.98389        0.98445
# H         -3.85276        3.38964        1.91098
# H         -2.54705        2.97307        0.92691
# O         -6.21202        4.25745        0.90282
# H         -6.15086        3.49213       -0.93685
#                                '''


