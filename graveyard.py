#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 17:16:20 2024

@author: benjaminlear
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from plotly import graph_objects as go
import numpy as np

# Define the SMILES string
smiles = 'CC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)N'  # Ethanol as an example

# Convert the SMILES string to a molecule
mol = Chem.MolFromSmiles(smiles)

# Add hydrogens to the molecule
mol = Chem.AddHs(mol)

# Embed the molecule in 3D space
AllChem.EmbedMolecule(mol, randomSeed=42)

# Optimize the 3D coordinates using force field optimization
AllChem.UFFOptimizeMolecule(mol)

# Get the 3D coordinates
conf = mol.GetConformer()

atoms_colors = dict(
    C = "darkgrey",
    H = "lightgrey",
    N = "blue",
    O = "red",
    )
atoms_sizes = dict( # these are just van der waals radii
    H = 1.2,
    C = 1.7,
    N = 1.55,
    O = 1.52,
    )





















    


    
def draw_mol_3D (mol, resolution = 100):
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()
    
    fig = go.Figure()
    
    #add atoms
    for i, atom in enumerate(atoms):
        
        asym = atom.GetSymbol()
        atom_trace = make_atom_mesh_trace(
           [
            conf.GetAtomPosition(i).x,
            conf.GetAtomPosition(i).y,
            conf.GetAtomPosition(i).z
            ], 
           color =  atoms_colors[asym],
           radius = atoms_sizes[atom.GetSymbol()]*0.2, # change this for ball and stick, stick, van der waals, etc. 
           resolution = resolution,
            )
        fig.add_trace(atom_trace)
        
        
       
    #add bonds
    for bond in bonds:
        atom1 = bond.GetBeginAtomIdx()
        atom2 = bond.GetEndAtomIdx()
        
        a1x = conf.GetAtomPosition(atom1).x
        a1y = conf.GetAtomPosition(atom1).y
        a1z = conf.GetAtomPosition(atom1).z
        
        a2x = conf.GetAtomPosition(atom2).x
        a2y = conf.GetAtomPosition(atom2).y
        a2z = conf.GetAtomPosition(atom2).z
        
        midx = (a1x + a2x)/2
        midy = (a1y + a2y)/2
        midz = (a1z + a2z)/2
        
        bond_trace = make_bond_mesh_trace(
                                     [a1x, a1y, a1z], 
                                     [midx, midy, midz],
                                     color = atoms_colors[atoms[atom1].GetSymbol()],
                                     resolution = resolution)
        fig.add_trace(bond_trace)
        
        bond_trace = make_bond_mesh_trace(
                                     [midx, midy, midz],
                                     [a2x, a2y, a2z],  
                                     color = atoms_colors[atoms[atom2].GetSymbol()],
                                     resolution = resolution)
        fig.add_trace(bond_trace)
        
    fig.update_layout(
    scene=dict(
        xaxis=dict(
            visible=False,
            showbackground=False,
            showgrid=False,
            zeroline=False
        ),
        yaxis=dict(
            visible=False,
            showbackground=False,
            showgrid=False,
            zeroline=False
        ),
        zaxis=dict(
            visible=False,
            showbackground=False,
            showgrid=False,
            zeroline=False
        ),
        #aspectmode='data',  # Ensure the aspect ratio is based on the data
    ),
    width=800,
    height=600,
    margin=dict(l=0, r=0, t=0, b=0)  # Reduce margins to focus on the data
)
    fig.show("browser")
    
    return fig

fig = draw_mol_3D(mol, resolution = 64) 



#%% Old functions that are not needed now. 



def generate_surface_cylinder(point1, point2, radius=0.1, resolution=50):
    # Vector from point1 to point2
    v = np.array(point2) - np.array(point1)
    mag = np.linalg.norm(v)
    v = v / mag  # Normalize vector

    # Create orthogonal vectors to v
    not_v = np.array([1, 0, 0]) if (v == np.array([0, 1, 0])).all() else np.array([0, 1, 0])
    n1 = np.cross(v, not_v)
    n1 /= np.linalg.norm(n1)
    n2 = np.cross(v, n1)
    
    # Create the circle in 2D
    theta = np.linspace(0, 2 * np.pi, resolution)
    circle_x = radius * np.cos(theta)
    circle_y = radius * np.sin(theta)

    # Extrude circle along the vector v
    t = np.linspace(0, mag, resolution)
    t_grid, theta_grid = np.meshgrid(t, theta)
    X = point1[0] + t_grid * v[0] + radius * np.outer(np.cos(theta), n1[0]) + radius * np.outer(np.sin(theta), n2[0])
    Y = point1[1] + t_grid * v[1] + radius * np.outer(np.cos(theta), n1[1]) + radius * np.outer(np.sin(theta), n2[1])
    Z = point1[2] + t_grid * v[2] + radius * np.outer(np.cos(theta), n1[2]) + radius * np.outer(np.sin(theta), n2[2])

    return X, Y, Z



def generate_mesh_cylinder(point1, point2, radius=0.1, resolution=32):
    # Ensure point1 and point2 are numpy arrays
    point1 = np.array(point1)
    point2 = np.array(point2)

    # Calculate the vector from point1 to point2
    v = point2 - point1
    height = np.linalg.norm(v)
    v = v / height  # Normalize the vector

    # Create orthogonal vectors to the cylinder axis
    if np.all(v == np.array([0, 0, 1])) or np.all(v == np.array([0, 0, -1])):
        # Special case where the axis is along the z direction
        not_v = np.array([1, 0, 0])
    else:
        not_v = np.array([0, 0, 1])

    n1 = np.cross(v, not_v)
    n1 /= np.linalg.norm(n1)
    n2 = np.cross(v, n1)

    # Generate the angles for the circular cross-section
    theta = np.linspace(0, 2 * np.pi, resolution, endpoint=True)
    circle = np.array([np.cos(theta), np.sin(theta)]) * radius

    # Generate the points for the bottom and top circles of the cylinder
    bottom_circle = point1[:, None] + n1[:, None] * circle[0] + n2[:, None] * circle[1]
    top_circle = bottom_circle + v[:, None] * height

    # Stack the bottom and top circles
    points = np.hstack([bottom_circle, top_circle])

    x, y, z = points

    return x, y, z



def generate_surface_sphere(center, radius=0.1, resolution=50):
    u = np.linspace(0, 2 * np.pi, resolution)
    v = np.linspace(0, np.pi, resolution)
    x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
    y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
    z = center[2] + radius * np.outer(np.ones(np.size(u)), np.cos(v))
    return x, y, z


def generate_mesh_sphere (center, radius = 0.1, resolution = 32):
    d = np.pi/resolution # sets resolution

    theta, phi = np.mgrid[0:np.pi+d:d, 0:2*np.pi:d]
    # Convert to Cartesian coordinates
    x = np.sin(theta) * np.cos(phi) * radius + center[0]
    y = np.sin(theta) * np.sin(phi) * radius + center[1]
    z = np.cos(theta) * radius + center[2]
    # print(x.shape, y.shape, z.shape)  # (33, 64) (33, 64) (33, 64)
    points = np.vstack([x.ravel(), y.ravel(), z.ravel()])
    # print(points.shape)  # (3, 2112)
    x, y, z = points
    return x, y, z


def make_bond_surface_trace(fig, point1, point2, radius = 0.1, resolution = 50, color = "grey"):
    x, y, z = generate_surface_cylinder(point1, point2, radius = radius, resolution = resolution)
    bond_trace = go.Surface(
        x=x, 
        y=y, 
        z=z, 
        surfacecolor=np.full_like(x, 0),  # Assign a single value to make the sphere a solid color
        colorscale=[[0, color], [1, color]],  # Define the colorscale as a single color
        lighting=dict(
            ambient=0,
            diffuse=1,
            specular=0,
            roughness=1,
            fresnel=0
        ),
        lightposition=dict(
            x=1000,
            y=1000,
            z=1000
        ),
        showscale=False  # Hide the individual color scales
        )   
    return bond_trace


def make_atom_surface_trace(fig, center, radius = 0.1, resolution = 50, color = "grey" ):
    x, y, z = generate_surface_sphere(center, radius = radius, resolution = resolution)
    atom_trace = go.Surface(
        x=x, 
        y=y, 
        z=z, 
        surfacecolor=np.full_like(x, 0),  # Assign a single value to make the sphere a solid color
        colorscale=[[0, color], [1, color]],  # Define the colorscale as a single color
        lighting=dict(
            ambient=0,
            diffuse=1,
            specular=0,
            roughness=1,
            fresnel=0
        ),
        lightposition=dict(
            x=1000,
            y=1000,
            z=1000
        ),
        showscale=False  # Hide the individual color scales
    )
    return atom_trace
