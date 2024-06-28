from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np

from dataclasses import dataclass, field
from typing import List

from rdkit import Chem
from rdkit.Chem import AllChem
from atomProperties import *

DEFAULT_RESOLUTION = 32
DEFAULT_RADIUS = 0.1

@dataclass
class Atom:
    atom_id: int
    atom_number: int = field(default=0)
    atom_symbol: str = field(default="unknown")
    atom_xyz: List[float] = field(default_factory=list)
    atom_vdw: float = field(default=1.70)  # default value

@dataclass
class Bond:
    a1_id: int
    a2_id: int 
    a1_number: int
    a2_number: int
    a1_xyz: List[float] = field(default_factory=list)
    a2_xyz: List[float] = field(default_factory=list)
    a1_vdw: float = field(default=1.70)
    a2_vdw: float = field(default=1.70)



###
# PROCESSING INPUTS
###

# def process_xyz_coords(xyz, charge = 0): # should be a block of coords, expected by rdkit
#     raw_mol = Chem.MolFromXYZBlock(ind)
#     from rdkit.Chem import rdDetermineBonds <-- hangs on import
#     conn_mol = Chem.Mol(raw_mol)
    
#     rdDetermineBonds.DetermineConnectivity(conn_mol)
#     rdDetermineBonds.DetermineBondOrders(conn_mol, charge=charge)
    
#     atoms = conn_mol.GetAtoms()
#     bonds = conn_mol.GetBonds()
    
#     return atoms, bonds


def process_smiles(smiles):
    # Convert the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)

    # Add hydrogens to the molecule
    mol = Chem.AddHs(mol)
    
    # Get atom and bond information
    atoms = mol.GetAtoms()
    bonds = mol.GetBonds()

    # Embed the molecule in 3D space
    AllChem.EmbedMolecule(mol, randomSeed=42)

    # Optimize the 3D coordinates using force field optimization
    AllChem.UFFOptimizeMolecule(mol)

    # Get the 3D coordinates
    conf = mol.GetConformer()
    
    
    # this could be its own function, that returns the atomList and bondList, working from rdkit's atoms bonds
    atomList = []
    for a in atoms:
        atomList.append(
            Atom(
                atom_id = a.GetIdx(),
                atom_number = a.GetAtomicNum(), 
                atom_symbol = a.GetSymbol(),
                atom_xyz = [
                    conf.GetAtomPosition(a.GetIdx()).x,
                    conf.GetAtomPosition(a.GetIdx()).y,
                    conf.GetAtomPosition(a.GetIdx()).z
                    ],
                atom_vdw = vdw_radii[a.GetAtomicNum()],
                ))
        
    
    bondList = []
    for b in bonds:
        bondList.append(
            Bond(
                a1_id = b.GetBeginAtomIdx(),
                a2_id = b.GetEndAtomIdx(),
                a1_number = b.GetBeginAtom().GetAtomicNum(),
                a2_number = b.GetEndAtom().GetAtomicNum(),
                a1_xyz = atomList[b.GetBeginAtomIdx()].atom_xyz,
                a2_xyz = atomList[b.GetEndAtomIdx()].atom_xyz,
                a1_vdw = atomList[b.GetBeginAtomIdx()].atom_vdw,
                a2_vdw = atomList[b.GetEndAtomIdx()].atom_vdw,
                ))
    
    return atomList, bondList

###
# ATOM DRAWING STUFF
###
DEFAULT_RADIUS = 0.1
DEFAULT_RESOLUTION = 32
a_res_scale = 10
def make_fibonacci_sphere(center, radius=DEFAULT_RADIUS, resolution = DEFAULT_RESOLUTION):
    
    num_points = resolution
    indices = np.arange(0, num_points, dtype=float) + 0.5
    phi = np.arccos(1 - 2*indices/num_points)
    theta = np.pi * (1 + 5**0.5) * indices
    
    x = radius * np.sin(phi) * np.cos(theta) + center[0]
    y = radius * np.sin(phi) * np.sin(theta) + center[1]
    z = radius * np.cos(phi) + center[2]

    return x, y, z


def make_atom_mesh_trace(atom, radius = DEFAULT_RADIUS, resolution = DEFAULT_RESOLUTION, color = "grey"):
    
    # we are now at the stage of drawing a single trace for a single atom. 
    
    # first, we make sure we are using the correct radius for the atom...
    if radius == "vdw": # check to see if we are going to use the vdw radii
        radius = atom.atom_vdw # if so, then reassign the value
    elif radius == "ball":
        radius = atom.atom_vdw * 0.2

    x, y, z = make_fibonacci_sphere (atom.atom_xyz, radius = radius, resolution = resolution*a_res_scale)
    
    atom_trace = go.Mesh3d(
        x=x, 
        y=y, 
        z=z, 
        color=atom_colors[atom.atom_number], 
        opacity=1, 
        alphahull=0,
        name=f'{atom.atom_symbol}{atom.atom_id}',
        hoverinfo='name',  # Only show the name on hover
        )
    return atom_trace

def draw_atoms (fig, atomList, resolution = DEFAULT_RESOLUTION, radius = DEFAULT_RADIUS):
    for a in atomList: # go through each atom...
        a_trace = make_atom_mesh_trace(a, resolution = resolution, radius = radius)
            
        fig.add_trace(a_trace)
    return fig


####
# DRAWING BONDS STUFF
####

# at each atom position, add a sphere of a given size and color
def generate_cylinder_mesh_rectangles(point1, point2, radius=DEFAULT_RADIUS, resolution=DEFAULT_RESOLUTION):
    point1 = np.array(point1)
    point2 = np.array(point2)
    
    v = point2 - point1
    height = np.linalg.norm(v)
    v = v / height  # Normalize the vector

    # Find two vectors orthogonal to the axis of the cylinder
    if np.allclose(v, np.array([0, 0, 1])) or np.allclose(v, np.array([0, 0, -1])):
        not_v = np.array([1, 0, 0])
    else:
        not_v = np.array([0, 0, 1])
        
    n1 = np.cross(v, not_v)
    n1 /= np.linalg.norm(n1)
    n2 = np.cross(v, n1)
    
    # Generate the angles for the circular cross-section
    theta = np.linspace(0, 2 * np.pi, resolution, endpoint=False)
    circle = np.array([np.cos(theta), np.sin(theta)])
    
    # Generate the points for the bottom and top circles of the cylinder
    bottom_circle = point1[:, None] + radius * (n1[:, None] * circle[0] + n2[:, None] * circle[1])
    top_circle = point2[:, None] + radius * (n1[:, None] * circle[0] + n2[:, None] * circle[1])
    
    x = np.concatenate([bottom_circle[0], top_circle[0]])
    y = np.concatenate([bottom_circle[1], top_circle[1]])
    z = np.concatenate([bottom_circle[2], top_circle[2]])
    
    return x, y, z


def make_bond_mesh_trace(point1, point2, radius = DEFAULT_RADIUS, resolution = DEFAULT_RESOLUTION, color = "grey"):
    
    x, y, z = generate_cylinder_mesh_rectangles(point1, point2, radius, resolution)
    
    # Create the faces for the cylinder using rectangles
    i, j, k, l = [], [], [], []
    num_vertices = resolution
    for n in range(num_vertices):
        next_n = (n + 1) % num_vertices
        i.extend([n, next_n, next_n, n])
        j.extend([n, n, n + num_vertices, n + num_vertices])
        k.extend([n + num_vertices, n + num_vertices, next_n + num_vertices, next_n + num_vertices])
        l.extend([next_n, next_n + num_vertices, next_n + num_vertices, next_n])
    
    bond_trace = go.Mesh3d(
        x=x,
        y=y,
        z=z,
        i=i,
        j=j,
        k=k,
        color=color,
        opacity=1,
        hoverinfo="skip",
    )
    return bond_trace


def draw_bonds(fig, bondList, resolution = DEFAULT_RESOLUTION, radius = DEFAULT_RADIUS):
    for bond in bondList:
        
        # eventually weight by atom size
        midx = (bond.a1_xyz[0] + bond.a2_xyz[0])/2
        midy = (bond.a1_xyz[1] + bond.a2_xyz[1])/2
        midz = (bond.a1_xyz[2] + bond.a2_xyz[2])/2
        
        bond_trace = make_bond_mesh_trace(
                                     bond.a1_xyz, 
                                     [midx, midy, midz],
                                     color = atom_colors[bond.a1_number],
                                     resolution = resolution,
                                     )
        fig.add_trace(bond_trace)
        
        bond_trace = make_bond_mesh_trace(
                                     [midx, midy, midz],
                                     bond.a2_xyz,  
                                     color = atom_colors[bond.a2_number],
                                     resolution = resolution,
                                     )
        fig.add_trace(bond_trace)
        
    return fig



def format_lighting (fig, ambient=0, diffuse = 1, specular = 0, roughness = 1, fresnel = 0,
                     lightx = 1000, lighty = 1000, lightz = 1000):
    
    fig.update_traces(
        lighting=dict(
            ambient=ambient,
            diffuse=diffuse,
            specular=specular,
            roughness=roughness,
            fresnel=fresnel
        ),
        lightposition=dict(
            x=lightx,
            y=lighty,
            z=lightz
        ))
    
    return fig

def format_figure (fig):
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
    margin=dict(l=0, r=0, t=0, b=0)  # Reduce margins to focus on the data
    )
    
    return fig

def draw_3D_mol (smiles = None, resolution = DEFAULT_RESOLUTION, radius = DEFAULT_RADIUS, mode = "ball+stick",
                 ambient=0, diffuse = 1, specular = 0, roughness = 1, fresnel = 0,
                 lightx = 1000, lighty = 1000, lightz = 1000):
    fig = make_subplots()
    
    fig = format_figure(fig)
    
    if smiles != None:
        atomList, bondList = process_smiles(smiles)
        
    if "ball" in mode:
        fig = draw_atoms(fig, atomList, resolution = resolution, radius = 'ball')
        if "stick" in mode:
            fig = draw_bonds(fig, bondList, resolution = resolution, radius = 'ball')
    elif "stick" == mode:
        fig = draw_atoms(fig, atomList, resolution = resolution, radius = radius)
    elif "vdw" == mode:
        fig = draw_atoms(fig, atomList, resolution = resolution*4, radius = "vdw")

    fig = format_lighting (fig, 
                           ambient=ambient, diffuse = diffuse, specular = specular, roughness = roughness, fresnel = fresnel,
                         lightx = lightx, lighty = lighty, lightz = lightz)
    
    fig.show("browser")
    
    return fig


          
       