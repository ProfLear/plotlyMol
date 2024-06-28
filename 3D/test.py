#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 11:41:46 2024

@author: benjaminlear
"""
from plotlyMol3D import draw_3D_mol

moldraw = draw_3D_mol(smiles = "CCNCO", 
            mode = "ball+stick", 
            ambient = 0.1
            ) 
