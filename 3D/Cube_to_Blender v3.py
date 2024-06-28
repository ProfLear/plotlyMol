#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 22 13:34:16 2022

@author: benjaminlear

Following this: https://github.com/alvin-yang68/Marching-Cubes/blob/main/Marching_Cubes.ipynb
"""
import collections
import numpy as np
import csv
import copy
from pathlib import Path

# #Build Marching Cubes Lookup table
#This table is straight from Paul Bourke's site (Source: http://paulbourke.net/geometry/polygonise/)
# the "-1" entriews are how we tell the program to stop eventually. 
triTable =[
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1],
            [3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1],
            [3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1],
            [3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1],
            [9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1],
            [1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1],
            [9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1],
            [2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1],
            [8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1],
            [9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1],
            [4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1],
            [3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1],
            [1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1],
            [4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1],
            [4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1],
            [9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1],
            [1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1],
            [5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1],
            [2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1],
            [9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1],
            [0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1],
            [2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1],
            [10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1],
            [4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1],
            [5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1],
            [5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1],
            [9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1],
            [0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1],
            [1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1],
            [10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1],
            [8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1],
            [2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1],
            [7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1],
            [9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1],
            [2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1],
            [11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1],
            [9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1],
            [5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1],
            [11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1],
            [11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1],
            [1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1],
            [9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1],
            [5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1],
            [2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1],
            [0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1],
            [5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1],
            [6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1],
            [0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1],
            [3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1],
            [6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1],
            [5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1],
            [1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1],
            [10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1],
            [6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1],
            [1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1],
            [8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1],
            [7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1],
            [3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1],
            [5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1],
            [0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1],
            [9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1],
            [8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1],
            [5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1],
            [0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1],
            [6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1],
            [10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1],
            [10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1],
            [8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1],
            [1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1],
            [3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1],
            [0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1],
            [10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1],
            [0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1],
            [3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1],
            [6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1],
            [9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1],
            [8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1],
            [3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1],
            [6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1],
            [0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1],
            [10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1],
            [10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1],
            [1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1],
            [2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1],
            [7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1],
            [7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1],
            [2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1],
            [1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1],
            [11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1],
            [8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1],
            [0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1],
            [7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1],
            [10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1],
            [2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1],
            [6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1],
            [7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1],
            [2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1],
            [1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1],
            [10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1],
            [10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1],
            [0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1],
            [7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1],
            [6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1],
            [8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1],
            [9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1],
            [6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1],
            [1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1],
            [4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1],
            [10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1],
            [8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1],
            [0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1],
            [1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1],
            [8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1],
            [10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1],
            [4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1],
            [10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1],
            [5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1],
            [11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1],
            [9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1],
            [6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1],
            [7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1],
            [3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1],
            [7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1],
            [9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1],
            [3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1],
            [6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1],
            [9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1],
            [1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1],
            [4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1],
            [7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1],
            [6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1],
            [3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1],
            [0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1],
            [6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1],
            [1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1],
            [0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1],
            [11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1],
            [6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1],
            [5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1],
            [9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1],
            [1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1],
            [1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1],
            [10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1],
            [0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1],
            [5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1],
            [10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1],
            [11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1],
            [0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1],
            [9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1],
            [7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1],
            [2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1],
            [8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1],
            [9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1],
            [9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1],
            [1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1],
            [9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1],
            [9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1],
            [5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1],
            [0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1],
            [10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1],
            [2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1],
            [0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1],
            [0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1],
            [9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1],
            [5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1],
            [3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1],
            [5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1],
            [8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1],
            [0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1],
            [9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1],
            [0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1],
            [1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1],
            [3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1],
            [4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1],
            [9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1],
            [11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1],
            [11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1],
            [2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1],
            [9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1],
            [3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1],
            [1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1],
            [4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1],
            [4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1],
            [0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1],
            [3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1],
            [3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1],
            [0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1],
            [9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1],
            [1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]];

def interpolation_value(v1,v2,t):
 if (v1==v2 and t==v1):
  return 0
 elif (t > v1 and t > v2):
  return None
 elif (t < v1 and t < v2):
  return None
 else:
  return (v1-t)/(v1-v2)

#Calculate x,y and z coordinates using linear interpolation
cache={}
def linear_interpolation(edge,cells,top,left,depth,thres):
 tval = 0
 point = None
  #edge 0
 if (edge == 0):
    if (((left,top,depth),(left+1,top,depth)) in cache):
        point = cache[((left,top,depth),(left+1,top,depth))]
    else:
        tval = interpolation_value (cells[left,top,depth],cells[left+1,top,depth],thres)
        if (tval is None):
            return None
        point = (left+tval,top,depth)
        cache[((left,top,depth),(left+1,top,depth))] = point
    return point

#edge 1
 if (edge == 1):
    if (((left+1,top,depth),(left+1,top+1,depth)) in cache):
        point = cache[((left+1,top,depth),(left+1,top+1,depth))]
    else:
        tval = interpolation_value (cells[left+1,top,depth],cells[left+1,top+1,depth],thres)
        if (tval is None):
            return None
        point = (left+1,top+tval,depth)
        cache[((left+1,top,depth),(left+1,top+1,depth))] = point
    return point

#edge 2
 if (edge == 2):
    if (((left,top+1,depth),(left+1,top+1,depth)) in cache):
        point = cache[((left,top+1,depth),(left+1,top+1,depth))]
    else:
        tval = interpolation_value (cells[left,top+1,depth],cells[left+1,top+1,depth],thres)
        if (tval is None):
            return None
        point = (left+tval,top+1,depth)
        cache[((left,top+1,depth),(left+1,top+1,depth))] = point
    return point

#edge 3
 if (edge == 3):
    if (((left,top,depth),(left,top+1,depth)) in cache):
        point = cache[((left,top,depth),(left,top+1,depth))]
    else:
        tval = interpolation_value (cells[left,top,depth],cells[left,top+1,depth],thres)
        if (tval is None):
            return None
        point = (left,top+tval,depth)
        cache[((left,top,depth),(left,top+1,depth))] = point
    return point

#edge 4
 if (edge == 4):
    if (((left,top,depth+1),(left+1,top,depth+1)) in cache):
        point = cache[((left,top,depth+1),(left+1,top,depth+1))]
    else:
        tval = interpolation_value (cells[left,top,depth+1],cells[left+1,top,depth+1],thres)
        if (tval is None):
            return None
        point = (left+tval,top,depth+1)
        cache[((left,top,depth+1),(left+1,top,depth+1))] = point
    return point

#edge 5
 if (edge == 5):
    if (((left+1,top,depth+1),(left+1,top+1,depth+1)) in cache):
        point = cache[((left+1,top,depth+1),(left+1,top+1,depth+1))]
    else:
        tval = interpolation_value (cells[left+1,top,depth+1],cells[left+1,top+1,depth+1],thres)
        if (tval is None):
            return None
        point = (left+1,top+tval,depth+1)
        cache[((left+1,top,depth+1),(left+1,top+1,depth+1))] = point
    return point

#edge 6
 if (edge == 6):
    if (((left,top+1,depth+1),(left+1,top+1,depth+1)) in cache):
        point = cache[((left,top+1,depth+1),(left+1,top+1,depth+1))]
    else:
        tval = interpolation_value (cells[left,top+1,depth+1],cells[left+1,top+1,depth+1],thres)
        if (tval is None):
            return None
        point = (left+tval,top+1,depth+1)
        cache[((left,top+1,depth+1),(left+1,top+1,depth+1))] = point
    return point

#edge 7
 if (edge == 7):
    if (((left,top,depth+1),(left,top+1,depth+1)) in cache):
        point = cache[((left,top,depth+1),(left,top+1,depth+1))]
    else:
        tval = interpolation_value (cells[left,top,depth+1],cells[left,top+1,depth+1],thres)
        if (tval is None):
            return None
        point = (left,top+tval,depth+1)
        cache[((left,top,depth+1),(left,top+1,depth+1))] = point
    return point

#edge 8
 if (edge == 8):
    if (((left,top,depth),(left,top,depth+1)) in cache):
        point = cache[((left,top,depth),(left,top,depth+1))]
    else:
        tval = interpolation_value (cells[left,top,depth],cells[left,top,depth+1],thres)
        if (tval is None):
            return None
        point = (left,top,depth+tval)
        cache[((left,top,depth),(left,top,depth+1))] = point
    return point

#edge 9
 if (edge == 9):
    if (((left+1,top,depth),(left+1,top,depth+1)) in cache):
        point = cache[((left+1,top,depth),(left+1,top,depth+1))]
    else:
        tval = interpolation_value (cells[left+1,top,depth],cells[left+1,top,depth+1],thres)
        if (tval is None):
            return None
        point = (left+1,top,depth+tval)
        cache[((left+1,top,depth),(left+1,top,depth+1))] = point
    return point

#edge 10
 if (edge == 10):
    if (((left+1,top+1,depth),(left+1,top+1,depth+1)) in cache):
        point = cache[((left+1,top+1,depth),(left+1,top+1,depth+1))]
    else:
        tval = interpolation_value (cells[left+1,top+1,depth],cells[left+1,top+1,depth+1],thres)
        if (tval is None):
            return None
        point = (left+1,top+1,depth+tval)
        cache[((left+1,top+1,depth),(left+1,top+1,depth+1))] = point
    return point

#edge 11
 if (edge == 11):
    if (((left,top+1,depth),(left,top+1,depth+1)) in cache):
        point = cache[((left,top+1,depth),(left,top+1,depth+1))]
    else:
        tval = interpolation_value (cells[left,top+1,depth],cells[left,top+1,depth+1],thres)
        if (tval is None):
            return None
        point = (left,top+1,depth+tval)
        cache[((left,top+1,depth),(left,top+1,depth+1))] = point
    return point

#Reference Marching Cubes: Cases Slide (Lecture Notes)
#Build the 8-bit conversion code (Modified it to generate the decimal value)
#This tells us which case we are dealing with, and we can lookup in the lookup table
def getContourCase(top,left,depth, thres,cells):
   x = 0  
   if (thres < cells[left,top+1,depth+1]):
        x = 128
   if (thres < cells[left+1,top+1,depth+1]):
        x = x + 64
   if (thres < cells[left+1,top,depth+1]):
        x = x + 32
   if (thres < cells[left,top,depth+1]):
        x = x + 16
   if (thres < cells[left,top+1,depth]):
        x = x + 8
   if (thres < cells[left+1,top+1,depth]):
        x = x + 4
   if (thres < cells[left+1,top,depth]):
        x = x + 2
   if (thres < cells[left,top,depth]):
        x = x + 1
   case_value = triTable[x] 
   return case_value  

#this is the main look, we pass an iso value (thres) and the volume data (cells)
# We get back vertices and faces, which we can use in blender to build our volume. 
def getContourSegments(thres,cells):
    rows =  cells.shape[0]  # get the number of rows
    cols = cells.shape[1]   # get the number of cols
    zcols   = cells.shape[2]# get the number of z-columns
    vertex_counter = 0 # get a rounter so we can iterate where we are.
    vertex_array = collections.OrderedDict()
    face_array = []
    #t1 = time.time()
    for left in range(0, rows-1): # since we are dealing with cubes, we only need to go to the 1 less than end. 
        for top in range(0, cols-1):
            for depth in range(0, zcols-1):
                case_val = getContourCase(top,left,depth,thres,cells)
                k = 0
                while (case_val [k] != -1):
                    v1 =  linear_interpolation(case_val [k],cells,top,left,depth,thres)
                    if v1 is None:
                        k = k + 3
                        continue
                    v2  = linear_interpolation(case_val [k+1],cells,top,left,depth,thres)
                    if v2 is None:
                        k = k + 3
                        continue
                    v3 =  linear_interpolation(case_val [k+2],cells,top,left,depth,thres)
                    if v3 is None:
                        k = k + 3
                        continue

                    k = k + 3
                    tmp = [3, 0, 0, 0]
                    if v1 not in vertex_array:
                        vertex_array[v1] = [vertex_counter, v1[0], v1[1], v1[2]]
                        tmp[1] = vertex_counter
                        vertex_counter += 1
                    else:
                        tmp[1] = vertex_array[v1][0]
                    if v2 not in vertex_array:
                        vertex_array[v2] = [vertex_counter, v2[0], v2[1], v2[2]]
                        tmp[2] = vertex_counter
                        vertex_counter += 1
                    else:
                        tmp[2] = vertex_array[v2][0]
                    if v3 not in vertex_array:
                        vertex_array[v3] = [vertex_counter, v3[0], v3[1], v3[2]]
                        tmp[3] = vertex_counter
                        vertex_counter += 1
                    else:
                        tmp[3] = vertex_array[v3][0]
                    face_array.append(tmp)
    #t2 = time.time()
    #print("\nTime taken by algorithm\n"+'-'*40+"\n{} s".format(t2-t1))
    vertex_array = np.array(list(vertex_array.values()))
    return vertex_array[:,1:], np.array(face_array)

#%% This is mine to handle the cube file...

# this will read a cube file
# it will return a res x res x res piece of information that has the values in their x, y, z coordinates. 
# but it also screens by the iso value given 
def cubefile_to_density(cubefile):  #$ turn this into a simple reader...
    # open up the file, make a file to get info from...
    with open(cubefile, mode ='r')as file:
        # reading the CSV file
        csvFile = csv.reader(file, delimiter = " ")

        read_file = []
        for line in csvFile:  
            read_file.append(list(filter(None, line))) #remove any blank spaces...
        
        resx = abs(int(read_file[3][0]))  # get the resolution we will need
        resy = abs(int(read_file[4][0]))  # get the resolution we will need
        resz = abs(int(read_file[5][0]))  # get the resolution we will need

        to_skip = abs(int(read_file[2][0])) # find out how many atoms there are...
        origin = [float(read_file[2][1]), float(read_file[2][2]), float(read_file[2][3])] # the origin for the coordinate sysutem in the cube file
        
        if int(read_file[2][0]) < 0: # check to see if we are dealing with angtroms or bohr units
            scale = 10**-10 #angstroms
        else:
            scale = 10**-9
            
    vector_scaling = [float(read_file[3][1]), float(read_file[4][2]), float(read_file[5][3])]
        
    # the numpy array that will hold our densities
    density_cube = np.zeros((resx, resy, resz))
    
    count = 0  
    for line in read_file[to_skip+7:]:  # skip to the coordinates. 
        #print(line)
        for entry in line:   # this is how we count what the indexes are...
        
        # IF THERE IS TROUBLE, LOOK HERE. 
        #this ordering makes it so that the coordinates for the orbitals agrees with the .xyz file, when both are generated by ORCA
            z = int(count % resz)
            y = int((count/resx) % resy)
            x = int((count/(resx*resy)) % resx)

            density_cube[x][y][z] = float(entry)
            
            count = count+1 #advance the count
    
    return density_cube, [origin, vector_scaling, scale] 

#takes the cube format, and returns collections of [x, y, z, and density] points
def cube_to_points(cube):
    resx = len(cube[0][0]) # get the resolution (just need one dimension, since this is a cube...)
    resy = len(cube[0])
    resz = len(cube)
    
    points = [] # to hold the columns for each phase...
    for phase in cube:
        new_points = [] # store the x, y, z, and density information
        
        count = 0
        while count < resx*resy*resz:
            x = int(count % resx)
            y = int((count/resx) % resy)
            z = int((count/(resx*resy)) % resz)
            
            if phase[x][y][z] != 0:
                new_points.append([x, y, z, phase[x][y][z]])
            
            count = count +1
         
        points.append(new_points)
    return points

# this separates out the cube file into lobes
def points_to_lobes(points):
    lobes = []
    for phase in points:
        phase_lobes = []
        to_check = [] # this will hold the points we have net to check
        
        while len(phase) > 0: # we are done, when we have found a lobe for every point
            
            #print(len(to_check))
            #print(len(phase))
            to_check.append(phase.pop(0)) #take the first point in phase, and move it to to_check
            #print(len(to_check))
            #print(len(phase))
            #int("stop by error")
            
            new_lobe = [] # this will hold the points we find for a lobe
            while len(to_check) > 0: # as long as there are entries, we are working on the same lobe
                new_lobe.append(to_check.pop(0)) #take the first one we need to check, and move it to our lobe
                
                # we can use the recently added point to see what values of x, y, and z are allowed. For some reason, needed to look two steps away... not totally sure why..
                possible_xs = [new_lobe[-1][0] + 2, new_lobe[-1][0] + 1, new_lobe[-1][0], new_lobe[-1][0] - 1, new_lobe[-1][0] - 2]  #these are the values of x neighbors could have.  The -1 index means we will check the most recently added point to the lobe
                possible_ys = [new_lobe[-1][1] + 2, new_lobe[-1][1] + 1, new_lobe[-1][1], new_lobe[-1][1] - 1, new_lobe[-1][1] - 2] # same for y
                possible_zs = [new_lobe[-1][2] + 2, new_lobe[-1][2] + 1, new_lobe[-1][2], new_lobe[-1][2] - 1, new_lobe[-1][2] - 2] # same fo z
                
                temp_neighbors = [] #keep track of the neighbors we found.  Need this because we cannot remove items fom list as we iterate through them
                for point in phase: #go thorugh the list of possible points
                    if point[0] in possible_xs and point[1] in possible_ys and point[2] in possible_zs: #check to see if this particular point is a neighbor
                        temp_neighbors.append(point)
                
                #remove dublicates from list of list, using: https://stackoverflow.com/questions/12198468/how-to-remove-duplicate-lists-in-a-list-of-list
                temp_set = set(map(tuple,temp_neighbors))  #need to convert the inner lists to tuples so they are hashable
                temp_neighbors = map(list,temp_set) #Now convert tuples back into lists (maybe unnecessary?)
                for temp in temp_neighbors:
                    if temp not in to_check and temp not in new_lobe:  # make sure this isn't already something we know about
                        to_check.append(temp)
                        phase.remove(temp)
            
            phase_lobes.append(new_lobe)        
 
        lobes.append(phase_lobes)
    
    return lobes  #will return a list with two lists, one for each phase (each of these is a list of lobes for the phase)

#this will generate a cube file similar to the one that we started with
def points_to_cube(points, res):
    cube = np.zeros((res, res, res))
    for point in points: 
        cube[point[0]] [point[1]] [point[2]] = point[3]
    return cube

def faces_to_edges (faces):
    edges = []
    for face in faces:
        if [face[0], face[1]] not in edges and [face[1], face[0]] not in edges:
            edges.append([face[0], face[1]])
            
        if [face[0], face[2]] not in edges and [face[2], face[0]] not in edges:
            edges.append([face[0], face[2]])
            
        if [face[2], face[1]] not in edges and [face[1], face[2]] not in edges:
            edges.append([face[2], face[1]])
    return edges

#need a function to trim the first entry of the faces. Not totally sure why it was included. 
def trim_faces (old_faces):
    new_faces = []
    for face in old_faces:
        new_faces.append(face[1:])
    
    return np.array(new_faces)
        
    
#%%
iso = 0.0 # this is the iso value we want to plot. 
cube_file = Path("/Users/benjaminlear/Documents/GitHub/plotlyMol/3D/anto_occ_1-min2.cube")

read_cube, scalings = cubefile_to_density(cube_file)


#all_points = cube_to_points(copy.deepcopy(read_cube))
#pos_lobes = points_to_lobes(copy.deepcopy(all_points))
#pos_cubes = points_to_cube(pos_lobes[0])

phase_verts = []
phase_faces = []

verts,faces=getContourSegments(0.02, read_cube)
#applying scaling and moving to the origin
scaled_verts = []
for v in verts: 
    temp_vert = []
    for c, s,o in zip(v, scalings[1], scalings[0]):
        temp_vert.append(0.5*(c*s +o)) #not sure why, but a factor of 1/2 is needed to get this to agree with the xyz module in blender. 
    scaled_verts.append(np.array(temp_vert))
t_faces = trim_faces(copy.deepcopy(faces))
phase_verts.append(np.array(scaled_verts))
phase_faces.append(t_faces)

verts,faces=getContourSegments(-0.02, read_cube)
#applying scaling and moving to the origin
scaled_verts = []
for v in verts: 
    temp_vert = []
    for c, s,o in zip(v, scalings[1], scalings[0]):
        temp_vert.append(0.5*(c*s +o))#not sure why, but a factor of 1/2 is needed to get this to agree with the xyz module in blender. 
    scaled_verts.append(np.array(temp_vert))
t_faces = trim_faces(copy.deepcopy(faces))
phase_verts.append(np.array(scaled_verts))
phase_faces.append(t_faces)

#edges = faces_to_edges(t_faces)


#%% for Blender...
import bpy
name = "Positive"
mesh = bpy.data.meshes.new(name)
obj = bpy.data.objects.new(name, mesh)
col = bpy.data.collections.get("Collection")  # select a collection that already exists...
col.objects.link(obj)
bpy.context.view_layer.objects.active = obj
mesh.from_pydata(phase_verts[0], [], phase_faces[0])

name = "Negative"
mesh = bpy.data.meshes.new(name)
obj = bpy.data.objects.new(name, mesh)
col = bpy.data.collections.get("Collection")  # select a collection that already exists...
col.objects.link(obj)
bpy.context.view_layer.objects.active = obj
mesh.from_pydata(phase_verts[1], [], phase_faces[1])

#%% for plotly... this works

import plotly.graph_objects as go

test = go.Figure()


# do I need this, or can I generate these more directly?
for i in range(2):
    posx, posy, posz = [], [], []
    for point in phase_verts[i]:
        posx.append(point[0])
        posy.append(point[1])
        posz.append(point[2])
    
    posi, posj, posk = [], [], []
    for face in phase_faces[i]:
        posi.append(face[0])
        posj.append(face[1])
        posk.append(face[2])

    mesh_trace = go.Mesh3d( # as long as you supply the points and connections correctly, you can put all things as one, and have it separate out. 
        x = posx, y = posy, z = posz, 
        i = posi, j = posj, k = posk
        )
    test.add_trace(mesh_trace)
    
test.update_layout(
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
test.show("browser")


