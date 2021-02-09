# -*- coding:utf-8 -*-
"""
@author:yoghurt
@file:PS6_1.py
@@time: 2020/11/25 10:34
"""

import numpy as np
import matplotlib.pyplot as plt 


nx = 20
ny = 20
Lx = 2
Ly = 2
deltax = int(2/nx)
deltay = int(2/ny)
deltat = 0.05

c = np.empty((nt, nx, ny), dtype=float)
c[0, :, :] = 0.
for i in range(0, nx):
    

