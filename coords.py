import math
import numpy as np

angles = np.asarray(list(range(-45,50,5))) 
print(angles)
radius = 12
x = np.cos(np.deg2rad(angles))*radius*16
y = np.sin(np.deg2rad(angles))*radius*16

x = np.round(x)/16
y = np.round(y)/16 

for i in range(len(angles)):
    print(f"Angle:{angles[i]} - X:{x[i]} - Y:{y[i]}")