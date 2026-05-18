import numpy as np

vcannon = np.asarray([271.92222,0,122.0461])
vprobe1 = np.asarray([224.1179,92.7537,647.428])
vlaunch1 = vprobe1 - vcannon
print(np.linalg.norm(vlaunch1))
vprobe2 = np.asarray([-54.3396,-222.2607,517.5892])
vprobe3 = np.asarray([-473.6798,-71.2808,108.0052])
vprobe4 = np.asarray([327.1968,-175.6138,449.6992])
vprobe5= np.asarray([-101.9103, -169.5816, -343.5451])
vlaunch2 = vprobe2 - vcannon
vlaunch3 = vprobe3 - vcannon
vlaunch4 = vprobe4 - vcannon
vlaunch5 = vprobe5 - vcannon

print(np.linalg.norm(vlaunch2))
print(np.linalg.norm(vlaunch3))
print(np.linalg.norm(vlaunch4))
print(np.linalg.norm(vlaunch5))
print("Probe final speeds")
print(np.linalg.norm(vprobe1))
print(np.linalg.norm(vprobe2))
print(np.linalg.norm(vprobe3))
print(np.linalg.norm(vprobe4))
print(np.linalg.norm(vprobe5))