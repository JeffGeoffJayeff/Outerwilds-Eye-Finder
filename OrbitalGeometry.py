# Explanation of terms and such: https://sangillee.com/2025-01-05-elliptical-orbit-mechnics/#how-to-compute-a-planets-position-in-reference-coordinates
# Future Improvements:
# - Make it so the orbit class isn't the parent, I only did this because I developed the orbit class first and had strange ideas
#   I would undo this but that seems like a lot of effort for something that isn't going to get used that often


import math
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
G = 1*10**(-3) #Gravitational Constant
sun_mass = 4*10**11

class point:
    def __init__(self,init_x:float=None,init_y:float=None,init_z:float=None,init_pos:list = [None,None,None]):
        if (((init_x is not None) and (init_y is not None)) and (init_z is not None)):
            self.coords = [init_x,init_y,init_z]
        elif (len(init_pos) == 3 and (init_pos[0] is not None)):
            self.coords = init_pos
        else:
            print("self.coords should be a list with three entries, or there is some other issue")
    @property
    def x(self):
        return self.coords[0]
    @x.setter
    def x(self, new_x):
        self.coords[0] = new_x
    @property
    def y(self):
        return self.coords[1]
    @y.setter
    def y(self, new_y):
        self.coords[1] = new_y
    @property
    def z(self):
        return self.coords[2]
    @z.setter
    def z(self, new_z):
        self.coords[2] = new_z
    @property
    def xyz(self):
        return np.asarray(self.coords)
    @xyz.setter
    def xyz(self,new_xyz:list):
        if len(new_xyz) == 3:
            self.x = new_xyz[0]
            self.y = new_xyz[1]
            self.z = new_xyz[2]
        else:
            print("ERROR: Length isn't proper")
    @property
    def displacementMag(self):
        return math.sqrt(self.x**2+self.y**2+self.z**2)
    @property #Returns the unit vector of the displacement vector
    def unitvec(self): 
        mag = self.displacementMag
        if self.displacementMag == 0:
            return np.asarray([0,0,0])
        else:
            x = self.x/mag
            y = self.y/mag
            z = self.z/mag
            return np.asarray([x,y,z])
class orbit:
    def __init__(self,SemiMajorAxis,e,parentmass,foci:point = None,Omega:float = 0, i:float = 0, omega:float = 0, name:str = "Default", parent = None):
        if e < 0:
            print("ERROR: Eccentrcity should be greater than or equal to 0")
        elif e >= 1:
            print("ERROR: Eccentrcity cannot be 1 or greater")
        else:
            self._e = e #Eccentricity of the orbit
        self._a = SemiMajorAxis #Semi-major axis
        self._b = self._a*math.sqrt(1-self._e**2) #Semi-minor axis
        self._m = parentmass #Mass of the parent
        self._T = math.sqrt(4*(math.pi)**2*self._a**3/(G*self._m)) #Orbital period
        self._M_0 = 0 #Mean anomaly Offset, the planetary mean anomaly at the start of the loop
        self.precision = 0.000001 #Controls how precise eccentric anomaly calculator should be
        if foci is None:
            self.foci = point(0,0,0)
        else:
            self.Foci = foci
        self._Omega:float = Omega #Longitude of the ascending node
        self._i:float = i #Inclination
        self._omega:float = omega #Argument of periapsis
        self.DF = None
        self.name = name
        self.parent = parent #The parent of the object
        self.createRotationMatrix()
    def createOmegaMatrix(self): #Creates the rotation matrix for the Longitude of the Ascending node
        self.OmegaMatrix = np.asarray([[math.cos(self.AscendNodeLong), -math.sin(self.AscendNodeLong), 0],
                                       [math.sin(self.AscendNodeLong), math.cos(self.AscendNodeLong), 0],
                                       [0, 0, 1]])
    def createiMatrix(self): #Creates the rotation matrix for the inclination rotation
        self.iMatrix = np.asarray([[1,0,0],
                                   [0, math.cos(self.Inclination),-math.sin(self.Inclination)],
                                   [0, math.sin(self.Inclination),math.cos(self.Inclination)]])
    def createomegaMatrix(self): #Creates the matrix for rotating the ellipse within its orbital plane
        self.omegaMatrix = np.asarray([[math.cos(self._omega),-math.sin(self._omega),0],
                                       [math.sin(self._omega),math.cos(self._omega),0],
                                       [0,0,1]])
    def createRotationMatrix(self): #Creates a matrix that is the combination of all of them
        self.createOmegaMatrix()
        self.createiMatrix()
        self.createomegaMatrix()
        self.rotationMatrix = self.OmegaMatrix @ self.iMatrix @ self.omegaMatrix
    def findM(self,t:float): #Finding mean anomaly
        return 2*(t)*math.pi/self._T + self._M_0
    def findE(self,t:float,M:float=None):#Finding eccentric anomaly
        if M is None: 
            M = self.findM(t)
        E_old = M + (self._e*math.sin(M))/(1-self._e*math.cos(M)) #E_0 = M
        deltaE = abs(E_old - M)
        i = 0
        while((deltaE > self.precision) and (i < 10)):
            i += 1 
            E_new = E_old - (E_old-self._e*math.sin(E_old)-M)/(1-self._e*math.cos(E_old))
            deltaE = abs(E_old - E_new)
            E_old = E_new
        return E_old
    def findTheta(self,t,E:float=None): #Finding true anomaly
        if E is None:
            E = self.findE(t)
        return 2*math.atan(math.sqrt((1+self._e)/(1-self._e))*math.tan(E/2))
    def findR(self,t:float,theta:float=None):
        if theta is None:
            theta = self.findTheta(t)
        return self._a*(1-self._e**2)/(1+self._e*math.cos(theta))
    def findLocalXY(self,t): #Find xy at a specific time on the orbital plane, relative to the parent
         theta = self.findTheta(t)
         r = self.findR(t)
         x = r*math.cos(theta)
         y = r*math.sin(theta)
         return point(init_x=x,init_y=y,init_z=0)
    def findLocalXYZ(self,t): #Find xyz at a specific time in xyz coordinates relative to the parent
        orbitalPoint = self.findLocalXY(t).xyz
        XYZ = np.matmul(self.rotationMatrix,orbitalPoint) + self.foci.xyz
        return point(init_pos=XYZ.tolist())
    def useStartingPoint(self,init_X,init_Z): #This takes that starting X,Z points from the outerwilds system and calculates the properties, except for the 
        #X is the X-axis, Z is the y-axis on the orbital plane
        orbtialplanepoints = self.rotationMatrix.T @ np.asarray([init_X,0,init_Z])
        self.Eccentricity = 0
        self.SemiMajorAxis = math.sqrt(orbtialplanepoints[0]**2+orbtialplanepoints[2]**2)
        self._M_0 = math.atan2(orbtialplanepoints[2],orbtialplanepoints[0])
    def useStartingPointInterloper(self,init_X,init_Z): #Hardcoding this because this is a whole another can of worms to do
        #interloper starts at Ap at x:-24100 z:0, Pe is reported to be , which would be at positive x: 2500
        self.Eccentricity = (24100-2500)/(24100+2500)
        self.SemiMajorAxis = (24100+2500)/2
        self._M_0 = math.pi
    def createDataFrame(self,stepsize:float=1/24,endminute:float=1): #Step size is in seconds
        print(f"Creating {self.name}'s Dataframe...",end="")
        a = np.arange(0,endminute*60,stepsize)#np.linspace(0,TimberHearth._T,num=points)
        points = len(a)
        x = np.zeros(points)
        y = np.zeros(points)
        z = np.zeros(points)
        for i in range(points):
            time = a[i]
            currentPoint = self.findLocalXYZ(time)
            x[i] = currentPoint.x
            y[i] = currentPoint.y
            z[i] = currentPoint.z
        self.DF = pd.DataFrame({
            "time":a,
            "x":x,
            "y":y,
            "z":z,
            "body": np.full(points,self.name)
        })
        if (self.parent is not None):
            if (self.parent.DF is None):
                self.parent.createDataFrame(stepsize,endminute)
            self.DF["x"] = self.DF["x"] + self.parent.DF["x"]
            self.DF["y"] = self.DF["y"] + self.parent.DF["y"]
            self.DF["z"] = self.DF["z"] + self.parent.DF["z"]
        print(f"{self.name}'s Dataframe made!")
    def setParent(self,parent): #Change how this works if you ever want the system to move
        self.parent = parent
    
    @property
    def Dataframe(self):
        if self.DF is None:
            print("Warning! No dataframe has been made yet!")
            return self.DF
        else:
            return self.DF
    @property
    def SemiMajorAxis(self):
        return self._a
    @SemiMajorAxis.setter
    def SemiMajorAxis(self,new_a:float):
        if new_a > 0:
            self._a = new_a
            self._b = self._a*math.sqrt(1-self._e**2) #Semi-minor axis
            self._T = math.sqrt(4*(math.pi)**2*self._a**3/(G*self._m))
        else:
            print("ERROR: Semi-Major Axis must be positive")
    @property
    def SemiMinorAxis(self):
        return self._b
    @SemiMinorAxis.setter
    def SemiMinorAxis(self,new_b:float):
        if new_b > 0:
            self._b = new_b
            self._a = math.sqrt((self._b**2)/(1-self._e**2))
            self._T = math.sqrt(4*(math.pi)**2*self._a**3/(G*self._m))
        else:
            print("ERROR: Semi-Minor Axis must be positive")
    @property
    def Eccentricity(self):
        return self._e
    @Eccentricity.setter
    def Eccentricity(self,new_e:float):
        if new_e < 0:
            print("ERROR: New eccentricity must be atleast 0!")
        elif new_e >= 1:
            print("ERROR: New eccentricity cannot be 1 or greater!")
        else:
            self._e = new_e
    @property
    def ParentMass(self):
        return self._m
    @ParentMass.setter
    def ParentMass(self,new_mass:float):
        if new_mass < 0:
            print("ERROR: New mass cannot be less than 0")
        else:
            self._m = new_mass
            self._T = math.sqrt(4*(math.pi)**2*self._a**3/(G*self._m)) # Need to recalculate period
    @property
    def Period(self):
        return self._T
    @Period.setter
    def Period(self,new_T:float):
        if new_T < 0:
            print("ERROR: New period cannot be less than 0")
        else:
            if (self.parent.isGravityLinear): #This is some jank to calculate the appropiate mass to have the same period in actual gravity as the linear gravity
                Velocity = math.sqrt(G*self.parent.mass)
                new_T = 2*math.pi*self.SemiMajorAxis/Velocity
            self._T = new_T
            self._m = 4*self.SemiMajorAxis**3*math.pi**2/(G*self._T**2) #Recalculate mass to create desired period. I hope I did my algebra correctly!
    @property
    def Foci(self): # Foci where the parent body is
        return self.foci
    @Foci.setter
    def Foci(self,newpoint:point):
        coords = newpoint.xyz
        self.foci = point(init_x=coords[0],init_y=coords[1],init_z=coords[2])
    @property
    def DistantFoci(self): # The other foci
        coords = self.Foci.xyz - np.asarray([2*self.Eccentricity*self.SemiMajorAxis,0,0])
        return point(init_pos=coords.tolist()) #Don't know what this is doing
    @property
    def InitialMeanAnomaly(self):
        return self._M_0 
    @InitialMeanAnomaly.setter
    def InitialMeanAnomaly(self,new_M_0:float):
        self._M_0 = new_M_0 #This could be moduloed by period or something but that makes things more complicated
    @property
    def AscendNodeLong(self):
        return self._Omega
    @AscendNodeLong.setter
    def AscendNodeLong(self,new_Omega:float):
        self._Omega = new_Omega
        self.createRotationMatrix()
    @property
    def Inclination(self):
        return self._i
    @Inclination.setter
    def Inclination(self,new_i:float):
        self._i = new_i
        self.createRotationMatrix()
    @property
    def ArgofPer(self):
        return self._omega
    @ArgofPer.setter
    def ArgofPer(self,new_omega):
        self._omega = new_omega
        self.createRotationMatrix()
class planet(orbit):
    def __init__(self, SemiMajorAxis, e, parentmass, foci = None, Omega = 0, i = 0, omega = 0, 
                 name = "Default", parent=None,
                 mass:float = None,linearGravity = False,surfradius:float = None,atmoradius:float=None, air_density:float = 0, water_radius = None,visit_radius:float = 1): #Me with my giant initialization function
        super().__init__(SemiMajorAxis, e, parentmass, foci, Omega, i, omega, name, parent)
        self.mass = mass
        self.isGravityLinear = linearGravity
        self.surface_radius = surfradius #Radius of the surface, if the probe ever goes within this radius of the point of the planet it is just said to have had crashed
        self.air_radius = atmoradius #Radius of atmosphere
        self.visit_radius = visit_radius
        if atmoradius is not None:
            self.has_atmosphere = True #This should probably have some function to update it but these are short lived objects anyways
        else:
            self.has_atmosphere = False
        self.water_radius = water_radius #Basically only for giant's deep
        if water_radius is not None:
            self.has_water = True
        else:
            self.has_water = False
        self.air_density = air_density #This is based off of the value found in the game files
class StationaryPlanet(planet): #This is really only going to be used for the white hole and the sun but whatever
    def __init__(self, SemiMajorAxis=1, e=0, parentmass=1, foci=None, Omega=0, i=0, omega=0, name="Default", parent=None, mass = 10, linearGravity=False, surfradius = 1, atmoradius:float = 0, air_density = 0,anchor:point=point(0,0,0)):
        super().__init__(SemiMajorAxis, e, parentmass, foci, Omega, i, omega, name, parent, mass, linearGravity, surfradius, atmoradius, air_density)
        self.anchor = anchor
    def createDataFrame(self, stepsize = 1 / 24, endminute = 1):
        print(f"Creating {self.name}'s Dataframe...",end="")
        a = np.arange(0,endminute*60,stepsize)
        points = len(a)
        x = np.full(points,self.anchor.x)
        y = np.full(points,self.anchor.y)
        z = np.full(points,self.anchor.z)
        self.DF = pd.DataFrame({
            "time":a,
            "x":x,
            "y":y,
            "z":z,
            "body": np.full(points,self.name)
        })
        if (self.parent is not None):
            if (self.parent.DF is None):
                self.parent.createDataFrame(stepsize,endminute)
            self.DF["x"] = self.DF["x"] + self.parent.DF["x"]
            self.DF["y"] = self.DF["y"] + self.parent.DF["y"]
            self.DF["z"] = self.DF["z"] + self.parent.DF["z"]
        print(f"{self.name}'s Dataframe made!")
class MagicallyMovingPlanet(planet): #This is for planets moving at a constant velocity, none of the planets do this but I thought it would be neat to test
    # The foci in this class is where the planet is at time = 0
    def __init__(self, SemiMajorAxis=1, e=0, parentmass=1, foci=None, Omega=0, i=0, omega=0, name="Default", parent=None, mass = 10, linearGravity=False, surfradius = 1, atmoradius = 0, air_density = 0,init_vel_vec:list=[0,0,0]):
        super().__init__(SemiMajorAxis, e, parentmass, foci, Omega, i, omega, name, parent, mass, linearGravity, surfradius, atmoradius, air_density)
        self.vel_vec = init_vel_vec
    def createDataFrame(self,stepsize:float=1/24,endminute:float=1): #Step size is in seconds
        print(f"Creating {self.name}'s Dataframe...",end="")
        a = np.arange(0,endminute*60,stepsize)#np.linspace(0,TimberHearth._T,num=points)
        points = len(a)
        x = np.zeros(points)
        y = np.zeros(points)
        z = np.zeros(points)
        for i in range(points):
            time = a[i]
            currentPoint = self.findLocalXYZ(time)
            x[i] = self.Foci.x + self.vel_vec[0]*time
            y[i] = self.Foci.y + self.vel_vec[1]*time
            z[i] = self.Foci.z + self.vel_vec[2]*time
        self.DF = pd.DataFrame({
            "time":a,
            "x":x,
            "y":y,
            "z":z,
            "body": np.full(points,self.name)
        })
        if (self.parent is not None):
            if (self.parent.DF is None):
                self.parent.createDataFrame(stepsize,endminute)
            self.DF["x"] = self.DF["x"] + self.parent.DF["x"]
            self.DF["y"] = self.DF["y"] + self.parent.DF["y"]
            self.DF["z"] = self.DF["z"] + self.parent.DF["z"]
        print(f"{self.name}'s Dataframe made!")
class StrangerMotion(planet): #This is just for the stranger, which based on my understanding of the RingWorldController script, departs 405 seconds after the loop begins, and is constantly accelerating at units of 0.2 m/s^2
    def __init__(self, SemiMajorAxis, e, parentmass, foci=None, Omega=0, i=0, omega=0, name="Default", parent=None, mass = None, linearGravity=False, surfradius = None, atmoradius = None, air_density = 0, water_radius=None, visit_radius = 1):
        super().__init__(SemiMajorAxis, e, parentmass, foci, Omega, i, omega, name, parent, mass, linearGravity, surfradius, atmoradius, air_density, water_radius, visit_radius)
    def createDataFrame(self,departtime = 405,accelmag = 0.2, stepsize = 1 / 24, endminute = 1):
        print(f"Creating {self.name}'s Dataframe...",end="")
        a = np.arange(0,endminute*60,stepsize)
        points = len(a)
        x = np.zeros(points)
        y = np.zeros(points)
        z = np.zeros(points)
        self.DF = pd.DataFrame({
            "time":a,
            "x":x,
            "y":y,
            "z":z,
            "body": np.full(points,self.name)
        })
        if (self.parent is not None):
            if (self.parent.DF is None):
                self.parent.createDataFrame(stepsize,endminute)
            self.DF["x"] = self.DF["x"] + self.parent.DF["x"]
            self.DF["y"] = self.DF["y"] + self.parent.DF["y"]
            self.DF["z"] = self.DF["z"] + self.parent.DF["z"]
        print(f"{self.name}'s Dataframe made!")