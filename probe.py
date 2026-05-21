import numpy as np
import math
import pandas as pd
from pathlib import Path
import scipy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import plotly.express as px
import plotly.graph_objects as go

#np.seterr(all='raise')
Properties = pd.read_pickle("Properties.pkl")
bodiesfolder = Path("Bodies")
files = list(bodiesfolder.glob("*.npy"))
G = 10**-3


class Body:
    def __init__(self,timeandpos,propertiesDataframe = None):
        self.array = timeandpos
        self.timestep = self.array[1,0] - self.array[0,0] #Find the timestep of this array
        self.timestart = self.array[0,0] #Start time
        self.timeend = self.array[-1,0] #End time
        self.surface_radius = 1
        self.isGravityLinear = False
        self.mass = None
        self.has_atmosphere = False
        self.air_radius = 1
        self.air_density = 1
        self.has_water = False
        self.water_radius = 1
        self.visit_radius = 2
        self.name = "DefaultName"
        if propertiesDataframe is not None:
            self.surface_radius = propertiesDataframe["surface_radius"]
            self.isGravityLinear = propertiesDataframe["isGravityLinear"]
            self.mass = propertiesDataframe["mass"]
            self.has_atmosphere = propertiesDataframe["has_atmosphere"]
            self.air_radius = propertiesDataframe["air_radius"]
            self.air_density = propertiesDataframe["air_density"]
            self.has_water = propertiesDataframe["has_water"]
            self.water_radius = propertiesDataframe["water_radius"]
            self.visit_radius = propertiesDataframe["visit_radius"]
            self.name = propertiesDataframe["name"]
        print(self)
    def __str__(self):
        string = ""
        for k,v in self.__dict__.items():
            string += f"{str(k)}: {str(v)}\n"
        return string
    def getXYZ(self,time:float):
        estimatedindex = time/self.timestep 
        np.clip(estimatedindex,0,(len(self.array)-1)) #Keep index in range
        index = math.trunc(estimatedindex) #Just going to round down
        return np.array(self.array[index,1:4])
    def getVel(self,time:float):
        start = self.getXYZ(time)
        end = self.getXYZ(time+self.timestep)
        vel = end - start
        return vel
    def converttoRealGravity(self):
        self.isGravityLinear = False
        self.mass = self.mass*self.surface_radius

class probe:
    def __init__(self,launchbody:Body,launchvel:float,launchunitvector:np.ndarray=[1,0,0],launchtime=0,endtime:float=(22+2/3),timestep:float=1/24):
        self.launchbody = launchbody
        self.launch_velcoity_mag = launchvel #Magnitude of launch velocity
        self.direction = launchunitvector #Unit vector of direction of launch velocity
        self.launchvector = self.direction*self.launch_velcoity_mag #Launch velocity vector
        self.initialvel = self.findLaunchVelVec(launchtime,self.launchvector) #Launch velocity vector accounting for the initial motion of the launch body
        self.currentlyVisiting = None #Acts as a latch to track when the probe is 'visiting' a certain body, when it gets close to a body but doesn't hit it

        self.launchtime = launchtime #When the probe is launched
        self.endtime = endtime #When to end simulation, in minutes
        self.timestep = timestep 
    def findLaunchVelVec(self,t,launchvel:np.ndarray=[0,0,0]): #Find the global cartesian vector components for launching from a body at a specific time
        return self.launchbody.getVel(t) + launchvel
    def netAcceleration(self,t,S):
        shippos = np.array([S[0],S[2],S[4]]) #This may be very perfomant but I'm not thinking about that rn
        shipvel = np.array([S[1],S[3],S[5]])
        a = np.zeros(3) #initialize acceleration
        for body in Bodies:
            diff = body.getXYZ(t) - shippos
            dist = np.linalg.norm(diff)
            #Acceleration due to gravity
            if np.isnan(body.mass): #These bodies shouldn't have gravity, the None gets turned into Nan, 
                continue
            else:
                a += G * body.mass* diff / dist**3
            #Acceleration due to drag
            if np.isnan(body.air_radius): 
                continue
            elif (body.air_radius > dist): #Inside the atmosphere
                fluidvelocity = body.getVel(t)
                relativefluidvel = fluidvelocity - shipvel
                if np.isnan(body.has_water): #Don't do anything about water if it doesn't have any
                    #Calculate air drag
                    a += calculateDrag(relativefluidvel,body.air_density,self.timestep)
                elif (body.water_radius > dist):
                    #Do water drag instead, the 30 is because water is defined to be 30
                    a += calculateDrag(relativefluidvel,30,self.timestep)
                else:
                    #Do air drag
                    a += calculateDrag(relativefluidvel,body.air_density,self.timestep)
        return a
                    
    def dSdt(self,t,S):
        x, vx, y, vy, z, vz = S
        acc = self.netAcceleration(t,S)
        return [vx,acc[0],vy,acc[1],vz,acc[2]]
    def runSimulation(self):
        initXYZ = self.launchbody.getXYZ(self.launchtime)
        
        self.path = solve_ivp(self.dSdt, [self.launchtime,(self.endtime)*60],y0 = [initXYZ[0],self.initialvel[0],initXYZ[1],self.initialvel[1],initXYZ[2],self.initialvel[2]],t_eval=np.arange(0,(self.endtime)*60,self.timestep),events=hitBody)
        print("Simulation done!")
        
def calculateDrag(relativeFluidVelocity,fluidDensity:float,dt):
    advectionmagnitude = np.linalg.norm(relativeFluidVelocity)
    dragmagintude = 0.5*fluidDensity*(advectionmagnitude)**2*0.00392*dt
    a = min(advectionmagnitude,dragmagintude)
    return a*relativeFluidVelocity/advectionmagnitude


def findProbeInitialVelocity(t=0,vel:list=[0,0,0]):# Eventually make this work with cannon launch direction
    cannonvel = Bodies[9].getVel(t)
    return cannonvel + np.asarray(vel)


def netAcceleration(t,x,y,z,vx,vy,vz):
    shippos = np.array([x,y,z]) #This may be very perfomant but I'm not thinking about that rn
    shipvel = np.array([vx,vy,vz])
    a = np.zeros(3) #initialize accerlation
    for body in Bodies:
        diff = body.getXYZ(t) - shippos
        dist = np.linalg.norm(diff)
        if np.isnan(body.surface_radius): #Gravity
            continue
        else:
            a += G * body.mass* diff / dist**3
        if np.isnan(body.air_radius): #Drag
            continue
        elif (body.air_radius > dist): #Inside the atmosphere
            fluidvelocity = body.getVel(t)
            relativefluidvel = fluidvelocity - shipvel
            if np.isnan(body.has_water): #Don't do anything about water if it doesn't have any
                #Calculate air drag
                a += calculateDrag(relativefluidvel,body.air_density,Test.timeste)
            elif (body.water_radius > dist):
                #Do water drag instead
                print("hi")
            else:
                #Do air drag
                print("air drag")
            
    return a
def sunRadiusCalculator(t):
    if t < 10*60:
        return 2000
    elif ((10*60 <= t) and (t<19*60)): #Sun goes from radius of 2000 to 4000 over 9 minutes, assuming linear growth
        return 7.407407*t-2444.4442
    else:
        return 4000
def hitBody(t,S):
    shippos = np.asarray([S[0],S[2],S[4]])
    Didnothit = 1 #1 Means it didn't hit anything, 0 means it hit something
    for body in Bodies:
        if np.isnan(body.surface_radius): #Not getting hit
            continue
        else:
            if body.name == "Sun":
                bodyradius = sunRadiusCalculator(t)
            else:
                bodyradius = body.surface_radius
            distance = np.linalg.norm(body.getXYZ(t) - shippos)
            if distance < bodyradius: #If it gets too close to the body its assumed its hit it
                Didnothit = 0
            else:
                continue #Didn't hit any body
    return Didnothit 

def random_3d_unit_vector():
    z= np.random.uniform(-1,1) #Change to -1,1
    phi = np.random.uniform(0,2*np.pi)
    r = np.sqrt(1-z**2)
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    return np.array([x,y,z])
    
def make_visit_event(body):
    def visit_event(t, S):
        shippos = np.asarray([S[0],S[2],S[4]])
        #if np.isnan(body.visi)


hitBody.terminal = True
Bodies = [] #Create list to store bodies into
Names = []
for i in range(0,len(files)): #Load in bodies
    Bodies.append(Body(np.load(files[i]),Properties.iloc[i]))
    Names.append(Bodies[i].name)
    if Bodies[i].name == "Cannon":
        Cannon = Bodies[i]
    if Properties.iloc[i]["isGravityLinear"]:
        Bodies[i].converttoRealGravity()
unitvec = random_3d_unit_vector()
print(unitvec)
mag = 500 
print(mag)
#unitvec = np.asarray([0.92224781,-0.06100891,0.38175502])
#mag = 217.4160386926305
#unitvec = np.asarray([0.93382793,0.11072035,0.34015645])
#mag = 181.33680249823882
#unitvec = np.asarray([0.22402356,0.97440414,0.01870879]) #Hit ember twin
#mag = 55.57546653762513
Test = probe(Cannon,mag,np.asarray(unitvec),0,timestep=1/60,endtime=22)
Test.runSimulation()


### Plotting
sun_radius = 2000
spherephi, spheretheta = np.mgrid[0.0:np.pi:20j, 0.0:2.0 * np.pi:20j] #Change the 20j to somethingelsej if you want different resolution on the sphere
    
    # Get Cartesian mesh grid
spherex = sun_radius*np.sin(spherephi) * np.cos(spheretheta)
spherey = sun_radius*np.sin(spherephi) * np.sin(spheretheta)
spherez = sun_radius*np.cos(spherephi)
probepath = np.zeros((len(Test.path.t),4))
probepath[:,0] = Test.path.t
probepath[:,1:4] = Test.path.y[[0,2,4],:].T

range = [-800000,800000]
step = 60*3
fig = px.scatter_3d(x=probepath[:,1][::step],y=probepath[:,2][::step],z=probepath[:,3][::step],animation_frame=probepath[:,0][::step],range_x=range,range_y=range,range_z=range) #
fig.add_trace(go.Scatter3d(
            x=probepath[:,1][::step],
            y=probepath[:,2][::step],
            z=probepath[:,3][::step],
            mode='lines',
            name="probe trajectory",
        ))
fig.add_surface(x=spherex, y=spherey, z=spherez, opacity=1.0,showscale=False)
fig.update_scenes(aspectmode='cube') #Making the axes be a cube
fig.show()
np.save("probepath.npy",probepath)