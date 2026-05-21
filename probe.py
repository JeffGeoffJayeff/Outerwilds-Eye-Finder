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
eye_distance = 286500 #Distance of the eye from the sun in meters https://www.reddit.com/r/outerwilds/comments/t7mxcy/how_far_away_is_the_eye_base_game_spoilers/
sunBodyIndex = 0 #Index that is the Sun in the Bodies list
NormalGravityforAll = True #This controls whether gravity is calculated using Newtonian gravity, or if it uses the so called linear gravity https://www.youtube.com/watch?v=dpKUoWgRBSU
n_sim = 10
Mass_Simulation_Mode = True #Whether or not you are simulating one or multiple launches
# If True then the mass for each planet is changed to produce the same gravity at the surface in both systems

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
    def __str__(self):
        string = ""
        for k,v in self.__dict__.items():
            string += f"{str(k)}: {str(v)}\n"
        return string
    def getXYZ(self,time:float): 
        estimatedindex = time/self.timestep 
        np.clip(estimatedindex,0,(len(self.array)-1)) #Keep index in range
        index = math.trunc(estimatedindex) #Just going to round down
        if index < (len(self.array[:,0])-1): #Do some linear interpolation as long as it isn't the last entry
            Y_2 = np.array(self.array[index+1,1:4])
            Y_1 = np.array(self.array[index,1:4])
            position = Y_1 + (Y_2-Y_1)/(self.timestep)*(time-self.array[index,0])
        else:
            position = np.array(self.array[index,1:4])
        return position
    def getVel(self,time:float):
        start = self.getXYZ(time)
        end = self.getXYZ(time+self.timestep)
        vel = end - start
        return vel
    def converttoRealGravity(self):
        oldMass = self.mass
        self.mass = self.mass*self.surface_radius
        self.isGravityLinear = False
        print(f"{self.name} mass has been changed from {oldMass} to {self.mass}")

class probe:
    def __init__(self,launchbody:Body,launchvel:float,launchunitvector:np.ndarray=[1,0,0],launchtime=0,endtime:float=(22+2/3),timestep:float=1/24):
        self.launchbody = launchbody
        self.launch_velcoity_mag = launchvel #Magnitude of launch velocity
        self.direction = launchunitvector #Unit vector of direction of launch velocity
        self.launchvector = self.direction*self.launch_velcoity_mag #Launch velocity vector
        self.initialvel = self.findLaunchVelVec(launchtime,self.launchvector) #Launch velocity vector accounting for the initial motion of the launch body
        self.currentlyVisiting = None #Acts as a latch to track when the probe is 'visiting' a certain body, when it gets close to a body but doesn't hit it
        self.path = None
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
            elif body.isGravityLinear: #If gravity is linear than we remove the bottom part
                a += G * body.mass * diff / dist**2
            else:
                a += G * body.mass * diff / dist**3
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
        
        self.path = solve_ivp(self.dSdt, [self.launchtime,(self.endtime)*60],y0 = [initXYZ[0],self.initialvel[0],initXYZ[1],self.initialvel[1],initXYZ[2],self.initialvel[2]],t_eval=np.arange(0,(self.endtime)*60,self.timestep),events=events)
        print("Simulation done!")
    def getXYZ(self,time:float): #TODO: Some more input handling should be added to this
        if self.path == None:
            print("ERROR: Can't get XYZ as simulation has not been run yet!")
            return
        estimatedindex = time/self.timestep 
        np.clip(estimatedindex,0,(len(self.path.t)-1)) #Keep index in range
        index = math.trunc(estimatedindex) #Just going to round down
        if index < (len(self.path.t)-1): #Do some linear interpolation as long as it isn't the last entry
            Y_2 = self.path.y[[0,2,4],index+1].T
            Y_1 = self.path.y[[0,2,4],index].T
            position = Y_1 + (Y_2-Y_1)/(self.timestep)*(time-self.path.t[index])
        else:
            position = self.path.y[[0,2,4],index].T
        return position
    def printSimulationEvents(self):
        if self.path == None:
            print("ERROR: Simulation has not been run yet!")
            return
        for i in range(len(events)):
            if i == (len(events)-1):
                print(f"Hitting Event: {self.path.t_events[i]}") #Better way to write this but who cares
            elif i == (len(events) - 2):
                print(f"Reached Eye Distance: {self.path.t_events[i]}")
            else:
                print(f"{Bodies[i].name} Visit Times: {self.path.t_events[i]}")
    def Results(self):
        output = []
        return
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
    
def make_visit_event(body:Body): 
    def visit_event(t, S):
        shippos = np.asarray([S[0],S[2],S[4]])
        if np.isnan(body.visit_radius):
            return 1 #Return 1 if there is no visit radius for the body
        else:
            distance = np.linalg.norm(body.getXYZ(t) - shippos) - body.visit_radius #negative if in visit radius, positive if not
            return distance
    visit_event.direction = -1
    visit_event.terminal = False
    return visit_event
def eyeDistance(t, S):
    shippos = np.asarray([S[0],S[2],S[4]])
    return np.linalg.norm(Bodies[sunBodyIndex].getXYZ(t)-shippos) - eye_distance #Negative if closer than the Eye is 
def cartToSpherical(coordinates:np.array): #Convert cartesian coordinates to spherical in radial, azimuthal, polar coordinates https://mathworld.wolfram.com/SphericalCoordinates.html
    coordinates = np.asarray(coordinates)
    r = np.linalg.norm(coordinates)
    theta = math.atan2(coordinates[1],coordinates[0]) #arctan(y/x)
    phi = math.acos(coordinates[2]/r) #acos(z/r)
    return np.asarray([r,theta,phi])
## Program actually starts here
hitBody.terminal = True
eyeDistance.terminal = False
eyeDistance.direction = 1 #Trigger when the probe leaves the sphere that represents the eye
Bodies = [] #Create list to store bodies into 
#Making this a global variable is a bit messed up but whatever
Names = []
for i in range(0,len(files)): #Load in bodies
    Bodies.append(Body(np.load(files[i]),Properties.iloc[i]))
    Names.append(Bodies[i].name)
    if Bodies[i].name == "Cannon":
        Cannon = Bodies[i]
events = [make_visit_event(body) for body in Bodies] #Make visiting events
events.append(eyeDistance) #Add event for reaching the distance that the eye is from the Sun
events.append(hitBody) #Add hitting event

if NormalGravityforAll:
    print("Changing masses for Newtonian gravitation...")
    for i in range(len(Bodies)):
        CurrentBody = Bodies[i]
        if CurrentBody.isGravityLinear == True: #Just making this explicit here
            if np.isnan(CurrentBody.mass):
                print(f"{CurrentBody.name} has NAN mass, skipping")
                continue
            if np.isnan(CurrentBody.surface_radius):
                print(f"{CurrentBody.name} has no surface, skipping")
                continue #Avoiding issues with NANs
            elif CurrentBody.surface_radius == 0:
                print(f"{CurrentBody.name} has a surface radius of 0, skipping")
                continue
            else:
                CurrentBody.converttoRealGravity()
        else:
            print(f"{CurrentBody.name} already has Newtonian gravity")
            continue
else:
    print("Using In-Game gravity")
## Probe settings
unitvec = random_3d_unit_vector()
print(unitvec)
mag = 500 
print(mag)

## Probe Simulation
Test = probe(Cannon,mag,np.asarray(unitvec),0,timestep=1/60,endtime=22)
Test.runSimulation()
Test.printSimulationEvents()
## Plotting
sun_radius = 2000
spherephi, spheretheta = np.mgrid[0.0:np.pi:20j, 0.0:2.0 * np.pi:20j] #Change the 20j to somethingelsej if you want different resolution on the sphere
    
# Get Cartesian mesh grid
spherex = sun_radius*np.sin(spherephi) * np.cos(spheretheta)
spherey = sun_radius*np.sin(spherephi) * np.sin(spheretheta)
spherez = sun_radius*np.cos(spherephi)
probepath = np.zeros((len(Test.path.t),4))
probepath[:,0] = Test.path.t
probepath[:,1:4] = Test.path.y[[0,2,4],:].T
print(f"Time: {Test.path.t_events[15]}, Cartesian Coordinates: {Test.getXYZ(Test.path.t_events[15][0])}, Spherical {cartToSpherical(Test.getXYZ(Test.path.t_events[15][0]))}")
range = [-800000,800000]
step = 60*1 #Step in stepsizes
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

# Create res
results = np.empty(
    n_sim,
    dtype = [
        ("Launch x",np.float32),
        ("Launch y",np.float32),
        ("Launch z",np.float32),
        ("Launch Velocity",np.float32),
        ("Reached Eye",np.bool_),
        ("Eye Shell Time",np.float32), #The time the probe reaches 286 km or whatever it is
        ("Eye Shell Polar",np.float32), #The polar of the above point
        ("Eye Shell Azimuth",np.float32), #Azimuth of the above point
        ("Sun Visits",np.int16),
        ("Sun Station Visits",np.int16),
        ("Ember Twin Visits",np.int16),
        ("Ash Twin Visits", np.int16),
        ("Timber Hearth Visits",np.int16),
        ("Attlerock Visits", np.int16),
        ("Brittle Hollow Visits", np.int16),
        ("Hollow's Lantern Visits", np.int16),
        ("Giant's Deep Visits",np.int16),
        ("Cannon Visits",np.int16),
        ("Dark Bramble Visits",np.int16),
        ("White Hole Visits", np.int16),
        ("Random Eye Visits", np.int16),
        ("Hit Something", np.bool_), #Whether or not the probe hit something
        ("Body Hit", "U2")
    ]
)
print(results)