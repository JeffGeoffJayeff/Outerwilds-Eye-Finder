

import numpy as np
import math
import pandas as pd
from pathlib import Path
from datetime import datetime
import time
import scipy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import plotly.express as px
import plotly.graph_objects as go
from multiprocessing import Pool
import multiprocessing as mp

#np.seterr(all='raise')
Properties = pd.read_pickle("Properties.pkl")
bodiesfolder = Path("Bodies")
files = list(bodiesfolder.glob("*.npy"))
G = 10**-3
eye_distance = 286500 #Distance of the eye from the sun in meters https://www.reddit.com/r/outerwilds/comments/t7mxcy/how_far_away_is_the_eye_base_game_spoilers/
sunBodyIndex = 0 #Index that is the Sun in the Bodies list
NormalGravityforAll = True #This controls whether gravity is calculated using Newtonian gravity, or if it uses the so called linear gravity https://www.youtube.com/watch?v=dpKUoWgRBSU
n_sim_per_pikmin = 250 #number of simulations to run per pikmin, where a pikmin is a multiprocessing worker, multiple launches is done per worker to reduce the overhead of starting a new process for each launch
total_n_pikmin_to_make = 800 #Total number of pikmin to make, this is the total number of processes that will be made, each pikmin will run n_sim_per_pikmin simulations
pikmin_on_field = None #Number of pikmin to run at once, this is the number of processes that will be running at once, if this is set to 1 then it will run in serial, if it is set to 4 then it will run 4 simulations at once, and so on, based on cores or something
Mass_Simulation_Mode = True #Whether or not you are simulating one or multiple launches
# If True then the mass for each planet is changed to produce the same gravity at the surface in both systems
plotPath = False #Whether to plot or not

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
    def __init__(self,launchbodyindex:int,launchvel:float,Bodies:list[Body],launchunitvector:np.ndarray=np.asarray([1,0,0]),launchtime=0,endtime:float=(22+2/3),timestep:float=1/24):
        self.launchbody = Bodies[launchbodyindex]
        self.launch_velocity_mag = launchvel #Magnitude of launch velocity
        self.direction = launchunitvector #Unit vector of direction of launch velocity
        self.launchvector = self.direction*self.launch_velocity_mag #Launch velocity vector
        self.initialvel = self.findLaunchVelVec(launchtime,self.launchvector) #Launch velocity vector accounting for the initial motion of the launch body
        self.path = None
        self.launchtime = launchtime #When the probe is launched
        self.endtime = endtime #When to end simulation, in minutes
        self.timestep = timestep 
        self.Bodies = Bodies
        self.events = [] #List of event functions to be used in the solve_ivp function
    ## Setting Up Simulation
    def findLaunchVelVec(self,t,launchvel:np.ndarray=[0,0,0]): #Find the global cartesian vector components for launching from a body at a specific time
        return self.launchbody.getVel(t) + launchvel
    def ChangeLaunchConditions(self,launchbodyindex:int,launchvel:float,launchunitvector:np.ndarray=np.asarray([1,0,0]),launchtime=0):
        self.launchbody = self.Bodies[launchbodyindex]
        self.launchtime = launchtime #When the probe is launched
        self.launch_velocity_mag = launchvel #Magnitude of launch velocity
        self.direction = launchunitvector #Unit vector of direction of launch velocity
        self.launchvector = self.direction*self.launch_velocity_mag #Launch velocity vector
        self.initialvel = self.findLaunchVelVec(launchtime,self.launchvector) #Launch velocity vector accounting for the initial motion of the launch body
    ## Simulation Functions
    def netAcceleration(self,t,S):
        shippos = np.array([S[0],S[2],S[4]]) #This may be very perfomant but I'm not thinking about that rn
        shipvel = np.array([S[1],S[3],S[5]])
        a = np.zeros(3) #initialize acceleration
        for body in self.Bodies:
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
                    a += calculateDrag(relativefluidvel,body.air_density)
                elif (body.water_radius > dist):
                    #Do water drag instead, the 30 is because water is defined to be 30
                    a += calculateDrag(relativefluidvel,30)
                else:
                    #Do air drag
                    a += calculateDrag(relativefluidvel,body.air_density)
        return a   
    def make_hitBody_event(self):            
        def hitBody(t,S):
            shippos = np.asarray([S[0],S[2],S[4]])
            Didnothit = 1 #1 Means it didn't hit anything, 0 means it hit something
            for body in self.Bodies:
                if np.isnan(body.surface_radius): #Not getting hit
                    continue
                else:
                    if body.name == "Sun":
                        bodyradius = calculateSunRadius(t)
                    else:
                        bodyradius = body.surface_radius
                    distance = np.linalg.norm(body.getXYZ(t) - shippos)
                    if distance < bodyradius: #If it gets too close to the body its assumed its hit it
                        Didnothit = 0
                    else:
                        continue #Didn't hit any body
            return Didnothit 
        hitBody.terminal = True
        return hitBody
    def make_visit_event(self,body:Body): 
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
    def make_eyeDistance_event(self):
        def eyeDistance(t, S):
            shippos = np.asarray([S[0],S[2],S[4]])
            return np.linalg.norm(self.Bodies[sunBodyIndex].getXYZ(t)-shippos) - eye_distance #Negative if closer than the Eye is 
        eyeDistance.terminal = False
        eyeDistance.direction = 1 #Trigger when the probe leaves the sphere that represents the eye
        return eyeDistance
    def dSdt(self,t,S):
        x, vx, y, vy, z, vz = S
        acc = self.netAcceleration(t,S)
        return [vx,acc[0],vy,acc[1],vz,acc[2]]
    def runSimulation(self,printoutput:bool=False):
        initXYZ = self.launchbody.getXYZ(self.launchtime)
        if len(self.events) == 0:
        # Setting up events
            self.events = [self.make_visit_event(body) for body in self.Bodies] #Make visiting events
            self.events.append(self.make_eyeDistance_event()) #Add event for reaching the distance that the eye is from the Sun
            self.events.append(self.make_hitBody_event()) #Add hitting event
        if printoutput:
            print("Running simulation...")
        self.path = solve_ivp(self.dSdt, [self.launchtime,(self.endtime)*60],y0 = [initXYZ[0],self.initialvel[0],initXYZ[1],self.initialvel[1],initXYZ[2],self.initialvel[2]],t_eval=np.arange(0,(self.endtime)*60,self.timestep),events=self.events)
        if printoutput:
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
        for i in range(len(self.events)):
            if i == (len(self.Bodies)+1): #Body visit events are always first
                print(f"Hitting Event: {self.path.t_events[i]}") #Better way to write this but who cares
            elif i == (len(self.Bodies)):
                print(f"Reached Eye Distance: {self.path.t_events[i]}")
            else:
                print(f"{self.Bodies[i].name} Visit Times: {self.path.t_events[i]}")
    
    def Results(self):
        output = [] #23
        output.extend(self.direction) #Adding Launch XYZ
        output.append(self.launch_velocity_mag) #Adding launch velocity
        output.append(self.arrivedAtEye) #Eye Tracking stuff
        if self.arrivedAtEye: 
            output.append(self.eyeArrivalTime)
            cartcoords =  self.getXYZ(self.eyeArrivalTime)
            sphericalcoords = cartToSpherical(cartcoords)
            output.extend([sphericalcoords[1],sphericalcoords[2]])
        else:
            output.extend([np.nan,np.nan,np.nan])
        #Visting Bodies stuff
        for i in range(len(self.Bodies)):
            output.append(self.path.t_events[i].size)
        output.append(self.hitSomething)#Hit something
        if self.hitSomething:
            output.append(self.closestObjectIndex(self.hitTime))
        else:
            output.append(np.nan)
        return tuple(output)
    def closestObjectIndex(self,time:float):
        distances = np.zeros(len(self.Bodies))
        for i in range(len(self.Bodies)):
            if i == 0: #there is a way tomake this better but I won't
                SurfaceRadius = calculateSunRadius(time)
            else:
                SurfaceRadius = self.Bodies[i].surface_radius
            if np.isnan(SurfaceRadius):
                distances[i] = np.linalg.norm(self.Bodies[i].getXYZ(time) - self.getXYZ(time))
            else:
                distances[i] = max(np.linalg.norm(self.Bodies[i].getXYZ(time) - self.getXYZ(time)) - SurfaceRadius,0) #Keep values above 0, turn negative values in 0
        return np.argmin(distances)
    @property
    def eyeArrivalTime(self):
        if self.path == None:
            print("ERROR: Trying to get eye arrival time before the probe has been simulated!")
            return
        elif self.arrivedAtEye:
            return self.path.t_events[len(self.Bodies)][0]
        else: #Didn't reach eye
            return np.nan
    @property
    def arrivedAtEye(self): #Probably a better way to write this
        if self.path == None:
            print("ERROR: 001 Simulate the probe first")
            return np.nan
        else:
            return self.path.t_events[len(self.Bodies)].size > 0
    @property
    def hitSomething(self):
        if self.path == None:
            print("ERROR: 002 Simulate the probe first")
            return np.nan
        else:
            return self.path.t_events[len(self.Bodies)+1].size > 0
    @property
    def hitTime(self):
        if self.path == None:
            print("ERROR: 003 Simulate the probe first")
            return np.nan
        elif self.hitSomething:
            return self.path.t_events[len(self.Bodies)+1][0]
        else: #Didn't hit anything
            return np.nan 
def calculateDrag(relativeFluidVelocity,fluidDensity:float):
    advectionmagnitude = np.linalg.norm(relativeFluidVelocity)
    if advectionmagnitude < 1e-12: #If the relative fluid velocity is too small, then we can just return 0 drag
        return np.array([0,0,0])
    dragmagnitude = 0.5*fluidDensity*(advectionmagnitude)*0.00392 #originally advection magnitude was squared but since when we return the vector we divide by the magnitude to get the unit vector, we can just use the magnitude instead of the squared magnitude
    return dragmagnitude*relativeFluidVelocity 
def calculateDragTest():
    #Zero relative velocity test
    print(f"Zero relative velocity test: {calculateDrag(np.array([0,0,0]),1)}")
    #Direction test
    v = np.array([10,0,0])
    a = calculateDrag(v,1)
    print(f"Direction test: {a}, should be in the same direction as {v} : {np.dot(a,v) < 0}")
    #Scaling test
    v1 = np.array([1,0,0])
    v2 = np.array([2,0,0])
    v3 = np.array([3,0,0])
    print(f"Scaling test: {calculateDrag(v1,1)}, {calculateDrag(v2,1)}, {calculateDrag(v3,1)}")
    #Extreme test
    v = np.array([1e6,0,0])
    print(f"Extreme test: {calculateDrag(v,1)}")
    #Symmetry test
    v = np.array([1,1,1])
    print(f"Symmetry test: {-calculateDrag(v,1)}, {calculateDrag(-v,1)}")
    #Rotation invariance test
    v1 = np.array([1,0,0])
    v2 = np.array([0,1,0])
    v3 = np.array([0,0,1])
    v4 = np.array([0.707,0.707,0])
    print(f"Rotation invariance test: {calculateDrag(v1,1)}, {calculateDrag(v2,1)}, {calculateDrag(v3,1)}, {calculateDrag(v4,1)} | Should have equal magnitudes")
def calculateSunRadius(t):
    if t < 10*60:
        return 2000
    elif ((10*60 <= t) and (t<19*60)): #Sun goes from radius of 2000 to 4000 over 9 minutes, assuming linear growth
        return 7.407407*t-2444.4442
    else:
        return 4000
def random_3d_unit_vector():
    phi = np.random.uniform(0,np.pi*2)
    costheta = np.random.uniform(-1,1)

    theta = np.arccos( costheta )
    x = np.sin( theta) * np.cos( phi )
    y = np.sin( theta) * np.sin( phi )
    z = np.cos( theta )
    return np.array([x,y,z]) 
def cartToSpherical(coordinates:np.array): #Convert cartesian coordinates to spherical in radial, azimuthal, polar coordinates https://mathworld.wolfram.com/SphericalCoordinates.html
    coordinates = np.asarray(coordinates)
    r = np.linalg.norm(coordinates)
    theta = math.atan2(coordinates[1],coordinates[0]) #arctan(y/x)
    phi = math.acos(coordinates[2]/r) #acos(z/r)
    return np.asarray([r,theta,phi])
## Program actually starts here
def makeResultsTemplate(length:int):
    # Create results template for storing simulation results
    resultsTemplate = np.empty(
        length,
        dtype = [
            ("Launch x",np.float64),
            ("Launch y",np.float64),
            ("Launch z",np.float64),
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
            ("Interloper Visits", np.int16),
            ("White Hole Visits", np.int16),
            ("Stranger Visits", np.int16),
            ("Random Eye Visits", np.int16),
            ("Spacey Visits", np.int16),
            ("Hit Something", np.bool_), #Whether or not the probe hit something
            ("Body Hit", np.float16) #index of whatever body it hit
        ]
    )
    return resultsTemplate
def simulationPikmin(cannonIndex:int,launchMag:float,bodiesList:list[Body],launchUnitVector:np.ndarray,launchTime:float,timestep:float,endtime:float,n_sims:int,outputdir:str,printoutput:bool=False):
    starttime = time.time()
    if printoutput: 
        print(f"[{mp.current_process().name}] Starting simulation with {n_sims} runs...")
    results = makeResultsTemplate(n_sims) 
    pikmin = probe(launchbodyindex=cannonIndex,launchvel=launchMag,Bodies=bodiesList,launchunitvector=launchUnitVector,launchtime=launchTime,endtime=endtime,timestep=timestep)
    for i in range(n_sims):
        #How to change the simulation each run
        newUnitVector = random_3d_unit_vector()
        pikmin.ChangeLaunchConditions(cannonIndex,launchMag,newUnitVector,launchTime)
        #Run the simulation
        pikmin.runSimulation(False)
        results[i] = pikmin.Results()
    if printoutput: 
        print(f"[{mp.current_process().name}] Simulation took {time.time() - starttime} seconds.")
    #Generatng a unique filename for the results and save
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    filename = f"{outputdir}/results_{timestamp}"
    np.save(f"{filename}.npy",results)
    np.savetxt(f"{filename}.csv",results,delimiter=",")
    return 

if Mass_Simulation_Mode:
    if __name__ == "__main__":
        starttime = time.time()
        print(f"Starting to compute {total_n_pikmin_to_make*n_sim_per_pikmin} simulations at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        Bodies = [] #Create list to store bodies into 
        #Making this a global variable is a bit messed up but whatever
        Names = []
        for i in range(0,len(files)): #Load in bodies
            xyzarray = np.load(files[i])
            xyzarray.flags.writeable = False #Make the array read-only to prevent accidental modification
            Bodies.append(Body(xyzarray,Properties.iloc[i]))
            Names.append(Bodies[i].name)
            if Bodies[i].name == "Cannon":
                CannonIndex = i

        if NormalGravityforAll: #Change masses to have the same surface gravity as in the linear gravity system
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

        outputdir = "Outputs"
        Path(outputdir).mkdir(parents=True, exist_ok=True)

        print(f"Starting {total_n_pikmin_to_make} pikmin with {n_sim_per_pikmin} runs each...")
        if pikmin_on_field is None:
            pikmin_on_field = mp.cpu_count() #If not specified, use all available cores
        unitvec = random_3d_unit_vector()
        
        with mp.Pool(processes=pikmin_on_field) as pool:
            for _ in range(total_n_pikmin_to_make):
                pool.apply_async(simulationPikmin, args=(CannonIndex, 500, Bodies, unitvec, 0, 1/60, 22, n_sim_per_pikmin, outputdir, True))
            pool.close()
            pool.join()
        print("All pikmin have finished their simulations.")
        print(f"Computation of {total_n_pikmin_to_make * n_sim_per_pikmin} simulations ended at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}, total time: {time.time() - starttime} seconds.")
else:
    Bodies = [] #Create list to store bodies into 
    #Making this a global variable is a bit messed up but whatever
    Names = []
    for i in range(0,len(files)): #Load in bodies
        xyzarray = np.load(files[i])
        xyzarray.flags.writeable = False #Make the array read-only to prevent accidental modification
        Bodies.append(Body(xyzarray,Properties.iloc[i]))
        Names.append(Bodies[i].name)
        if Bodies[i].name == "Cannon":
            CannonIndex = i

    if NormalGravityforAll: #Change masses to have the same surface gravity as in the linear gravity system
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
    unitvec = [0.066543, 0.997782, 0.001502]#random_3d_unit_vector()
    print(unitvec)
    mag = 500
    print(mag)
    calculateDragTest()
    ## Probe Simulation
    Test = probe(CannonIndex,mag,Bodies,np.asarray(unitvec),0,timestep=1/60,endtime=22)
    Test.runSimulation()
    Test.printSimulationEvents()
    print(Test.Results())
## Plotting
    if plotPath:
    # Get Cartesian mesh grid
        sun_radius = 2000
        spherephi, spheretheta = np.mgrid[0.0:np.pi:20j, 0.0:2.0 * np.pi:20j] #Change the 20j to somethingelsej if you want different resolution on the sphere
        spherex = sun_radius*np.sin(spherephi) * np.cos(spheretheta)
        spherey = sun_radius*np.sin(spherephi) * np.sin(spheretheta)
        spherez = sun_radius*np.cos(spherephi)
        probepath = np.zeros((len(Test.path.t),4))
        probepath[:,0] = Test.path.t
        probepath[:,1:4] = Test.path.y[[0,2,4],:].T
        print(f"Time: {Test.path.t_events[15]}, Cartesian Coordinates: {Test.getXYZ(Test.eyeArrivalTime)}, Spherical {cartToSpherical(Test.getXYZ(Test.eyeArrivalTime))}")
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

    
    