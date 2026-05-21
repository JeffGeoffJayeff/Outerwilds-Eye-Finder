import OrbitalGeometry as OG
import math
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

# Purpose: The point of this file is to calculate the position of each body at a specific time, and output them to a file 
# These files are the ones in "Bodies"
# It also graphs their positions as refactoring it into a separate program is an entire thing


G = 1*10**(-3) #Gravitational Constant
sun_mass = 4*10**11
sun_radius = 2000
### Config Stuff
Stepsize = 1 #In seconds
EndMinute = 23 #Timeloop ends at 22:40 in the games, in minutes btw
graphresults = True #Make a graph of the planet positions? (Doesn't do well with lots of points)
Savemotion = False #Save the position of all bodies to file?
Path = graphresults #Shows the trajectory of probe path
BodyPaths = graphresults #Whether to show the path of the bodies in the system or not

Twin_period = math.pi*2*250/28.28427
Twin_mass = 1.6*10**6

# Yes the dataframe step could've been done with a loop but I added these planets slowly over time so the need didn't arrise until later
Sun = OG.StationaryPlanet(mass=sun_mass,anchor=OG.point(0,0,0),name="Sun")
Sun.createDataFrame(endminute=EndMinute,stepsize=Stepsize)
Sun.surface_radius = sun_radius
Sun.visit_radius = 4500

SunStation = OG.planet(50,0,sun_mass,name="Sun Station",parent=Sun)
SunStation.surface_radius = 20
SunStation.visit_radius = 50
SunStation.useStartingPoint(-0.001220703,-2296)
SunStation.createDataFrame(endminute=EndMinute,stepsize=Stepsize)


FocalBody = OG.planet(50,0,sun_mass,name="Focal Body",parent=Sun)
FocalBody.useStartingPoint(-2867.882,4095.762)
FocalBody.isGravityLinear = False
FocalBody.mass = Twin_mass
FocalBody.createDataFrame(endminute=EndMinute,stepsize=Stepsize,)
FocalBody.visit_radius = 1000

CaveTwin = OG.planet(250,0,FocalBody.mass,name="Ember Twin",parent=FocalBody,mass=1.6*10**6,linearGravity=True)
CaveTwin.Inclination = math.pi
CaveTwin.useStartingPoint(-204.788,-143.394)
CaveTwin.Period = Twin_period
CaveTwin.createDataFrame(endminute=EndMinute,stepsize=Stepsize)
CaveTwin.surface_radius = 169
CaveTwin.air_radius = 250
CaveTwin.air_density = 1.2
CaveTwin.visit_radius = 300

TowerTwin = OG.planet(250,0,FocalBody.mass,name="Ash Twin",parent=FocalBody,mass=1.6*10**6,linearGravity=True)
TowerTwin.Inclination = math.pi
TowerTwin.useStartingPoint(204.788,143.394)
TowerTwin.Period = Twin_period
TowerTwin.createDataFrame(endminute=EndMinute,stepsize=Stepsize)
TowerTwin.surface_radius = 170
TowerTwin.air_radius = 250
TowerTwin.air_density = 1.2
TowerTwin.visit_radius = 300

TimberHearth = OG.planet(50,0.7,sun_mass,foci=OG.point(init_pos=[0,0,0]),name="Timber Hearth",parent=Sun)
TimberHearth.isGravityLinear = True
TimberHearth.mass = 3*10**6
TimberHearth.useStartingPoint(1492.172,-8462.538)
TimberHearth.createDataFrame(endminute=EndMinute,stepsize=Stepsize)
TimberHearth.surface_radius = 254
TimberHearth.air_radius = 380
TimberHearth.air_density = 1.2
TimberHearth.visit_radius = 600

Attlerock = OG.planet(900,0,sun_mass,name="Attlerock",parent=TimberHearth,mass=5*10**7)
Attlerock.Period = 50
Attlerock.surface_radius = 80
Attlerock.visit_radius = 160
Attlerock.useStartingPoint(886.327,156.283)
Attlerock.createDataFrame(endminute=EndMinute,stepsize=Stepsize)

BrittleHollow = OG.planet(50,0,sun_mass,foci=OG.point(init_pos=[0,0,0]),name="Brittle Hollow",parent=Sun)
BrittleHollow.mass = 3*10**6
BrittleHollow.isGravityLinear = True
BrittleHollow.useStartingPoint(11513.28,-2030.102)
BrittleHollow.createDataFrame(endminute=EndMinute,stepsize=Stepsize)
BrittleHollow.surface_radius = 272
BrittleHollow.air_radius = 500
BrittleHollow.air_density = 1.2
BrittleHollow.visit_radius = 600

HollowLantern = OG.planet(1000,0,sun_mass,parent=BrittleHollow,name="Hollow's Lantern",mass=9.10*10**5,linearGravity=True)
HollowLantern.Period = 84
HollowLantern.useStartingPoint(984.81,-173.648) #All of the moon starting points were calculated by taking the absolute starting position of the moon and subtracting the absolute position of the planet
HollowLantern.createDataFrame(endminute=EndMinute,stepsize=Stepsize)
HollowLantern.surface_radius = 97.3
HollowLantern.air_radius = 150
HollowLantern.air_density = 0.6
HollowLantern.visit_radius = HollowLantern.surface_radius*2

GiantsDeep = OG.planet(10,0,sun_mass,name="Giant's Deep",parent=Sun,mass=2.18*10**7)
GiantsDeep.mass = 2.18*10**7
GiantsDeep.isGravityLinear = True
GiantsDeep.useStartingPoint(3421.723,-16097.95)
GiantsDeep.createDataFrame(endminute=EndMinute,stepsize=Stepsize)
GiantsDeep.surface_radius = 200#Where the core starts
GiantsDeep.air_radius = 950
GiantsDeep.water_radius = 500
GiantsDeep.visit_radius = 1000

ORP = OG.planet(10,0,GiantsDeep.mass,parent=GiantsDeep,name="Cannon")
ORP.useStartingPoint(-1006.406,653.57)
ORP.Period = 1
ORP.createDataFrame(endminute=EndMinute,stepsize=Stepsize)

DarkBramble = OG.planet(50,0,sun_mass,foci=OG.point(init_pos=[0,0,0]),name="Dark Bramble",parent=Sun,mass=3.25*10**6)
DarkBramble.useStartingPoint(-3472.959,19696.16)
DarkBramble.isGravityLinear = True
DarkBramble.createDataFrame(endminute=EndMinute,stepsize=Stepsize)
DarkBramble.surface_radius = 203.3
DarkBramble.visit_radius = 500

Interloper = OG.planet(50,0,sun_mass,foci=OG.point(init_pos=[0,0,0]),name="The Interloper",parent=Sun,mass=5.50*10**5)
Interloper.useStartingPointInterloper(-24100,0)
Interloper.isGravityLinear = True
Interloper.Inclination = math.pi
Interloper.createDataFrame(endminute=EndMinute,stepsize=Stepsize)
Interloper.surface_radius = 83
Interloper.visit_radius = Interloper.surface_radius*2


WhiteHole = OG.StationaryPlanet(anchor=OG.point(-23000,0,0), name="White Hole",parent=Sun)
WhiteHole.createDataFrame(endminute=EndMinute,stepsize=Stepsize)
WhiteHole.visit_radius = 500

RingWorld = OG.StrangerMotion(foci=OG.point(8168.197,2049.528,8400),name="The Stranger",parent=Sun)
RingWorld.createDataFrame(endminute=EndMinute,stepsize=Stepsize)
RingWorld.surface_radius = 300
RingWorld.visit_radius = 800

TheEye = OG.planet(286500,0,sun_mass,foci=OG.point(init_pos=[0,0,0]),name="The Eye",parent=Sun,mass = 9*10**6)
TheEye.createDataFrame(endminute=EndMinute,stepsize=Stepsize)
#All bodies to be displayed/saved
BodiesList = [Sun,
              SunStation, 
              CaveTwin,
              TowerTwin,
              TimberHearth,
              Attlerock,
              BrittleHollow,
              HollowLantern,
              GiantsDeep,
              ORP,
              DarkBramble,
              Interloper,
              WhiteHole,
              RingWorld,
              TheEye]
if Savemotion:
    numofbody = len(BodiesList)
    for i in range(numofbody):
        currentbody = BodiesList[i]
        np.save(f"Bodies/{i:03}",currentbody.Dataframe.filter(["time","x","y","z"],axis=1).to_numpy())
        print(f" \"{currentbody.name}\" saved as Bodies/{i:03}.npy")
    propertiestable = pd.DataFrame({"surface_radius":[i.surface_radius for i in BodiesList],
                "isGravityLinear":[i.isGravityLinear for i in BodiesList],
                "mass":[i.mass for i in BodiesList],
                "has_atmosphere":[i.has_atmosphere for i in BodiesList],
                "air_radius":[i.air_radius for i in BodiesList],
                "air_density":[i.air_density for i in BodiesList],
                "has_water":[i.has_water for i in BodiesList],
                "water_radius":[i.water_radius for i in BodiesList],
                "visit_radius":[i.visit_radius for i in BodiesList],
                "name":[i.name for i in BodiesList]})
    propertiestable.to_pickle("Properties.pkl")
if graphresults:
    # Create mesh grid for spherical coordinates
    spherephi, spheretheta = np.mgrid[0.0:np.pi:20j, 0.0:2.0 * np.pi:20j] #Change the 20j to somethingelsej if you want different resolution on the sphere
    
    # Get Cartesian mesh grid
    spherex = sun_radius*np.sin(spherephi) * np.cos(spheretheta)
    spherey = sun_radius*np.sin(spherephi) * np.sin(spheretheta)
    spherez = sun_radius*np.cos(spherephi)
    df = pd.concat([i.Dataframe for i in BodiesList])

    df["x"] = pd.to_numeric(df["x"], errors="coerce")
    df["y"] = pd.to_numeric(df["y"], errors="coerce")
    df["z"] = pd.to_numeric(df["z"], errors="coerce")
    range = [-300000,300000]
    colormap = {
        Sun.name:"#FFDF22",
        SunStation.name:"#9D00FF",
        FocalBody.name:"#D3D3D3",
        CaveTwin.name:"#f05f44",
        TowerTwin.name:"#e19d49",
        TimberHearth.name: "#5a8240",
        Attlerock.name: "#8e7263",
        BrittleHollow.name: "#647297",
        HollowLantern.name: "#d86f23",
        GiantsDeep.name: "#31a174",
        ORP.name:"#676e4c",
        DarkBramble.name: "#55503a",
        Interloper.name: "#349fb9",
        WhiteHole.name: "#D3D3D3",
        RingWorld.name: "#22a185",
        "Probe":"#FF8C00",
        TheEye.name: "#8673A1"
    }
    if Path:
        step = Stepsize*24
        path = np.load("probepath.npy")[::int(step)]
        probeDF = pd.DataFrame({
            "time":path[:,0],
            "x":path[:,1],
            "y":path[:,2],
            "z":path[:,3],
            "body": np.full(len(path[:,0]),"Probe")
        })
        df = pd.concat([df,probeDF])
    fig = px.scatter_3d(df,x="x",y="y",z="z",animation_frame="time",color="body",range_x=range,range_y=range,range_z=range,color_discrete_map=colormap)
    if BodyPaths:
        for body, group in df.groupby('body'): #Honestly this code is ChatGPT :(
            fig.add_trace(go.Scatter3d(
                x=group['x'],
                y=group['y'],
                z=group['z'],
                mode='lines',
                name=f"{body} trajectory",
                line=dict(width=3,color=colormap[body]),
                legendgroup=body,   # ties it to the same legend entry
                showlegend=False    # prevents duplicate legend items
            ))
    elif Path:
        fig.add_trace(go.Scatter3d(
            x=path[:,1],
            y=path[:,2],
            z=path[:,3],
            mode='lines',
            name="probe trajectory",
        ))
    fig.add_surface(x=spherex, y=spherey, z=spherez, opacity=1.0,showscale=False)

    # Speeding up animation, don't know if this works
    fig.layout.updatemenus[0].buttons[0].args[1]['frame']['duration'] = 10   # Speed of frame display
    fig.layout.updatemenus[0].buttons[0].args[1]['transition']['duration'] = 0 # Speed of motion between frames
    fig.update_scenes(aspectmode='cube') #Making the axes be a cube

    fig.show()