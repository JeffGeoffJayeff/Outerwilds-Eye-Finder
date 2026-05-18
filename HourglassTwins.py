import OrbitalGeometry as OG
import math
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

G = 1*10**(-3) #Gravitational Constant
sun_mass = 4*10**11
Stepsize = 1 #In seconds
EndMinute = 1 

Twin_mass = 500000

FocalBody = OG.planet(50,0,sun_mass,name="Focal Body")
FocalBody.useStartingPoint(10,0)
FocalBody.isGravityLinear = True
FocalBody.mass = 800000
FocalBody.createDataFrame(endminute=EndMinute,stepsize=Stepsize)

CaveTwin = OG.planet(250,0,Twin_mass,name="Ember Twin",parent=FocalBody)
CaveTwin.Period = 125 #The period number here doesn't matter
CaveTwin.useStartingPoint(-204.788,-143.394)
CaveTwin.createDataFrame(endminute=EndMinute,stepsize=Stepsize)

TowerTwin = OG.planet(250,0,Twin_mass,name="Ash Twin",parent=FocalBody)
TowerTwin.Period = 125
TowerTwin.useStartingPoint(204.788,143.394)
TowerTwin.createDataFrame(endminute=EndMinute,stepsize=Stepsize)



range = [-1000,1000]
df = pd.concat([CaveTwin.Dataframe,TowerTwin.Dataframe])
fig = px.scatter_3d(df,x="x",y="y",z="z",animation_frame="time",color="body",range_x=range,range_y=range,range_z=range,color_discrete_sequence=["#f05f44","#5a8240","#647297","#31a174","#55503a","#349fb9"])

# Speeding up animation, don't know if this works
fig.layout.updatemenus[0].buttons[0].args[1]['frame']['duration'] = 10   # Speed of frame display
fig.layout.updatemenus[0].buttons[0].args[1]['transition']['duration'] = 0 # Speed of motion between frames
fig.update_scenes(aspectmode='cube') #Making the axes be a cube
fig.show()

