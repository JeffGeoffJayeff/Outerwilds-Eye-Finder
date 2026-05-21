import numpy as np
import math
import pandas as pd
from pathlib import Path
import scipy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import plotly.express as px
import plotly.graph_objects as go
import probe

## Purpose: This one will have many probe simulations in it instead of one
Properties = pd.read_pickle("Properties.pkl")
bodiesfolder = Path("Bodies")
files = list(bodiesfolder.glob("*.npy"))
G = 10**-3
eye_distance = 286500 #Distance of the eye from the sun in meters https://www.reddit.com/r/outerwilds/comments/t7mxcy/how_far_away_is_the_eye_base_game_spoilers/
sunBodyIndex = 0 #Index that is the Sun in the Bodies list
NormalGravityforAll = True #This controls whether gravity is calculated using Newtonian gravity, or if it uses the so called linear gravity https://www.youtube.com/watch?v=dpKUoWgRBSU
# If True then the mass for each planet is changed to produce the same gravity at the surface in both systems
n_sim = 1

Bodies = [] #Create list to store bodies into
Names = []
for i in range(0,len(files)): #Load in bodies
    Bodies.append(probe.Body(np.load(files[i]),Properties.iloc[i]))
    Names.append(Bodies[i].name)
    if Bodies[i].name == "Cannon":
        Cannon = Bodies[i]
events = [probe.make_visit_event(body) for body in Bodies] #Make visiting events
events.append(probe.eyeDistance) #Add event for reaching the distance that the eye is from the Sun
events.append(probe.hitBody) #Add hitting event

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

results = np.empty(
    n_sim,
    dtype = [
        ("Reached Eye",np.bool_),
        ("Sun Visits",np.int16),
        ("Sun Station Visits",np.int16),
        ("Ember Twin Visits",np.int16),
        ("Ash Twin Visits", np.int16)
    ]
)