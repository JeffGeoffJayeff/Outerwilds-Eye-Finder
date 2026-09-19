# Outerwilds Eye Finder

Spoilers for the game Outerwilds, also I'm pretty this is only going to really make sense if you've played it.

What this does is simulates launching the probe at a specific time, velocity, and direction. It records whether anything is hit, what bodies are visited, if it traveled a specific distance from the Sun, where it reached that specific position from the Sun, and the final position at the end of the loop.

Atmospheric drag is simulated, and there is the option to either use true Newtonian gravity, or use in-game gravity rules. [More information about what that means](https://www.youtube.com/watch?v=dpKUoWgRBSU)

# Brief Overview of Relevant Files
Right now this place is a bit of a mess, here are the files of interest.

- `OrbitalGeometry.py` - Library that makes possible to define planets
- `SystemMovementCalculator.py` - This is where the planets are defined, their positions at specific times are calculated and output to a file for the probe simulation to use. Also allows for visualization of the system using plotly
- `probe.py` - Simulates probe launches, currently uses python Multithreading to run in parallel. Results are saved within a designated folder, currently it takes around 60 seconds to run 1,000 simulations in a single processes, allowing for 14 processes allows for ~14,000 launches per minute. 
- `plot_spherical.py` - Takes a results folder and processes files inside. It tallies how many times a body was visited by a probe. It provides a 3D visualization of the ending state of a random subset of ~100,000 probes. A map of where all probes in the folder intersected with the Eye Shell (if they did) is made. In both equirectangular and cube-map form.
- `live_analyzer.py` - Allows for analysis of probe results in a terminal format. With the following commands 
  -  `FolderLoad` - Loads in the results of all files in a folder, notably this command can be run multiple times in a session to load data from multiple folders
  -  `VisitStats` - Provides statistics about how many times each body was visited/hit
  -  `SortbyVisits` - Sorts the loaded probes in descending order by any of the following: total visits, total visits excluding visits to the Sun (nosun), total visits excluding visits to the Sun and Giant's Deep (interesting), and by visits to a specific body
  -  `LookupbyID` - Prints out all of the information of a launch, such as launch direction, launch velocity, bodies visited, etc. Accepts either the index in the giant list of all probes, or the UUID of a launch
  -  `LaunchbyID` - Simulates a specific launch, using either an Index or UUID and then shows a plot of the trajectory, putting a True after the UUID or Index also shows the Planets
  -  `CreatePathbyID` - Simulates a specific launch, using either an Index or UUID, and just saves the result to the `ProbePaths` folder, does not show the trajectory on a plot
  -  `PlotAllPaths` - Plots the trajectory of every probe in `ProbePaths`, adding a True also shows the paths of every planet
 

