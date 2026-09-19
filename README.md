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
 

# Using This
Assuming all of the libraries have been installed, here is how this is used.
## Folder Structure
In the same directory as the scripts create these folders (case-sensitive) 
`Bodies`
`ImageOutputs`
`ProbePaths`
`Simulations`

## Generate Body Files
Run `SystemMovementCalculator.py` this script requires no modification

These "Body Files" contain the position of the planets at 1/60 second timesteps for the timespan of 0 - 23 minutes. These are stored in the `Bodies` folder as files `000.npy` to `015.npy` This also generates the `Properties.pkl` file which contains information about each planet such as mass, surface radius, etc.
## Specify Simulation Properties
This contains information only about running multiple simulations.
This section will only concern setting up `probe.py` 
The first variable to change would be `outputdir` set this to wherever the results of the simulation should be saved. I have them saved into a sub-directory of the simulations folder. `outputdir = "Simulations/0Front500-1500"`
If you want to change the Eye Shell radius, change `eye_distance`

Set NormalGravityforAll to True to make every planet have Newtonian gravity, set to false for game-accurate gravity.
`n_sim_per_pikmin` controls how many simulations happens before a results file is made and the process is killed. While it can be set to any number, I recommend something around 500-2000
`total_n_pikmin_to_make` This controls how many processes are made
The total number of launches made will be the product of `n_sim_per_pikmin` and `total_n_pikmin_to_make`. If you want 8,000,000 probe launches you can set `n_sim_per_pikmin` to 1,000 and `total_n_pikmin_to_make` to 8,000 For an analogy, imagine you are a teacher. You assign your class homework to do a certain number of practice problems out of the book. `n_sim_per_pikmin` Is how many practice problem each student will do, while `total_n_pikmin_to_make` is how many students you have in your class. And the total number of problems you must grade is how many probe launches are simulated.

`pikmin_on_field` Is how many processes will be made at once, can set to None and all available cores will be used, although your computer will be unable to do anything else until the simulation finishes. Unless you know how many cores you have I would just use None

If you want to change the conditions of the launch you must scroll down to line 575, and look for `pool.apply_async(simulationPikmin, args=(CannonIndex,None, Bodies, unitvec, 0, 1/60, 22, n_sim_per_pikmin, outputdir, True,500,1500))` above `#NOTE: This is where you change the settings for mass simulation mode` because of how multiprocessing is implemented its bit esoteric.
The current arguments (`pool.apply_async(simulationPikmin, args=(CannonIndex,None, Bodies, unitvec, 0, 1/60, 22, n_sim_per_pikmin, outputdir, True,500,1500))`) launching from CannonIndex, None for launch magnitude so it will be randomly selected, launch at time 0, go in timesteps of 1/60 seconds, end simulation at minute 22, do `n_sim_per_pikmin` simulations, output results to `outputdir` the minimum for the randomly selected launch magnitude is 500 m/s and the maximum is 1500 m/s

To instead make every probe launch at the same velocity (in this example 499m/s) relative to the probe cannon change the `None` to 499, so the line becomes  `pool.apply_async(simulationPikmin, args=(CannonIndex,499, Bodies, unitvec, 0, 1/60, 22, n_sim_per_pikmin, outputdir, True,500,1500))` the last two arguments aren't used anymore

With the configuration done, the script can now be run. The script will be closed at anytime with little data loss. If you only used one core, you set `n_sim_per_pikmin` to 1000, and 14,802 simulations were finished when you killed the python terminal, 14,000 simulations would be saved and the other 802 would be lost forever.
## Analyze Data
For `plot_spherical.py` just change `inputfolder` to whatever folder you're interested in viewing and run. When running subsequent analyses the image files will be overwritten so change the names or move them if you want to save them. The first tab that will pop up will show the fates of a probe launched at a certain direction on a unit sphere, relative to the Probe Cannon. The second tab will show the fates of probes launched at a certain direction after adding the velocity of the probe cannon.

The first tab shows what direction to launch if you were operating the Probe Cannon, and the second shows what direction to fire if the Probe Cannon was sitting stationary relative to the Sun.

The third tab shows the final position of a subset of probes. The largest sphere is the Eye Shell, the middle-sized sphere is the Sun and the smallest is Giant's Deep at the start of the loop.

For `live_analyzer.py` before doing anything such as `SortbyVisits` you should load in all the folders first to prevent any issues with differing numbers of columns. 
