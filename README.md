# Outerwilds Eye Finder

Spoilers for the game Outerwilds, also I'm pretty this is only going to really make sense if you've played it.

What this does is simulates launching the probe at a specific time, velocity, and direction. It records whether anything is hit, what bodies are visited, if it travelled a specific distance from the Sun, where it reached that specific position from the Sun, and the final position at the end of the loop.

Atmospheric drag is simulated, and there is the option to either use true Newtonian gravity, or use in-game gravity rules. 

# Brief Overview of Relevant Files
Right now this place is a bit of a mess, here are the files of interest.

- `OrbitalGeometry.py` - Library that makes possible to define planets
- `SystemMovementCalculator.py` - This is where the planets are defined, their positions at specific times are calculated and output to a file for the probe simulation to use. Also allows for visualization of the system using plotly
- `probe.py` - Simulates probe launches, currently uses python Multithreading to run in parallel. Results are saved within a folder, currently it takes around 16-90 seconds to run 250 simulations
- `plot_spherical.py` - Takes the results file and processes them. It tallies how many times a body was visited by a probe. It also provides a 3D visualization of the ending state of a subset of probe. It also creates a cube map and equirectangular scatter plot of where probes reached the "Eye Shell" which is a certain distance from the Sun where the Eye is known to orbit. This provides information on how many measurements a region of the sky gets. 
