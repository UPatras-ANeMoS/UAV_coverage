# UAV Coverage under positioning uncertainty
Simulation for planar area coverage of a convex region by a swarm of aerial robots. The aerial robots are assumed to remain at a constant altitude while the positioning uncertainty only affects their `x` and `y` coordinates.

The robot positioning can be exact, in which case a Voronoi partitioning is used, or uncertain in which case a Guaranteed Voronoi partitioning is used instead. The implemented control laws are:
- Cell centroid
- r-limited cell centroid
- Free arcs (optimal for Voronoi)
- GV complete (optimal for Guaranteed Voronoi)
- GV compromise (simpler suboptimal for Guaranteed Voronoi)

The robots are also allowed to have finite communication ranges, in which case they create their cells using only the nodes inside their communication radius. There are also options to restrict the movement of robots inside the region of interest and to prevent them from colliding when their positioning is uncertain.

### Usage
Before using the scripts, you must copy the the `gpcmex` file from `MATLAB/[Release]/toolbox/map/map/private/` into `Functions/Polygon`. The sample script `copy_gpcmex.sh` does this for a typical Linux installation of MATLAB R2016b.

`SIM_coverage_disks.m` Use this to run simulations. You can select the control law used by uncommenting the appropriate line at the beginning of the file. You can also set other simulation parameters there.

`PLOT_sim.m` Used to plot the data from a `.mat` file saved after the end of a simulation. Set the `.mat` file to load and the plots to show at the beginning of the file.

`PLOT_compare_sims.m` Used to compare the results of two simulations. Used the same way as `PLOT_sim.m`. Note that the simulations must have the same duration and time step.

### Screenshots
<img src="./Screenshots/state_comm.png" width="49%"> <img src="./Screenshots/state_uncert.png" width="49%">
<img src="./Screenshots/objective.png" width="49%"> <img src="./Screenshots/trajectories.png" width="49%">

## Relevant Publications
[1] S. Papatheodorou, A. Tzes, and Y. Stergiopoulos, [*Collaborative Visual Area Coverage*](https://doi.org/10.1016/j.robot.2017.03.005), Robotics and Autonomous Systems, ISSN 0921-8890, Volume 92, June 2017, Pages 126–138, Elsevier

[2] S. Papatheodorou and A. Tzes, *Cooperative Visual Convex Area Coverage using a Tessellation-free strategy*, In Proceedings of the 56th IEEE Conference on Decision and Control (CDC) 2017, December 12-15, 2017, Melbourne, Australia [Accepted]

### License
Distributed under the [Apache License Version 2.0](LICENSE.txt)
<br>
Copyright © 2016-2017 Sotiris Papatheodorou
