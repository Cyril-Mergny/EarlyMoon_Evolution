# EarlyMoon_Evolution
Master research project by Cyril Mergny, supervised by Gabriel Tobie.
Thermodynamic and orbital evolution of the early Moon

1. Download the moon_algo folder
2. Open a terminal and cd to the moon_algo folder location
3. Run the command line " ./main.sh " on your terminal 

You may need to change permissions to run the main.sh file, to do so use the chmod "+rwx main.sh" command.

You can change the parameters of the simulation by modifying the input_parameters.in file. Most common changes are:

- n_layers: how much points is used to discretize each layer. In order (Core, Mantle, Upper Mantle, Magma ocean, Crust). Note that you do not need to specify the discretisation of the upper mantle and the crust as they are calculated from the discretisation of the magma ocean. The core and mantle have low interest in tidal heating so few points are used in these regions.

- the nbr of loops: how much time the algorithm will iterate over one cycle of lunar evolution

- the initial timestep: I suggest not to change this one because the timestep is adaptative anyway. 

- a,e the initial orbit of the Moon. If values are chosen to close to the Earth, overheating of the crust may occurs and cause flaws in the simulation.

If you have any question, do not hesitate to contact me.


