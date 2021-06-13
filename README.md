#   Two-Dimensional Vectorized Plasma Particle-in-Cell Simulation Code With Hardware Accelerated Graphics Output

A two-dimensional plasma particle-in-cell simulation code in C/C++ has been developed and vectorized onto the graphics-processing-unit (GPU) with hardware-accelerated graphics output.  

The simulation can be run either interactively where real-time graphics showing the current state of the simulation is rendered using OpenGL or in batch mode where an extensive list of averaged data is calculated and output to text files for import into applications such as Matlab.  Additionally, datafiles saving the current state of the simulation can be dumped at the end of a run and optionally read in at the start of a later simulation to pick up where the previous simulation left off. This feature is available in both batch and interactive mode whereby the user has the flexibility of running a simulation in batch mode up to a time point of interest and then restarting from that state in graphics mode for real-time dynamic observation.

All sources referenced during code development are located as pdfs in the *References* folder.

![Demonstration](https://raw.githubusercontent.com/SelmaLeathem/Vectorized-Particle-in-Cell-Code-With-Graphics/main/PIC_Simulation2.gif)


#### *See the following video for an overview of the project and a demonstration:* https://media.oregonstate.edu/media/t/1_stl8aycr


## Instructions

###  i.	Requirements
*	Windows Operating System (preferably Windows 10)
*	A graphics processing unit (GPU).
*	Visual Studio 2019 
*	The ability to run OpenCL on the GPU.

The code has been tested on the following systems.

1)	Windows 10 with the following GPU:
Name    = 'NVIDIA CUDA'
	Vendor  = 'NVIDIA Corporation'
	Version = 'OpenCL 1.2 CUDA 11.1.114'
	Profile = 'FULL_PROFILE'
	Number of Devices = 1
	Device #0:
		Type = 0x0004 = CL_DEVICE_TYPE_GPU
		Device Vendor ID = 0x10de (NVIDIA)
		Device Maximum Compute Units = 28
		Device Maximum Work Item Dimensions = 3
		Device Maximum Work Item Sizes = 1024 x 1024 x 64
		Device Maximum Work Group Size = 1024
		Device Maximum Clock Frequency = 1582 MHz

2)	Windows 10 with the following GPU:
Name    = 'Intel(R) OpenCL HD Graphics'
	Vendor  = 'Intel(R) Corporation'
	Version = 'OpenCL 2.1 '
	Profile = 'FULL_PROFILE'
	Number of Devices = 1
	Device #0:
		Type = 0x0004 = CL_DEVICE_TYPE_GPU
		Device Vendor ID = 0x8086 (Intel)
		Device Maximum Compute Units = 24
		Device Maximum Work Item Dimensions = 3
		Device Maximum Work Item Sizes = 256 x 256 x 256
		Device Maximum Work Group Size = 256
		Device Maximum Clock Frequency = 1200 MHz

###  ii.	Start the Simulation

1.	Right Click on the file Pic2D.sln and select open with -> Visual Studio 2019
1.	Select Build -> Build Solution
1.	Select Debug -> Start Without Debugging

*Note: There are example files in place so the simulation can be started out of the box.*

*If altering the user input file parameters.h (see section v) in between runs it might be necessary to additionally* 

	*Select Build -> Clean Solution* 

*before building the solution in step 2.*

###  iii.	Batch Mode

To run the simulation in batch mode set the variable INTERACTIVE to “NO” in the *parameters.h* file (more information on the *parameters.h* file is given in section v)  and then follow the directions in section ii above to compile and run the simulation.

When running in batch mode the current timestep number and number of particles of each species type is printed to the console. This information is useful if you are using natural boundaries and sources such as ionization and want to monitor the particle number. In addition, for a steady state boundary (explained in section v) where an electron-ion pair is *“plopped”* randomly into the simulation domain every time an ion arrives at the boundary the number of pairs *“plopped”* is printed to the console.

In batch mode averaged graph data is calculated and output to text files for import into external applications such as Excel or Matlab. The user can control when and for how long the averages are taken and more by setting values in the *parameters.h* file as is explained further in section v. 

**Graph data output files in batch mode includes:**

*  **Av_velocity_e.txt:** The average electron velocity Vx per cell taken along X where Y = Length/2 and the average electron velocity Vy per cell taken along Y where X = Length/2.
*  **Av_velocity_i.txt:** The average electron velocity Vx per cell taken along X where Y = Length/2 and the average electron velocity Vy per cell taken along Y where X = Length/2.
*  **Density2d_e.txt:** The average electron charge density on the two-dimensional grid.
*  **Density2d_i1.txt:** The average ion charge density on the two-dimensional grid.
*  **Egrid2dX.txt:** The average X component of the electric field on the two-dimensional grid.
*  **Egrid2dY.txt:** The average Y component of the electric field on the two-dimensional grid.
*  **First_v_distX_e.txt:** The average electron Vx velocity distribution at location 1.
*  **First_v_distX_i1.txt:** The average ion Vx velocity distribution at location 1.
*  **First_v_distY_e.txt:** The average electron Vy velocity distribution at location 1.
*  **First_v_distY_i1.txt:** The average ion Vy velocity distribution at location 1.
*  **Second_v_distX_e.txt:** The average electron Vx velocity distribution at location 2.
*  **Second_v_distX_i1.txt:** The average ion Vx velocity distribution at location 2.
*  **Second_v_distY_e.txt:** The average electron Vy velocity distribution at location 2.
*  **Second_v_distY_i1.txt:** The average ion Vy velocity distribution at location 2.
*  **Phi2d.txt:** The average electrostatic potential on the two-dimensional grid.
*  **Rho_phi_EgridX.txt:** The average net charge density, electrostatic potential and X component of the electric field along X at Y = Length/2.
*  **Rho_phi_EgridY.txt:** The average net charge density, electrostatic potential and Y component of the electric field along Y at X = Length/2.
*  **Rho2d.txt:** The average net charge density on the two-dimensional grid.
*  **timeData.txt:** The following points in time during the simulation: Total number of electrons, total number of ions, electrostatic potential at the left boundary, electrostatic potential in the center of the simulation, electrostatic potential at the right boundary, electron density at the right boundary, ion density at the right boundary

###  iv.	Interactive (Graphics) Mode

To run the simulation in interactive mode set the variable INTERACTIVE to “YES” in the *parameters.h* file (more information on the *parameters.h* file is given in section v) and then follow the directions in section ii above to compile and run the simulation. When the simulation starts the graphics output window and graphical user interface appear.

**Start**: Starts the simulation
**Pause**: Pauses the simulation. Press Pause or Start to start the simulation again.
**Quit**: Stops the simulation and performs a complete exit.
**Dump Data When Quit**: Select this checkbox to preserve the current state of the simulation and dump it to a text file when quitting.

*Note, if starting simulation from where a previous run left off, ensure the following dump files from the previous run are present in the solution folder:*

      dataDump_e.txt
      dataDump_i1.txt

**Rotate Graph**: Select the rotate wheel with the mouse and move to rotate two-dimensional field plots. *Note, graphs can also be rotated by selecting the plot window and dragging the mouse in the desired rotational direction.*

**Choosing a Graph:**

Select a graph from the list by selecting the radio button beside the graph. Graphs include:

*	Particle X- Y plot
*	Electron Vx-Vy plot, where V refers to the particle’s velocity
*	Ion Vx-Vy plot
*	The electrostatic potential on the two-dimensional grid.
*	The X component of the electric field on the two-dimensional grid.
*	The Y component of the electric field on the two-dimensional grid.
*	The electron charge density on the two-dimensional grid.
*	The ion charge density on the two-dimensional grid.

###   v.	User Input File

The file *parameters.h* is where users may select and alter input options. The user may specify which mode to run the simulation in (batch or interactive) or whether to start the simulation from a previous run and much more…

*Note, if starting simulation from where a previous run left off, ensure the following dump files from the previous run are present in the solution folder:*

dataDump_e.txt
dataDump_i1.txt

The user has significant flexibility with the physics accounted for in the simulations through the *parameter.h* file. For example, collisions can be turned on or off and secondary emission and reflection at the boundary can optionally be accounted for and the user can select the total run time for batch mode. All options are fully explained in the *parameters.h* file and I would like to highlight some of the more salient features:

**Collision Options Include:**

*	Electron impact ionization
*	Electron-atom scattering
*	Ion-Atom charge exchange

**Boundary Options Include:**

*	Reflection
*	Absorption
*	Back Scattering 
*	Electron and ion impact secondary electron emission

There is a special steady-state boundary where collisions and other boundary effects are turned off and an electron-ion pair is placed randomly in the simulation domain when an ion reaches a surface. This results in a rapid and consistent steady state.

**Field Options Include:**

A magnetic field at an angle on the X-Y plane

**Species Options Include:**

*	Argon 
*	Hydrogen

**General Simulation Options Include:**

*	Number of grid points
*	The cell size
*	The length of the domain
*	The time step
*	The microcharge size per particle
*	The value of the voltage at each boundary. For the right boundary can optionally select a period for a sinusoidal variation.

**Key Options in Batch Mode Include:**
*	Total time to run the simulation
*	How long to take averages over for graph output data
*	When to start taking the averages
*	Where abouts on the two-dimensional grid to take velocity distributions and how many cells should these distributions account for

*Note, this is not an exhaustive list. For all options refer to the parameters.h file.*



###  vi.	Examples

Three example *parameter.h* files are in the folder’s:

1.	“Parameter file for Ar steadystate”
This file starts the simulation from a previous run and employs a rapid steady-state boundary where an electron-ion pair are placed randomly in the simulation domain when an ion reaches a surface. The dump-files from the previous run are included in this folder.
1.	“Parameter file for Ar with coll-boundary”
This file starts the simulation from time zero and employs all collisions types and boundary effects for an Argon gas.
1.	“Parameter file for H with coll-boundary”
This file starts the simulation from time zero and employs all collisions types and boundary effects for a Hydrogen gas.

Example of simulation graph output in batch Mode is in the folder: *Example Output for Batch Mode*, where a steady-state boundary was employed. 


