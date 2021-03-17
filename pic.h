/*
 * File: pic.h
 * 
 * ------------
 * This class is the engine of the simulation. It is responsible for
 * initializing the data, running the simulation in a time loop and doing
 * closing up duties such as data dumps. The pic class initializes
 *  and runs a particle-in-cell simulation by declaring objects
 * and calling the methods defined in the boundary, particle, fields,  
 * collisions, e_collisions, ion_collisions, and graphs classes.
 *
 * 10/15
 * 
 * Date: 3/1/2021
 * 
 * Add functions to retreive graph averaged plot data and particle data
 * for use in graphics output displays.
 * 
 * Date: 1/28/2021
 * Description: Update file to account for a two-dimensional simulation. This 
 * file is pretty insular from what it calls and so the only change to 
 * the simulation runtime was the addition of a function call to generate
 * particle fields from grid fields and a few function calls to dump more
 * graph data to output files.
 * 
 * In addition the initialization of simulation data from user input was
 * added to the pic constructor.
 * 
 * The infrastructure is in place to account for more than one ion species
 * but this is not yet implemented in all classes.
 * 
 * Initialize, set up and destruct openCL variables in this file for other 
 * classes to use.
 * 
 */


#ifndef pic_h
#define pic_h

#define _CRT_SECURE_NO_WARNINGS 

#include "constants.h"
#include "parameters.h"
#include "boundary.h"
#include "particle.h"
#include "fields.h"
#include "collisions.h"
#include "e_collisions.h"
#include "ion_collisions.h"
#include "graphics.h"  // get opencl context
#include "graphs.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <string.h>
#include <cmath>  
#include <chrono>   /* use the high resolution clock */
#include "nr3.h"
//#ifdef WIN32
#include <windows.h>
//#else
//#include <unistd.h>
//#endif
#include <omp.h>

#include "CL/cl.h"
#include "CL/cl_platform.h"

#include "openCLutility.h"


using namespace std;

class pic{

   public:
  
   /*
    * Constructor: pic
    * Usage: pic pic_object;
    *        pic pic_object(true);
    *        pic pic_object(true,true,100000,90000,99999);
    * ------------------------------------------------------  
    */ 
    pic(bool restart=false,bool dump=true,int run_time=TOTAL_RUN_TIME,int start_graph=START_GRAPH_TIME,int end_graph=END_GRAPH_TIME);

   /* Destructor: ~pic */
   ~pic();

   /*
    * Method: pic_run
    * Usage:  pic_run();
    * -------------------
    * Implements the simulation cycle over the number of iterations defined by
    * the run_time variable.
    */
   void pic_run(void);

   /* The following methods set the values of private variables:
    *
    *   void set_dump(bool dump_in);
    *   void set_run_time(int run_time_in);
    *   void set_start_graph(int start_graph_in);
    *   void set_end_graph(int end_graph_in);
    *   void set_time_point1(int time_point);
    *   void set_time_point2(int time_point);
    *   void set_time_point3(int time_point);
    *
    */

   /*
    * Method: set_dump
    * Usage:  set_dump(false);
    * --------------------------
    * Pass false to turn off dumping simulation data into text files at the 
    * end of the simulation. The default is true. These files can be used
    * to start a new simulation where the previous one left off.
    */
   void set_dump(bool dump_in);

   /*
    * Method: set_run_time
    * Usage:  set_run_time(number_of_iterations); 
    * --------------------------------------------
    * An alternative way to set the number of iterations that the simulation 
    * steps through. 
    */
   void set_run_time(int run_time_in);

   /*
    * Method: set_start_graph
    * Usage:  set_start_graph(time_to_start_summing_av_data);
    * --------------------------------------------------------
    * An alternative way to set the iteration number from which to start
    * summing average data for the graphs class.
    */
   void set_start_graph(int start_graph_in);

   /*
    * Method: set_end_graph
    * Usage:  set_end_graph(time_to_stop_summing_av_data);
    * -----------------------------------------------------
    * An alternative way to set the iteration number to stop summing
    * average data for the graphs class.
    */
   void set_end_graph(int end_graph_in);

   /*
    * Method: set_time_point1, set_time_point2, set_time_point3
    * Usage:  set_time_point1(cell_position1),set_time_point2(cell_position2),
    *         set_time_point3(cell_position3)
    * -------------------------------------------------------------------------
    * The three cells at which to dump the variation in time of several
    * variables including the potential and charge density.
    */
   void set_time_point1(int time_point);
   void set_time_point2(int time_point);
   void set_time_point3(int time_point);

   /*
   * Dump the current state of each particle type to a separate output file.
   * The contents of the struct particleVec for each particle is output.
   */
   void dump_data(void);

   /* initialize openCL parameters and openCL parameters needed by other classes */
   void initialize_OpenCL(void);

   /*
    * Set the number of timesteps to run the simulation.
    */
   void setRuntime(int runtime);

   /*
    * Used to retrieve particle data for output plots by the graphics output
    * functions.
    */
   vector<particleVec> getElectrons(void);
   vector<particleVec> getIons(void);

   /*
    * Used to retrieve averaged graph class data for output plots rendered  
    * in the graphics output functions.
    */
   double** get_phiPlot(void);  // electrostatic potential
   double getPhiPlotMax(void);  // max value of the electrostatic potential

   double** get_ExPlot(void); // x component of the electric field
   double getExPlotMax(void); // max value 

   double** get_EyPlot(void); // y component of the electric field
   double getEyPlotMax(void); // max value 

   double** get_eQPlot(void); // electron density
   double get_eQPlotMax(void); // max value 

   double** get_ionQPlot(void); // ion density
   double get_ionQPlotMax(void); // max value 

   /* True if want to write data to files. Default is true for
   * batch mode and is always false for interactive mode.
   */
   void setWriteGraphDataToFiles(bool setFlag);

   /*
    * Initialize the number of timesteps over which to take
    * averaged data that is displayed in the graphics output
    * window.
    */
   void initializeGraphsForInteractive(int interval);

   private:

   particle electrons,*ions;   // for calculating particle quantities
   fields grid_fields;        // for calculating field quantites
   e_collisions *e_collision;  // for electron collisions
   ion_collisions *ion_collision;  // for ion collisions
   graphs graphs_fields,graphs_electrons,*graphs_ions;  // for graph data dump outputs
   boundary sec_e, *sec_ions; // for boundary calculations

   //OpenCL variables
   cl_context context;  // holds the OpenCL context
   cl_command_queue cmdQueue; //holds the current Queue
   cl_program program;  // used to build a program with the kernel function
   cl_kernel kernelCC;  //holds the OpenCL Kernel
   cl_kernel kernelE;  //holds the OpenCL Kernel
   cl_kernel kernelV;  //holds the OpenCL Kernel 
   cl_kernel kernelVg; //holds the OpenCL Kernel 
   cl_kernel kernelVd; //holds the OpenCL Kernel 

   // Memory where variable arrays are stored on the GPU
   cl_mem dP;  // for particle arrays
   cl_mem dExp;   // for Ex at the particles
   cl_mem dEyp;   //for Ey at the particles
   cl_mem dQ;  // for the electric charge on the grid
   cl_mem dEx; // for Ex on the grid
   cl_mem dEy; // for Ey on the grid
   cl_mem dCountpX; // number of particles per cell where y = length/2
   cl_mem dCountpY; // number of particles per cell where x = length/2
   cl_mem dv_grid_tempX; //holds temp values for summing average velocities
   cl_mem dv_grid_tempY; //holds temp values for summing average velocities
   cl_mem dtemp_distX1; //holds temp values for summing average vx distribution
   cl_mem dtemp_distX2; //holds temp values for summing average vx distribution
   cl_mem dtemp_distY1; //holds temp values for summing average vy distribution
   cl_mem dtemp_distY2; //holds temp values for summing average vy distribution

   int run_time,start_graph,end_graph;  // start and end times for capturing data for graph dumps
   int no_ions_inject[NO_ION_SPECIES];  // number of ions to inject per time step

   bool writeGraphDataToFiles;  //by default is true

   // initialize particle data with user input
   void initialize_ion_particles (void);
   void initialize_electron_particles(void);

   // initialize collision data with user input
   void initialize_ion_collision(void);
   void initialize_electron_collision(void);

   // initialize boundary data with user input
   void initialize_boundary(void);

   // initialize graph data with user input
   void initialize_graph_parameters(void);

   

  /*
   * Variables used to store the variation of select data in time, including
   * the total particle_number for each species, electric_potential at three
   * poistions, and the particle charge density of each species at position3.
   */
   int *time_plot,*no_e,**no_i; //number of electrons/ions at one position 
   // potential, and electron and ion densities at indicated positions
   // on the grid
   double *phi_mid,*phi_left,*phi_right,*eden_right,**iden_right;

   bool haveCollisions; // flag that is true if want collisions
  
   ofstream timedata;   // output some of the physical quantities in time

   bool dump,restart;   // true if want to start simulation from a previous dump 

   cl_int status;  // returned status from OpenCL calls which is tested against CL_SUCCESS
   cl_device_id device; // the GPU device ID

};

#endif /*pic*/
