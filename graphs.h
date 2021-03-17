/*
 * File: graphs.h
 * ---------------
 * 
 * 
 * This file defines the graphs class which stores variables and methods
 * that make text files of diagnostic data including the average fields, 
 * particle velocites, densities, and particle velocity distributions over
 * a period of time selected by the user.
 *
 * 4/15
 * 
 * Date: 2/28/2021
 * 
 * Functions that calculate the average velocity per cell and the velocity
 * distributions have been vectorized using OpenCL
 * 
 * Date: 2/6/2021
 * 
 * Altered some of the functions to gather two-dimensional data and added
 * extra functions for extra output files with two-dimensional data suitable
 * for contour and surface plots. Current output includes
 * 
 * -A slice of the fields (net charge density, potential, electric field in X) 
 * along a line in x halfway up the grid in y.
 * 
 * -A slice of the fields (net charge density, potential, electric field in Y) 
 * along a line in y halfway along the grid in x.
 * 
 * - Two dimensional output of the fields: net charge density, elecron charge denisty,
 *  ion charge density, potential, electric field in X, electric field in Y
 * 
 * - The velocity distribution in X at two locations on the grid
 * 
 * - The velocity distribution in Y at two locations on the grid
 * 
 * - The average ion and electron velocities along a line in x halfway up the grid
 *  in y.
 * 
 *  - The average ion and electron velocities along a line in y halfway up the grid
 *  in x.
 */
#define _CRT_SECURE_NO_WARNINGS   //Visual Studio Uses to Ignore VS related warnings

#include "constants.h"
#include "parameters.h"
#include "particleVec.h" 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <limits>


//#ifdef WIN32
#include <windows.h>
//#else
//#include <unistd.h>
//#endif
#include <omp.h>

#include "CL/cl.h"
#include "CL/cl_platform.h"

#include "mglin.h"

#include "openCLutility.h"

using namespace std;

class graphs{
 
  public:

   /*
    * Constructor: graphs
    * Usage: graphs graph_object;
    *        graphs graph_object(ION_TYPE,SPREAD,VT,QCI);
    */       
    graphs(char *type_inp=ION_TYPE,int spread_in=SPREAD,double vt_in=VTI,double qc_in=QCI);
    
   /* Destructor: ~graphs */
    ~graphs();

   /* These methods sum the indicated quantities over a number of time-steps
    * defined by the user:
    *
    * void sum_fields(double *rho,double *phi,double *Egrid);
    * void sum_density(double *q);
    * void sum_av_velocities(vector<particleVec> &particleObj);
    * void sum_av_dist1(vector<particleVec> &particleObj);
    * void sum_av_dist2(vector<particleVec> &particleObj);
    */

   /*
    * Method: sum_fields
    * Usage:  sum_fields(net_charge_density,electric_potential,electric_field);
    * -------------------------------------------------------------------------
    * Sums field quantities, including the net_charge density, electric 
    * potential and the electric field.
    */
    void sum_fields(double **rho, MatDoub_IO phi,double *EgridX, double *EgridY);

   /*
    * Method: sum_density
    * Usage:  sume_density(particle_charge_density);
    * -----------------------------------------------
    * Sums the charge density of each species.
    */
    void sum_density(cl_double *q);

   /* 
    * Method: sum_av_velocities
    * Usage:  sum_av_velocities(particle_vector);
    * --------------------------------------------
    * Sums the average velocity of each grid cell for the species passed to 
    * the arguement.
    */
    void sum_av_velocities(vector<particleVec> &particleObj);

   /*
    * Methods: sum_av_dist1, sum_av_dist2
    * Usage:   sum_av_dist1(particle_vector), sum_av_dist2(particle_vector);
    * -----------------------------------------------------------------------
    * Sums the particle velocity distribution at cell locations defined by the
    * dist_location and dist_location2 variables. Particle velocities are 
    * included from several cells around each location indicated by the 
    * spread variable.
    */ 
    void sum_av_dist(vector<particleVec> &particleObj);
  

   /*
    * The following functions take averages of the indicated quantity by
    * dividing it over the number of steps each item was summed up over. 
    * 
    * void average_fields(void);
    * void average_density(void);
    * void average_velocities(void);
    * void average_dist1(void);
    * void average_dist2(void);
    *
    */

    void average_fields(void);
    void average_velocities(void);
    void average_dist1(void);
    void average_dist2(void);
    void average_density(void);

   /* 
    * The following methods write the average of the indicated quantity to
    * a text file.
    *
    * void write_fields_to_file(void);
    * void write_density_to_file(void);
    * void write_v_to_file(void);
    * void write_dist1_to_file(void);
    * void write_dist2_to_file(void);
    *
    */

    void write_fields_to_file(void);
    void write_density_to_file(void);
    void write_v_to_file(void);
    void write_dist1_to_file(void);
    void write_dist2_to_file(void);

    void write_rho_to_file(void);
    void write_phi_to_file(void);
    void write_EgridX_to_file(void);
    void write_EgridY_to_file(void);

   /*
    * A group of methods that sets the values of private variables.
    *
    * void set_dist_location(int location);
    * void set_dist_location2(int location);
    * void set_type_in(char* str);
    * void set_spread(int val);
    * void set_vt(double vt_in);
    * void set_qc(double qc_in);
    *
    */

   /* 
    * Method: set_dist_location
    * Usage:  set_dist_location(cell_numberY, cell_numberX);
    * ----------------------------------------   
    * Pick the central cell from which particle velocities are picked to
    * determine the particle velocity distribution.
    */
    void set_dist_location(int cell_numberY, int cell_numberX);

   /* 
    * Method: set_dist_location2
    * Usage:  set_dist_location2(cell_number);
    * ----------------------------------------   
    * Pick a second central cell from which particle velocities are picked to
    * determine the particle velocity distribution.
    */
    void set_dist_location2(int cell_numberY, int cell_numberX);

   /*
    * Method: set_spread
    * Usage:  set_spread(number_of_cells);
    * -------------------------------------
    * Set the number of cells around the central one from which to pick
    * particle velocities from. These are used to calculate the particle
    * velocity distribution.
    */
    void set_spread(int val);
 
   /*
    * Method: set_type_in
    * Usage:  set_type_in("i1");
    * ----------------------------
    * Define a string that is used to distinguish which species each diagnostic
    * text file is reffering to.
    */
    void set_type_in(char* str);
   
   /*
    * Method: set_vt
    * Usage:  set_vt(thermal_velocity);
    * ----------------------------------
    * Set the value of the thermal velocity used to normalize the particle
    * velocity distribution.
    */
    void set_vt(double vt_in);

   /*
    * Method: set_qc
    * Usage:  set_qc(charge_of_macro_particle);
    * -------------------------------------------
    * Use to set the the value of the charge of each macro-particle that is
    * multiplied by the particle number density of each cell. To plot the number
    * density set this value to 1.0.
    */
    void set_qc(double qc_in);

    /* Sets parameters needed for openCL. The parameters are initialized
     * and destroyed in the PIC class and are sent here for use.
     */
    void init_openCL(cl_context context, cl_command_queue cmdQueue,
                    cl_kernel kernelVg, cl_kernel kernelVd,cl_mem dCountpX,
                    cl_mem dCountpY, cl_mem dv_grid_tempX, cl_mem dv_grid_tempY,
                    cl_mem dtemp_distX1, cl_mem dtemp_distX2,
                    cl_mem dtemp_distY1, cl_mem dtemp_distY2);

    /* Sets parameters needed for openCL for the case of a steadystate
     * boundary. In this instance the particle arrays don't change size
     * therefore memory need only be reserved during initialization.
     */
    void init_openCL_steady(cl_mem dP);

    /*
    * Used to retrieve averaged graph class data for output plots rendered  
    * in the graphics output functions.
    */
    double** get_phiPlot(void); // electrostatic potential
    double getPhiMax(void); // max value of the electrostatic potential
    double** getExPlot(void); // x component of the electric field
    double** getEyPlot(void); // y component of the electric field
    double getExMax(void); // max value of Ex
    double getEyMax(void); // max value of Ey
    double** getQPlot(void); // density
    double getQMax(void); // max value of density

    /* 
     * Set the number of steps over which to take averages for graph
     * output.
     */
    void set_no_graph_steps(int numSteps);

    private:

    // the location and distribution spread for velocity distribution plots
    int spread;
    int dist_locationX,dist_location2X;
    int dist_locationY,dist_location2Y;
    int no_graph_steps;

    double **rho_plot;  // for plotting the net charge density 
    double **q_plot;  // for plotting the electron or ion charge density
    double **phi_plot;  // for plotting the potential
    double **Egrid_plotX;   // for plotting the electric field in the x direction
    double **Egrid_plotY; // for plotting the electric field in the y direction

    double phiMax; // Maximum phi plot value for scaling
    double ExMax; // Maximum Ex plot value for scaling
    double EyMax; // Maximum Ey plot value for scaling
    double qMax; // Maximum density plot value for scaling

    // for plotting the electron or ion particle x-velocity distribution at location 1
    double *d_plot1X;

    // for plotting the electron or ion particle x-velocity distribution at location 2
    double *d_plot2X;

    // for plotting the electron or ion particle y-velocity distribution at location 1
    double *d_plot1Y;

    // for plotting the electron or ion particle y-velocity distribution at location 2
    double *d_plot2Y;

    // for plotting the average electron or ion x-velocity per cell along a slice of 
    // the grid in y
    double *v_plotX;

    // for plotting the average electron or ion y-velocity per cell along a slice of 
    // the grid in x
    double *v_plotY;

    // for holding temporary data when calculating the above velocity distributions
    cl_double *temp_distX;
    cl_double *temp_distY;
    cl_double *temp_distX2;
    cl_double *temp_distY2;

    // holds the number of particles per cell for the above average x-velocity plots
    cl_double *count_particlesX; 

    // holds the number of particles per cell for the above average y-velocity plots
    cl_double *count_particlesY; 
    //holds temp values for summing average velocities
    cl_double *v_grid_tempX;
    cl_double *v_grid_tempY;

    //OpenCL variables
    cl_int status; // holds the OpenCL status of OpenCL calls
    cl_device_id device;  // indicates the GPU device being used
    cl_context context;  // holds the OpenCL context
    cl_command_queue cmdQueue; //holds the current Queue
    cl_kernel kernelVg;  //holds the OpenCL Kernel 
    cl_kernel kernelVd;  //holds the OpenCL Kernel 

    // Memory where variable arrays are stored on the GPU
    cl_mem dP; // for particle arrays
    cl_mem dCountpX; // number of particles per cell where y = length/2
    cl_mem dCountpY; // number of particles per cell where x = length/2
    cl_mem dv_grid_tempX; //holds temp values for summing average velocities
    cl_mem dv_grid_tempY; //holds temp values for summing average velocities
    cl_mem dtemp_distX1; //holds temp values for summing average vx distribution
    cl_mem dtemp_distX2; //holds temp values for summing average vx distribution
    cl_mem dtemp_distY1; //holds temp values for summing average vy distribution
    cl_mem dtemp_distY2; //holds temp values for summing average vy distribution
    char *type_in;  // indicates if an electron or ion

    double dx;  // cell size
    double vt;  // thermal velocity
    double qc;  // macro charge

};
