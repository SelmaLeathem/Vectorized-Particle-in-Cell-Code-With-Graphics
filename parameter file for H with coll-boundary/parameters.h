/*
 * File: parameters.h
 * ------------------  
 * 
 * Date: 2/26/2021
 * Description: Holds all user defined input parameters for the particle-in-cell simulation.
 *              These values are not currently validated anywhere. 
 * 
 *              In addition to defining physical quantites these parameters define the
 *              boundaries of the simulation such as whether to include collisions, boundary
 *              effects, should cross sections be constant or energy dependent, should there
 *              be additional particle injection, where/if to take graph data for output etc...
 * 
 *              Use this file to set whether running in interactive (graphics) mode or batch mode.
 * 
 *              Currently, can choose between Argon and Hydrogen gases.
 * 
 *              See Plasma Physics Via Computer Simulation by Birdsall and Langdon for more
 *              information on how to define input parameters.
 * 
 */
#define _CRT_SECURE_NO_WARNINGS 

#ifndef parameters_h
#define parameters_h

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <string.h>
#include <cmath>

#include "constants.h"

using namespace std;


/* Random seed generator. Supply an integer to repeat a simulation. */
#define RAND_SEED time(NULL)

/* Select yes if running in interactive mode with graphics, otherwise choose "NO" */
#define INTERACTIVE "YES"  // "NO" // "YES"

/* At the end of each batch mode simulation and optionally for the interactive mode 
 * the user can dump data text files. If the user selects RESTART = true then instead
 * of starting from time zero the simulation picks up where the data in the dump files
 * left off. This option can be used in interactive and batch modes. 
 */
#define RESTART false  //true // false

/* Create an array of neutral velocities to randomly pick from during a collision.
 * The following parameter determines the size of this array.
 */
#define INIT_NO_NEUTRALS 10000

/*Number of ions and electrons to optionally inject uniformly throughout the domain, 
* per time-step. 
*/
#define NUMBER_TO_INJECT 10  

/* This is not required by the simulation directly. It is used to calculate later
 * quantities that are passed to the simulation:
 * 
 * Electron particle density is used to calculate electron plasma frequency
 * and the Debye length. This density is the order of magnitude that I expect
 * the simulation to attain after ionization. 
 */
#define ELECTRON_PARTICLE_DENSITY (double)1.0e15
#define ION_PARTICLE_DENSITY (double)1.0e15

/* Used directly by the simulation. Define the NEUTRAL_DENSITY. */
#define GAS_PRESSURE ((double)0.2)*((double)133.32) //Gas pressure in pascals.
#define NEUTRAL_DENSITY (GAS_PRESSURE/((double)0.025*CHARGE))

/* VTE, ION_MASS, ION_CHARGE, VTI AND VTN  are direct inputs to the simulation. */
#define TEMPERATURE  (double)1.5
#define VTE sqrt(TEMPERATURE*CHARGE/ELECTRON_MASS)  // electron thermal velocity
#define ION_MASS ((double)AMU)               // mass of the ion
#define ION_CHARGE CHARGE                           // charge per ion
#define VTI sqrt(TEMPERATURE*CHARGE/ION_MASS)       // ion thermal velocity
#define VTN sqrt((double)0.025*CHARGE/ION_MASS)     //neutral thermal velocity

/* Used directly by the simulation for the ionization threshold */ 
//#define THRESHOLD (double)15.75  //Argon ionization threshold
#define THRESHOLD (double)13.6  //Hydrogen ionization threshold if instead choose hydrogen

/* Used directly by the simulation: the work function of the surface */
#define WORK_FUNCTION (double)4.6 

/* For either electrons or ions if the there are no collisions then select collision probability 
 * is constant and make the constant value zero.
 */

/* Used directly by the simulation: a flag to determine if there are collisions */
#define HAVE_COLLISIONS true  // true  // true if want collisions

/* Used directly by the simulation to decide which collision species are used. 
 * Options are ARGON or HYDROGEN 
 */
#define SPECIES  "HYDROGEN" // "ARGON" // "HYDROGEN"

// ----- for ions ----

/* Used directly by the simulation to decide if ion collision probability (collision cross-section) 
 * is constant or energy dependent for charge exchange collisions
 */
#define SIGMA "ENERGY_DEPENDENT"   //"ENERGY_DEPENDENT" //"CONSTANT"

/* Is only used to define a constant value of the collision cross-section for ions if 
 * select "CONSTANT" above
 */
#define SIGMA_VALUE 0.0   

// ---- for electrons ----

/* Used directly by the simulation to decide if electron elastic collision probability (collision cross-section)
 * is constant or energy dependent
 */
#define SIGMA_ELASTIC "ENERGY_DEPENDENT"  // "ENERGY_DEPENDENT" // "CONSTANT"

/* Is only used to define a constant value of the collision cross-sectionif select "CONSTANT" above */
#define ELASTIC_SIGMA_VALUE 0.0 //double type

/* Used directly by the simulation to decide if electron ionization collision probability 
 * (collision cross-section) is constant or energy dependent
 */
#define SIGMA_IONIZATION "ENERGY_DEPENDENT"  // "ENERGY_DEPENDENT" // "CONSTANT"

/* Is only used to define a constant value of the collision cross-section  if select "CONSTANT" above */
#define IONIZATION_SIGMA_VALUE 0.0 //double type


//---- select boundary options ------

/* Used directly by the simulation: select  YES to include electron impact secondary 
 * emmision from the surface
 */
#define E_IMPACT_SEE "YES"  // "YES" // "NO"

/* Used directly by the simulation: Only available for Argon Ions, select YES to include 
 * ion impact secondary electron emmision
 */
#define ION_IMPACT_SEE "NO"  // "YES" // "NO"

//select YES for reflection from surface
#define REFLECTION  "YES"  // "YES" // "NO"

/* Used directly by the simulation: select YES for backscatter from the surface */
#define BACKSCATTER  "YES"  // "YES" // "NO"

/* Used directly by the simulation:
 * Select YES for steady state boundary where an electron-ion pair are inserted into 
 * the central plasma when an ion reaches a surface. Note when an electron reaches
 * the surface it is absorbed and deleted. This boundary results in a rapid
 * steady state being achieved. 
 * WARNING: BOUNDARY EFFECTS, IONIZATION AND PARTICLE INJECTION MUST BE TURNED OFF WITH THIS OPTION
 * as the number of particles is treated as a constant.
*/
#define STEADYSTATE "NO" //  "NO" //"YES"

/* Used directly by the simulation:
 * The simulation uses the multigrid method (mglin function) to calculate the electrostatic potential 
 * as is described in the book Numerical Recipes (The Art of Scientific Computing)by Press et al.
 * For multigrid need to pass a square matrix where DX = DY and the dimension n is of the form 2^j + 1
 * where j is the number of grid levels used in the solution. These are kept independently in the simulation
 * to leave room for other potential solvers in the future.
 * WARNING: X_NO_OF_CELLS AND Y_NO_OF_CELLS NEED TO CONFORM TO ABOVE DESCRIPTION
 * The grid size is (X_NO_OF_CELLS x Y_NO_OF_CELLS)
*/
#define X_NO_OF_CELLS 257
#define Y_NO_OF_CELLS 257
#define NO_OF_CELLS X_NO_OF_CELLS 

/* Used directly by the simulation:
 * LENGTH, XLENGTH, YLENGTH are the input parameters to the program giving the length of the  
 * system in meters. Roughly 4-5 cells per debye length
 * WARNING XLENGTH, YLENGTH NEED TO BE THE SAME SIZE. These are kept independently in the simulation
 * to leave room for other potential solvers in the future.
 */
#define NO_OF_DEBYE_LENGTHS (int)64  
#define DEBYE_LENGTH sqrt((EPSILON*TEMPERATURE)/(ELECTRON_PARTICLE_DENSITY*CHARGE))
#define XLENGTH ((double)NO_OF_DEBYE_LENGTHS*DEBYE_LENGTH)
#define YLENGTH XLENGTH
#define LENGTH XLENGTH

/* Used directly by the simulation: DX, DY are the cell lengths in x and y. WARNING: THESE NEED
 * TO BE THE SAME SIZE. These are kept independently in the simulation to leave room for other
 * potential solvers in the future.
 */
#define DX (XLENGTH/(double)(X_NO_OF_CELLS-1))
#define DY (YLENGTH/(double)(Y_NO_OF_CELLS-1))

/* Used directly by the simulation: initial particle number loaded into the domain
 * easier if power of two -- should be a squared value
 */
#define INIT_PARTICLE_NO 512*512 // 1024*1024 // 

/* Used directly by the simulation: The macro particle size of each particle type.
 * See Plasma Physics Via Computer Simulation by Birdsall and Langdon for more
 * information.
 */
#define E_MACROPARTICLE_SIZE 1.0e6 
#define ION_MACROPARTICLE_SIZE 1.0e6 

/* Not used directly by the simulation. Is used in defining the timesteps */
#define PLASMA_FREQUENCY sqrt((ELECTRON_PARTICLE_DENSITY*ELECTRON_CHARGE2)/(EPSILON*ELECTRON_MASS))

/* Used directly by the simulation: DT, DTI, DT_DTI */
#define DT ((double)0.1/PLASMA_FREQUENCY)  //Electron time-step
//A different DTI should not be used when implementing a steady_state_boundary.
#define DTI ((double)1.0*DT) //Ion time-step
#define DT_DTI (int)(DTI/DT) //Ratio of ion time-step over electron

/* Used directly by the simulation: TWO_PI_F_DT */
#define POTENTIAL_FREQUENCY (double)13.6e6  // common frequency if use osciallating boundary potential
#define TWO_PI_F_DT (2.0*DT*PI*POTENTIAL_FREQUENCY)  
// use: PHI_RIGHT*cos(TWO_PI_F_DT*timestep) to calculate potential at right boundary
// see PHI_RIGHT below

/* Not used directly by simulation but instead to define RUN TIME, and NUMBER OF GRAPH STEPS below */
#define PERIOD (1.0/TWO_PI_F_DT)  

/* Used directly by the simulation only during batch mode and NOT DURING INTERACTIVE MODE: 
 * NO_GRAPH_STEPS, TOTAL_RUN_TIME, START_GRAPH_TIME 
 */
#define NO_GRAPH_STEPS ((int)50*(int)PERIOD)  //number of steps to average of data over
#define TOTAL_RUN_TIME ((int)300*(int)PERIOD + 5) //total run-time of simulation
#define START_GRAPH_TIME ((int)250*(int)PERIOD) //when to start collecting graph data points

/* Used directly by the simulation: when to stop collecting graph data points */
#define END_GRAPH_TIME (START_GRAPH_TIME+NO_GRAPH_STEPS)

/* Used directly by the simulation for taking data every TIME_DATA_FREQUENCY number of 
 * timesteps at the same location in space. 
 * SIZE_TIME_DATA is input parameter that defines the length of the arrays used for this. 
 */
#define TIME_DATA_FREQUENCY (int)100
#define SIZE_TIME_DATA TOTAL_RUN_TIME/TIME_DATA_FREQUENCY +(int)1

/* Used directly by the simulation: DXPE, DXPI, DYPE, DYPI.
 * Dxp is the distance between particles during the initial load at the start
 * of the simulation
 */
#define DXYPE XLENGTH/(sqrt(INIT_PARTICLE_NO) + 1.0)
#define DXYPI YLENGTH/(sqrt(INIT_PARTICLE_NO) + 1.0)
#define DXPE  DXYPE
#define DXPI  DXYPI
#define DYPE  DXYPE
#define DYPI  DXYPI

/* Used directly by the simulation: QC is the charge of each macroparticle */
#define QCI (ION_MACROPARTICLE_SIZE*ION_CHARGE)
#define QCE (E_MACROPARTICLE_SIZE*ELECTRON_CHARGE)

/* Used directly by the simulation: q_over_m is independent of the size of the macroparticle 
 * and is used in the equations of motion
 */
#define Q_OVER_ME (ELECTRON_CHARGE/ELECTRON_MASS)
#define Q_OVER_MI (ION_CHARGE/ION_MASS)
 
/* Used directly by the simulation: optional magnetic field on x-y plane, Bx = Bcos(angle),
 * By=Bsin(angle). Set magnetic field =0 for no field
 */
#define  MAGNETIC_FIELD (double)0.0
#define  ANGLE ((double)30.0)*PI/((double)180.0)

/* Used directly by the simulation: The right boundary can oscillate */
#define PHI_LEFT (double)0.0  // left boundary voltage
// right boundary voltage can also take the form: PHI_RIGHT*sin(PERIOD*(double)time_step)
#define PHI_RIGHT (double)0.0 
#define PHI_PERIOD (double)0.0   // voltage period
#define PHI_TOP (double)0.0  // voltage of top boundary
#define PHI_BOTTOM (double)0.0  // voltage of bottom boundary

/* Used directly by the simulation: NO_ION_SPECIES is an input parameter describing the number of species
 * The original 1D code accounted for multiple species. Much of the 2D code has left things
 * in place so this can potentially be added in the future, but currently this value must be set
 * to one. 
 */
#define NO_ION_SPECIES (int)1

/* Used directly by the simulation: Parameters used to make text files of particle velocity distributions
 * Number of velocity bins currently corresponds to the max number of rows in excel workbooks
 */
#define NO_BINS 30001  // number of velocity bins to partition the particle velocities into
#define DIVIDE_BINS_BY ((NO_BINS-1)/((int)100))  // used for scaling

/* Used directly by the simulation: locations on the grid to take the particle velocity distribution */
// location one:
#define DIST_LOCATIONX  (int)(0.9*(float)X_NO_OF_CELLS)
#define DIST_LOCATIONY  (int)(0.5*(float)Y_NO_OF_CELLS)
// location two:
#define DIST_LOCATION2X (int)(0.1*(float)X_NO_OF_CELLS)
#define DIST_LOCATION2Y (int)(0.5*(float)Y_NO_OF_CELLS)

/* Used directly by the simulation: The number of cells around the dist_location from which 
 * velocity data is pulled
 */
#define SPREAD (int)(0.05*(float)X_NO_OF_CELLS)

/* Used directly by the simulation: Which cell numbers to take time-variable data from */
// location one:
#define POINT1X  (int)3
#define POINT1Y  (int)(Y_NO_OF_CELLS/2)
// location two:
#define POINT2X  (int)(X_NO_OF_CELLS/2)
#define POINT2Y  (int)(Y_NO_OF_CELLS/2)
// location three:
#define POINT3X  (int)(X_NO_OF_CELLS -3)
#define POINT3Y  (int)(Y_NO_OF_CELLS/2)

/* Used directly by the simulation: Names used to distinguish text files of plot data */
#define NO_TYPE (char *)"n"  // no type
#define ION_TYPE (char *)"i1"  // ions
#define ELECTRON_TYPE (char *)"e"  // electrons

/* Used directly by the simulation:
 * parameters used for reducing the particle number when the particle vectors
 * get too large:
 * Limit the size of the vector particle arrays by reducing the particle number
 * according to the scheme set out by E E Kunhardt et al:
 * -Journal of Computational Physics,vol 67,(1986) 279
 * -----------------------------------------------------
 * This scheme involves dividing the particles up into velocity groups,
 * and decreasing the size of each group once the upper group size is reached.
 * Each remaining particle is then rewieghted according to the amount that the 
 * group is reduced by. The number of particle velocity groups and the upper 
 * size limit of each group is determined by the user. 
 */

/* Number of velocity groups to account for when reducing the particle number.  
 * Choose more than one velocity group to preserve high velocity tails in the
 * random reduction of the number of particles. This is a mandatory constant 
 * that is also used to construct group arrays.
 */
#define NO_GROUPS (int)3  

//Maximum allowable particle number of each species in each velocity group.
//Note: Different groups can have different sizes.
#define UPPER_GROUP_LIMIT (int)400000
#define UPPER_GROUP_LIMIT2 (int)100000

//Desired group size after reduction
#define NEW_GROUP_SIZE (int)200000
#define NEW_GROUP_SIZE2 (int)50000

//The largest velocity value to include in the reduction. Note: can exclude tail//tails.
#define UPPER_V_LIMIT ((double)1000.0*VTE)
#define UPPER_V_LIMIT_I ((double)1000.0*VTI)

#define VI1 (double)2.0*VTI
#define VI2 (double)3.0*VTI
#define VE1 (double)2.0*VTE
#define VE2 (double)3.0*VTE

/* Used directly by the simulation:
 * Putting particles into groups to check each group size is time consuming
 * so this is only done when the number of particles reaches the below limit.
 * Set to one to check each group size every time-step.
 * NOTE: THIS CHECK IS NOT DONE FOR A STEADY STATE BOUNDARY
 */
#define SIZE_TO_CHECK (int)1000001
#define SIZE_TO_CHECK_E (int)1000001

#endif
