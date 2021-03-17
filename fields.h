/*
 * File: fields.h
 * ----------------
 * 
 * 
 * Field class definition.  Calculates field values on the grid including
 * the charge density (rho), the electrostatic potential (phi) and the
 * electric field (E).
 *
 * 4/15
 *
 * Date: 2/6/2021
 * 
 * Alter the fields class from one-dimensional calculations to two-dimensional
 * calculations on a square grid. The electrostatic potential is calculated 
 * using the multigrid method using a multigrid solver from:
 * Numerical Recipes: The Art of Scientific Computing by William H. Press et al.
 * Since the electric field is a vector quantity E has become Ex and Ey for the
 * x and y components of the electric field on a two-dimensional grid.
 * 
 *        Multigrid solves Poission equation for the potential on a two-dimensional grid 
 *        in O(n) as compared to O(n^3) for straight LU decomposition.
 */

#ifndef fields_h
#define fields_h

#include "constants.h"
#include "parameters.h"
#include <cmath>
#include "particle.h"
#include "mglin.h"

#ifdef WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif
#include <omp.h>

#include "CL/cl.h"
#include "CL/cl_platform.h"

#include "openCLutility.h"

using namespace std;

class fields
{
   public:
    
    /*
     * Constructor: fields
     * Usage: fields field_object;
     *
     */
     fields();

    /* Destructor: ~fields */
     ~fields();

    /*
     * Method: grid_rho
     * Usage:  grid_rho(ion_objects,e_grid_charge);
     * ----------------------------------------------
     * Calculates the net charge density on each grid point.
     *
     */
     void grid_rho(particle *ion, cl_double *qe);

    /*
     * Method: grid_potential
     * Usage:  grid_potential(lhs_boundary_potential,rhs_boundary_potential, timestep);
     * ------------------------------------------------------------------------
     * Calculates the the electrostatic potential at each grid point due to 
     * the net charge density as defined by Poisson's equation.
     *
     */
     void grid_potential(double phi0X, double phiLX, double phi0Y, double phiLY , double timestep, double phi_period);

    /*
     * Method: grid_electric_field
     * Usage:  grid_electric_field();
     * -------------------------------     
     * Calculates the electric field on the grid from the electrostatic  
     * potential.
     */
     void grid_electric_field(void);

    /* Return pointers to private field variables. These are called from 
     * within the graph  and particle classes.
     *
     *     double* get_elec_field(void){return electric_field;}
     *     double* get_rho(void){return rho;}
     *     double* get_phi(void){return phi;}
     *     double  get_phi(int elem){return phi[elem];}
     */

    /*
     * Method: get_elec_field
     * Usage: double *E_field = object.fields::get_elec_field();
     * ---------------------------------------------------------
     * Returns a pointer to the electric field array. This is called
     * from the particle and graph class methods.
     */ 
     double* get_elec_fieldX(void){return electric_fieldX;}

     /*
     * Method: get_elec_field
     * Usage: double *E_field = object.fields::get_elec_field();
     * ---------------------------------------------------------
     * Returns a pointer to the electric field array. This is called
     * from the particle and graph class methods.
     */ 
     double* get_elec_fieldY(void){return electric_fieldY;}
 
    /*
     * Method: get_rho
     * Usage: double *rho= object.fields::get_rho();
     * ----------------------------------------------
     * Returns a pointer to the net charge density array. This is called
     * by methods in the graph class.
     */
     double** get_rho(void){return rho;}

    /*
     * Method: get_phi
     * Usage: double *phi= object.fields::get_phi();
     * -----------------------------------------------
     * Returns a pointer to the electrostatic potential array. This is 
     * called by methods in the graph class.
     */
     //double get_phi(void){return phi;}
     MatDoub_IO get_phi(void){return phi;}  

    /*
     * Method: get_phi
     * Usage: double phi_element= object.fields::get_phi(element_number);
     * -------------------------------------------------------------------
     * Returns the electrostatic potential at the specified grid point. This
     * is used to graph the variation of the potential at a specific point in
     * time.
     */
     double  get_phi(int j, int i){return phi[j][i];}
   

     private:

     /* Instance variables */

     double **rho; // the charge density
     MatDoub_IO phi;  //the electrostatic potential
     cl_double* electric_fieldX; // the electric field in the x direction
     cl_double* electric_fieldY; // the electric field in the y direction

};

#endif
