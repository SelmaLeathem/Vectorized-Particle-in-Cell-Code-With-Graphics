/*
 * File: boundary.h
 * -------------------
 * 
 * 
 * Definition the boundary class, which stores variables and 
 * functions that implement the particle boundary. Electrons are either
 * reflected, backscattered, absorbed, or used to generate secondary   
 * electrons at the walls. All ions are absorbed, and some generate 
 * secondary electrons.
 *
 * Secondary yield functions are passed to a function pointer enabling
 * user defined functions to be used. 
 *
 * For a completely absorbing boundary, or to turn off a type of process, set
 * the corresponding yield(s) to zero.
 * 
 * In addition there is a special boundary used to implement a rapid steady-
 * state by placing an electron-ion pair in the plasma when an ion hits the
 * boundary.
 *
 * 10/15
 * 
 * Date: 2/28/21
 * Description: Adjust functions to make calculations on a two-dimensional
 *              grid instead of a one-dimensional line. The original code
 *              above only accounted for two boundaries at x = 0 and x = length.
 *              Now there are two additional boundaries at y = 0 and y = length.
 *              Both boundaries undergo the same events.
 * 
 */

#ifndef boundary_h
#define boundary_h

#define _CRT_SECURE_NO_WARNINGS   //Visual Studio Uses to Ignore VS related warnings

#include "constants.h"
#include "parameters.h"
#include "particleVec.h"
#include <vector>
#include <cmath>

using namespace std;

/* The particle boundary is implemented through the ion_boundary and e_boundary
 * functions. 
 *
 */

    
class boundary
{
  public:

 
     /*
      * Constructor: boundary
      * Usage: boundary e_boundary;
      *        boundary e_boundary(time_step_size);
      */
     boundary(double dt=DT);

     /* Destructor: ~boundary */
     ~boundary();

     /* 
      * Method: calc_energy
      * Usage:  calc_energy(x_velocity,y_velocity,z_velocity);
      * ------------------------------------------------------
      * Calculates the particle energy in eV for use in yield
      * calculations. Note, that user defined yield functions can access
      * this value.
      */ 
     void calc_energy(double vx,double vy,double vz);

     /* 
      * Method: calc_angle
      * Usage:  calc_angle(x_velocity);
      * -------------------------------
      * Calculates the electron impact angle for use in yield
      * calculations. Note, that user defined yield functions can access
      * this value.
      */ 
     void calc_angleX(double vx);
     void calc_angleY(double vy);

     /*
      * Method: set_e_impact_yield
      * Usage:  set_e_impact_yield(&boundary::see_yield_function);
      * ----------------------------------------------------------
      * Sets the electron e_impact_see_yield function pointer to point to a
      * function determined by the user, including a user defined one.
      * This function calculates the electron-impact yield.
      *
      */
     void set_e_impact_yield(void (boundary::*func1)(void));

     /*
      * Method: set_i_impact_yield
      * Usage:  set_i_impact_yield(&boundary::see_yield_function);
      * ----------------------------------------------------------
      * Sets the ion i_impact_see_yield function pointer to point to a
      * function determined by the user, including a user defined one.
      *
      */
     void set_i_impact_yield(void (boundary::*func1)(void));

     /*
      * Method: set_reflection_yield
      * Usage:  set_reflection_yield(&boundary::reflection_yield_function);
      * ----------------------------------------------------------
      * Sets the electron set_reflection_yield function pointer to point to a
      * function determined by the user, including a user defined one.
      *
      */
     void set_reflection_yield(void (boundary::*func1)(void));


     /*
      * Method: set_backscatter_yield
      * Usage:  set_backscatter_yield(&boundary::backscatter_yield_function);
      * ----------------------------------------------------------
      * Sets the electron set_backscatter_yield function pointer to point to a
      * function determined by the user, including a user defined one.
      *
      */
     void set_backscatter_yield(void (boundary::*func1)(void));

     /*
      * Sample yield functions 
      *
      * void e_impact_see_yield_v(void);
      * void i_impact_see_yield_Ar(void);
      * void calc_reflection_yield_v(void);
      * void calc_backscatter_yield_v(void);
      *
      */
     

     /*
      * Method: reflect_electron
      * Usage:  reflect_electron(particle_number,particle_vector)
      * --------------------------------------------------------------
      * Reflects the particle_vector element identified by the particle_number.
      *
      */
      
     void reflect_electron(int current_particle,vector<particleVec> &particleObj);
     
     /*
      * Method: secondary_vel
      * Usage:  secondary_vel(current_particle_number,new_velocity,
      *                                                    particle_vector)
      * -------------------------------------------------------------------
      * Calculates the new velocity of a backscattered or true secondary
      * electron due to electron impact.
      *
      */       
  
      void secondary_vel(int current_particle,double velocity,vector<particleVec> &particleObj);
     

     /*
      * Method: i_secondary_vel
      * Usage:  i_secondary_vel(last__e_particle_number,electron_vector,
      *                             ion_x_position,ion_y_position,ion_weight)
      * -------------------------------------------------------------------
      * Calculates the velocity of secondary electrons due to ion impact.
      * Before this function is called a new secondary electron is added to the 
      * end of the electron vector.
      *
      */
      void i_secondary_vel(int last_particle, vector<particleVec> &electrons, double x, double y, double weight);
     
     /*
      * Method: secondary_electrons
      * Usage:  secondary_electrons(current_particle_number,electron_vector)
      * --------------------------------------------------------------------
      * Calculates the position, wieght and velocity (via the secondary_vel  
      * function) of secondary electrons.
      *
      */
     void secondary_electrons(int current_particle,vector<particleVec> &particleObj);

     /*
      * Method: plop_particles
      * Usage:  plop_particles(elecElement,ionElement,vte,vti);
      * ------------------------------------------------
      * Plops an ion-electron pair into the simulation domain.
      *
      */
      void plop_particles(particleVec &elecElement,particleVec &ionElement,
                            double vte, double vti);

     /*
      * Method: get_particle_velocity
      * Usage:  double v = get_particle_velocity(thermal_velocity);
      * --------------------------------------------
      * Returns a random particle velocity taken off a Maxwellian distribution.
      * This method is a temporary place holder for a more accurate one.
      */
     double get_particle_velocity(double thermal_velocity);
    
     /*
      * Method: ion_boundary
      * Usage:  ion_boundary(e_particle_vector,ion_particle_vector)
      * --------------------------------------------------------------
      * Implements the ion particle boundary. All ions are absorbed at the
      * surface, some create secondary electrons prior to absorption.
      */
      void ion_boundary(vector<particleVec> &electrons,vector<particleVec> &ions); 
     
     /*
      * Method: e_boundary
      * Usage:  e_boundary(e_particle_vector,ion_particle_vector)
      * ------------------------------------------------------
      * Implements the electron particle boundary. Electrons are either
      * reflected, backscattered, absorbed, or create secondary electrons.
      *
      */
      void e_boundary(vector<particleVec> &electronObj,vector<particleVec> &ionObj);
  
     /*
      * Method: steady_e_boundary
      * Usage:  steady_e_boundary(electron_vector,ion_vector);
      * ---------------------------------------------------------
      * Implements a type of boundary which results in a rapid steady state.
      * When an ion arrives at a boundary it sits there. When an electron
      * arrives it and one of the boundary ions gets put back into the 
      * simulation domain. Ionization is switched off when this option is
      * selected. Boundary::ion_proportion needs to be set prior to using this.
      */
      void steady_e_boundary(vector<particleVec> &electronObj,vector<particleVec> &ionObj);
      // For boundaries at x = 0 and x = Xlength
      void steadyBoundaryX(vector<particleVec> &electronObj,vector<particleVec> &ionObj);
      // For boundaries at y = 0 and y = Ylength
      void steadyBoundaryY(vector<particleVec> &electronObj,vector<particleVec> &ionObj);

     /*
      * Method: steady_i_boundary
      * Usage:  steady_i_boundary(electron_vector,ion_vector);
      * ---------------------------------------------------------
      * Implements a type of boundary which results in a rapid steady state.
      * When an ion arrives at a boundary it sits there. When an electron
      * arrives it and one of the boundary ions gets put back into the 
      * simulation domain. Ionization is switched off when this option is
      * selected.
      */
      void steady_i_boundary(vector<particleVec> &electronObj,vector<particleVec> &ionObj);
 
     /*
      * Method: absorb
      * Usage:  absorb(particle_vector)
      * ---------------------------------
      * Implements an absorbing boundary, which is quicker than setting the
      * yields to zero. Currently only used for debugging.
      * 
      */ 
      void absorb(vector<particleVec> &particleObj);

     /*
      * Set parameters commonly used by functions that calculate yields.
      * Note that user defined functions are able to use set and access 
      * these values.
      *
      * void set_max_yield(double val); 
      * void set_E0(double val);     
      * void set_Emax0(double val); 
      * void set_threshold(double val);
      * void set_work_function(double val);
      * void set_incident_mass(double val);
      * void set_vt(double val);
      * void set_ion_proportion(double val);
      *
      */

     /*
      * Method: set_max_yield
      * Usage:  set_max_yield(max_yield_value)
      * ---------------------------------------
      * Sets the maximum secondary emmision coefficient at normal incidence.
      */
      void set_max_yield(double val); 
     
     /*
      * Method: set_Emax0
      * Usage:  set_Emax0(Emax0_value)
      * -------------------------------
      * Sets the energy corresponding to the maximum secondary emmision
      * coefficient at normal incidence max_yield.
      */
     void set_Emax0(double val); 

     /*
      * Method: set_E0
      * Usage:  set_E0(E0_value)
      * ------------------------
      * Sets the secondary emmision threshold in eV.
      */  
      void set_E0(double val);     

     /*
      * Method: set_threshold
      * Usage:  set_threshold(threshold_value)
      * ---------------------------------------
      * Sets the ionization threshold of the atom species in eV. 
      */
      void set_threshold(double val);

     /*
      * Method: set_work_function
      * Usage:  set_work_function(work_function_value)
      * -----------------------------------------------
      * Sets the surface work function in eV.
      */
      void set_work_function(double val);

     /*
      * Method: set_incident_mass
      * Usage:  set_incident_mass(incident_mass_value)
      * -----------------------------------------------
      * Sets the mass of the particle impinging upon the wall in kg.
      */
      void set_incident_mass(double val);

     /*
      * Method: set_vt
      * Usage:  set_vt(thermal_velocity_value)
      * -----------------------------------------------
      * Sets the thermal velocity.
      */
      void set_vt(double val);

     /*
      * Method: set_ion_proportion
      * Usage:  set_ion_proportion(50.0);
      * -----------------------------------
      * Set the loading proportion of each ion species. If only one species
      * then this value is 100. As an example an initial plasma of 60% species
      * A and 40% of species B would use (60) and (40) as arguements to 
      * separate calls.
      * This is used in the steady state boundary.
      */
      void set_ion_proportion(double val);

     /*
      * Method: set_implement_e_boundary
      * Usage:  set_implement_e_boundary(&boundary::boundary_function);
      * --------------------------------------------------------------
      * Sets the 'void (boundary::*implement_e_boundary)(...)' function pointer
      * to point to the function passed in the arguement.                  
      */  
      void set_implement_e_boundary(void (boundary::*func1)(vector<particleVec> &electronObj,vector<particleVec> &ionObj) );

     /*
      * Method: set_implement_i_boundary
      * Usage:  set_implement_i_boundary(&boundary::boundary_function);
      * --------------------------------------------------------------
      * Sets the 'void (boundary::*implement_i_boundary)(...)' function pointer
      * to point to the function passed in the arguement.                  
      */  
      void set_implement_i_boundary(void (boundary::*func1)(vector<particleVec> &electronObj,vector<particleVec> &ionObj) );

 /*
  * Method: (boundary::*implement_e_boundary)(vector<particleVec> &electrons, *vector<particleVec> &ions)
  * Usage:  (*implement_e_boundary)(electron_vector,ion_vector)
  * ------------------------------------------------
  * Implements the selected boundary type. Options include absorbing, boundaries  * with secondary electrons, and a rapid steady state boundary where ionization  * is switched off and an electron-ion pair is inserted in the simulation
  * domain everytime an electron arrives at the surface.
  */
   void (boundary::*implement_e_boundary)(vector<particleVec> &electronObj, 
                                            vector<particleVec> &ionObj);   

 /*
  * Method: (boundary::*implement_i_boundary)(vector<particleVec> &electrons,     *                                              vector<particleVec> &ions)
  * Usage:  (*implement_i_boundary)(electron_vector,ion_vector)
  * ------------------------------------------------
  * Implements the selected boundary type. Options include absorbing, boundaries  * with secondary electrons, and a rapid steady state boundary where ionization  * is switched off and an electron-ion pair is inserted in the simulation
  * domain everytime an electron arrives at the surface.
  */
   void (boundary::*implement_i_boundary)(vector<particleVec> &electronObj, 
                                            vector<particleVec> &ionObj);   

   /*
 * Implementation Notes: e_impact_see_yield_v
 * -------------------------------------------
 * Implement Vaughan model:
 *   -M Radmilovic-Radjenovic et al , Journal of Physics: Conference Series 71
 *    (2007) 012007
 *   -J Rodney M Vaughan, IEEE Transaction on Electron Devices,Vol 36,(1989)
 *    1963
 *
 */
   void e_impact_see_yield_v(void);

/*
 * Implementation Notes: i_impact_see_yield_Ar
 * ---------------------------------------------
 * Ion impact secondary electron yield for Argon ions impinging a dirty surface.
 * See Also:
 *   -M Radmilovic-Radjenovic et al, Eur. Phys. J. D,Vol 54,(2009),445
 *
 */
   void i_impact_see_yield_Ar(void);

  /* Prevent excess at high energies: reflection and backscatter references:
  * Ralf Krinke and H M Urbassek,  J. Phys. D: Appl. Phys., vol 29 (1996) 378
  * M Radmilovic-Radjenovic and Z L Petrovic, Eur. Phys. J. D, vol 54 (2009) 445
  */
   void calc_reflection_yield_v(void);

  /* Prevent excess at high energies: reflection and backscatter references:
  * Ralf Krinke and H M Urbassek,  J. Phys. D: Appl. Phys., vol 29 (1996) 378
  * M Radmilovic-Radjenovic and Z L Petrovic, Eur. Phys. J. D, vol 54 (2009) 445
  */
   void calc_backscatter_yield_v(void);


   /* Set yield functions to return zero in constructor. The user and program then
    * only sets the yields they want to use.
    */
   void set_backscatter_yield_to_zero(void);
   void set_reflection_yield_to_zero(void);
   void set_i_impact_see_yield_Ar_to_zero(void);
   void set_e_impact_see_yield_to_zero(void);


  private:

    //used for multispecies simulations in steady state to determine what percent
    //of the ions get paired up with an electron to be plopped back in. Typically
    // its fifty-fifty for two ion species
   double ion_proportion; 

   // mass of particle incident on the surface
   double incident_mass; 

   // incoming particle enery and angle to the surface
   double energy, angleX, angleY;

   // the maximum secondary yield corresponding to the maximum incoming
   // energy 
   double max_yield;  
   double Emax0; 

   double E0; //secondary emmision threshold

   // secondary yield
   double yield;
   double reflection_yield, backscatter_yield;

   // minimum threshold and work function of the surface 
   double threshold, work_function;

   double dt; //time-step
   double vt; //thermal velocity
   double wall_offset; // used to approximate == for a double

   bool onXboundary; // true if examining an x-directed boundary

/** The following function pointers point to functions selected by the user.
 **/

 /*
  * Method: (boundary::*e_impact_see_yield)(void)
  * Usage:  (this->*e_impact_see_yield)()
  * ----------------------------------------------
  * Calculates the electron impact secondary yield.
  */
   void (boundary::*e_impact_see_yield)(void);

 /*
  * Method: (boundary::*i_impact_see_yield)(void)
  * Usage:  (this->*i_impact_see_yield)()
  * ----------------------------------------------
  * Calculates the ion impact secondary yield.
  */
   void (boundary::*i_impact_see_yield)(void);

 /*
  * Method: (boundary::*calc_reflection_yield)(void)
  * Usage:  (this->*calc_reflection_yield)()
  * ----------------------------------------------
  * Calculates the electron reflection yield.
  */
   void (boundary::*calc_reflection_yield)(void);

 /*
  * Method: (boundary::*calc_backscatter_yield)(void)
  * Usage:  (this->*calc_backscatter_yield)()
  * ----------------------------------------------
  * Calculates the electron backscatter yield.
  */
   void (boundary::*calc_backscatter_yield)(void);


};


#endif
