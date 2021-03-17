/*
 * File: collisions.h
 * ------------------- 
 * 
 * This file defines the collisions class which serves as a base class for
 * the e_collisions and ion_collisions clases.
 *
 * 5/15
 * 
 * Date: 2/28/2021
 * 
 * The scattering function defined in this class currently only applies to electrons
 * although is expandable to include ions. This class is kept as a parent to 
 * e_collisions and ion_collisions since at some future date it might provide abstract
 * function headers should ion scattering additionaly be accounted for.
 */

#ifndef collisions_h
#define collisions_h

#define _CRT_SECURE_NO_WARNINGS

#include "constants.h"
#include "parameters.h"
#include "particleVec.h"
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>

using namespace std;

class collisions
{

 public:

  /*
   * Constructor: collisions
   * Usage: collisions coll_object;
   *        collisions coll_object(is_an_electron,dt);  
   */ 
   collisions(int if_electron_in, double dt_in);

  /* Destructor: ~collisions */
   ~collisions();

  /*
   * Method: calc_neutral_velocity
   * Usage:  calc_neutral_velocity();
   * ---------------------------------
   * Calculates the neutral velocities in the neutral velocity array for
   * that species before the simulation starts. Neutral velocities are randomly
   * chosen from the arrays, when collisions are implemented instead of 
   * making the calculations at the time they are needed. 
   */ 
   void calc_neutral_velocity(void);

  /*
   * Method: scattered_velocity
   * Usage: scatterd_velocity(element_number,particle_vector,particle_energy,
   *                                                    coll_factor,YES)
   * --------------------------------------------------------------------------
   * Calculates the post-collision velocities after a scattering event. Note:
   * the particle_energy is in eV.
   */
   void scattered_velocity(int R1, vector<particleVec> &particleObj, double energy, double coll_factor,int yes_or_no);

   protected:
   
   /*Determines if the particle is an electron or ion. This is used
    *to calculate post-collision scattering angles. 
    */
   int if_electron;  
   
   double dt;  // size of the timestep

   /*Neutral velocities are randomly chosed from predefined arrays */
   double *neutral_vx;
   double *neutral_vy;
   double *neutral_vz;

   double neutral_density; 
   double vtn; // neutral velocity
   double ion_mass; // ion mass
};

#endif /*collisions_h*/
