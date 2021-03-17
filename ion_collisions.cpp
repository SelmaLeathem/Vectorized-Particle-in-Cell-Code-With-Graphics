/*
 * File: ion_collisions.cpp
 * -----------------------
 * 
 * Date: 1/25/2021
 * 
 * Remove scattering collision since question accuracy after some research.
 * 
 * Change the collision probablity calculation. This is an ongoing work in
 * progress. 
 * 
 * This file implements the ion_collisions.h interface. User defined collision
 * cross-section methods are declared in "ion_coll_input_declare" and user
 * defined cross-section methods are implemented in "ion_coll_input_methods".
 *
 * 6/15
 *
 */

#include "ion_collisions.h"

/*
* Constuctor: ion_collisions
* Usage: ion_collisions ion_coll;
*        ion_collisions ion_coll(IS_ION,dt);
*/
ion_collisions::ion_collisions(int if_electron_in, double dt_in, bool CM_elastic, bool CM_CX):collisions(if_electron_in,dt_in),center_mass_El(CM_elastic),center_mass_CX(CM_CX)
{
    probabilityCX = new double[NO_OF_CELLS]();
    if (!probabilityCX) {
        printf("Error allocating memory for probabilityCX\n");
        exit(1);
    }

    for (int i = 0; i < NO_OF_CELLS; i++)
    {
        probabilityCX[i] = 0.0;
    }
}

/* Destructor: ~ion_collisions */
ion_collisions::~ion_collisions()
{
    delete[] probabilityCX;
}

/*
* Method: set_neutral_density
* Usage:  set_neutral_density(NEUTRAL_DENSITY);
* ----------------------------------------------
* Set the neutral density. This value is used to calculate the collision
* probability.
*/
void ion_collisions::set_neutral_density(double density){
      neutral_density=density;
}

/*
* Method: set_vtn
* Usage:  set_vtn(neutral_thermal_velocity);
* -------------------------------------------
* Set the neutral thermal velocity.
*/
void ion_collisions::set_vtn(double neutral_velocity){
      vtn=neutral_velocity;
}

/*
* Method: set_ion_mass
* Usage:  set_ion_mass(ION_MASS);
* --------------------------------
* Set the ion_mass of the species the electrons are interacting with.
*/
void ion_collisions::set_ion_mass(double mass){
     ion_mass= mass;
}

/*
* Method: set_CMtrue_CX
* Usage:  set_CMtrue_CX(true);
* ----------------------------------
* Set center_mass_CX to true if the ion-atom center of mass velocity is
* required to calculate the charge-exchange cross-section. The default 
* is false.
*/ 
void ion_collisions::set_CMtrue_CX(bool val){
     center_mass_CX=val;
}

/*
* Method: set_crossSec_CX
* Usage: set_crossSec_CX(&ion_collisions::charge_exchange_crossSec_func);
* ------------------------------------------------------------------------
* Sets the 'double (ion_collisions::*sigma_CX)(void)' function pointer
* to point to the function passed in the arguement. User defined functions
* can be passed. This function returns energy-dependent charge-exchange
* cross sections.
*/
void ion_collisions::set_crossSec_CX(double (ion_collisions::*func1)(void)){
     sigma_CX= func1;
}

/*
* Method: set_constant_CX
* Usage:  set_constant_CX(constant_cross_section_val);
* -----------------------------------------------------------
* Sets the constant charge-exchange cross-section value, if it is constant.
* This variable does not need to be set. If the cross-section is not constant
* do not set this value, otherwise, this must be set -in addition- to setting
* the cross-section function pointer to point to the
* sigma_const_cx method defined in this file.
*/
void ion_collisions::set_constant_CX(double val){
     constant_cx =val;
}

/*
* Method: calc_CM_energy
* Usage: calc_CM_energy(x_velocity,y_velocity,z_velocity)
* ------------------------------------------------------
* Calculates the particle center-of-mass energy for use in cross-section
* calculations.
*/
void ion_collisions::calc_CM_energy(double vx,double vy,double vz){
     double ux,uy,uz;

     ux=(ion_mass*vx+ion_mass*vtn);
     uy=(ion_mass*vy+ion_mass*vtn);
     uz=(ion_mass*vz+ion_mass*vtn);

     CMenergy=(0.5/CHARGE)*(1.0/(2.0*ion_mass))*(ux*ux +uy*uy +uz*uz);
}

/*
* Method: calc_energy
* Usage: calc_energy(x_velocity,y_velocity,z_velocity)
* ------------------------------------------------------
* Calculates the particle energy for use in cross-section and
* collision calculations.
*/
void ion_collisions::calc_energy(double vx,double vy,double vz){
     double u2;

     u2 = vx*vx + vy*vy +vz*vz;
     energy = (0.5/CHARGE)*ion_mass*u2;
}


/*
 * Implementation Notes: sigma_const_cx
 * --------------------------------------
 * For a constant value or no collision need to set constant_cx to a value
 * or zero before passing this method to the set_crossSec_CX method.
 *
 * e.g.
 * object.set_constant_CX(value);
 * object.set_crossSec_CX(&ion_collisions::sigma_const_cx);
 */
double ion_collisions::sigma_const_cx(void){
      return constant_cx;
}


/*
 * Implementation Notes: sigma_Ar_cx
 * ----------------------------------------
 * Calculates the ion-atom energy dependent cross-section with formulae
 * pinched from Berkeley's pdp1 code.
 */
double ion_collisions::sigma_Ar_cx(void){
   
   if(energy > 4.0)
      return(2.0 + 5.5/sqrt(energy))*10.0;

   return(-2.95*sqrt(energy) + 10.65)*10.0;
}

/*
 * Implementation Notes: sigma_H_cx
 * ---------------------------------------
 * Calculates the ion-atom energy dependent cross-section. The formulae are from
 * -R K Janev et al, Atomic and Plasma-Material Interaction Data for Fusion,
 * vol 4,(1993)
 */
double ion_collisions::sigma_H_cx(void) {
    double sigma;
    double A1 = 3.2345, A2 = 235.88, A3 = 0.038371, A4 = 3.8068e-6, A5 = 1.1832e-10, A6 = 2.3713;
    energy = energy / 1000.0;
    energy = energy / AMU;

    if (energy < 1.0e-4)
        sigma = 0.0;

    else
        sigma = (A1 * log((A2 / energy) + A6)) / (1.0 + A3 * energy + A4 * pow(energy, 3.5)
            + A5 * pow(energy, 5.4));

    return sigma; 
}


/*
 * Implementation Notes: charge_exchange_velocity
 * ------------------------------------------------
 * Implement charge_exhange collision by giving ions a post-collision velocity
 * equall to the neutral incident speed as described in:
 * -V Vahedi et al,Computer Physics Communications,vol 87,(1995)179
 */
void ion_collisions::charge_exchange_velocity(int i, vector<particleVec> &particleObj)
{
   int index;
   double R;

   R = (double)rand()/((double)RAND_MAX);
   index= rand()%(INIT_NO_NEUTRALS - 1);

   particleObj[i].particleVec::vx = neutral_vx[index];
   particleObj[i].particleVec::vy = neutral_vy[index];
   particleObj[i].particleVec::vz = neutral_vz[index];

}

/*
 * Implementation Notes: collision
 * ---------------------------------
 * Determines if a collision occurs, and calls the collision methods if it
 * does. This is done in a few steps.
 *
 *1)Determines if a collision occurs by comparing a random number to the
 *  total collision probability as is described in the following articles:
 *  -V. Vahedi et al,Computer Physics Communications,vol 87,(1995) 179
 *  -C K Birdsall, IEEE Trans. Plasma Sci.,vol 19,(1991) 65
 *  -Kenichi Nanbu et al,IEEE Trans. Plasma Sci.,vol 28,(2000) 971
 * 
 *   Use Birdall, IEEE article's method for the probability
 *
 *2)If a collision occurs the type of collision is determined by comparing
 *  each individual cross-section to random numbers as described in the
 *  articles listed in step one.
 *
 *3)If a collision occurs the relevant collision method to implement the
 *  the collision is called with the coll_factor set to -reduced mass
 *  
 *  See Kenichi Nanbu and Takizuka for more information about the coll_factor.
 *  -Tomonori Takizuka et al,Journal of Computational Physics,vol 25,(1977)205
 *
 * Need to verify that have implemented tests to find out if collision occurs
 * correctly.
 *
 * Currently implement an elastic collision with an atom by scattering the ion
 * using the scattering formulae in Nanbu and Takizuka which serves as a 
 * temporary place holder until this is properly looked into.
 *
 */
void ion_collisions::collision(vector<particleVec> &particleObj)
{
   int count_cx=0; // count number of collisions
   int i=0, particle_no, jp, kp;
   double sigma_iEl, sigma_iCx, sigma_total; // cross sections
   double R;
   double u;
   double reduced_mass,coll_factor;
   double Pmax=1.0; 

   particle_no = particleObj.size();

  /*treat ion mass as same wieght as atom */ 
   reduced_mass = 1.0/2.0;

   while(i < particle_no)
   {
       // identify location on the grid
       jp = (int)(particleObj[i].x / DX);
       kp = (int)(particleObj[i].y / DY);

    // if within the boundaries then is eligible for a consideration
    if (jp >= 0 && jp < X_NO_OF_CELLS && kp >= 0 && kp < Y_NO_OF_CELLS) {
             
        // calculate the particle energy 
      calc_energy(particleObj[i].vx,particleObj[i].vy,particleObj[i].vz);

      u= sqrt( (2.0*CHARGE*this->energy)/ion_mass);

        //if need center of mass energy then calculate it
      if( center_mass_CX==true)
         calc_CM_energy(particleObj[i].vx,particleObj[i].vy,particleObj[i].vz);

      sigma_iCx = (this->*sigma_CX)();

     /* Charge exchange collisions */

    // from Birdsall IEEE Trans Sci articles page 75
    // neutral density*sigma = 1/mean-free-path
    // sum the probability for each particle
    // when the probability > 1 a collision occurs and the probability is reset to 0

      probabilityCX[jp] += (1.0 -exp(-(u*dt* 1.0e-20 * neutral_density *sigma_iCx)));

      if(probabilityCX[jp] > 1.0 ){
               charge_exchange_velocity(i,particleObj);
               count_cx += 1;  //used to monitor the number of collisions
               probabilityCX[jp] = 0.0;
      } /*goes with if R > */     
    } /* if x within boundaries */
      i+=1;
   } /* goes with while loop */

}
