/*
 * File: e_collisions.cpp
 * -----------------------
 * 
 * This file implements the e_collisions class which stores data and defines
 * methods that calculate a variety of electron impact collisions including
 * scattering and ionization.
 *
 * 6/15
 * 
 * Date: 1/25/2021
 * Change ionization collisions to use:
 * 
 * -Kenichi Nanbu et al,IEEE Trans. Plasma Sci.,vol 28,(2000) 971
 *  Collision probabilities are a work in progress.
 */

#include "e_collisions.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

/*
* Constuctor: e_collisions
* Usage: e_collisions e_coll;
*        e_collisions e_coll(IS_ELECTRON,dt);
*/
e_collisions::e_collisions(int if_electron_in,double dt_in):collisions(if_electron_in,dt_in)
{
    countQualifying = 0;

    probabilityE = new double[NO_OF_CELLS]();
    if (!probabilityE) {
        printf("Error allocating memory for probabilityE\n");
        exit(1);
    }

    probabilityI = new double[NO_OF_CELLS]();
    if (!probabilityI) {
        printf("Error allocating memory for probabilityI\n");
        exit(1);
    }

    for (int i = 0; i < NO_OF_CELLS; i++)
    {
        probabilityE[i] = 0.0;
        probabilityI[i] = 0.0;
    }
}

/* Destructor: ~e_collisions */
e_collisions::~e_collisions()
{
    delete[] probabilityE;
    delete[] probabilityI;
}

/*
* Method: set_neutral_density
* Usage:  set_neutral_density(NEUTRAL_DENSITY);
* ----------------------------------------------
* Set the neutral density. This value is used to calculate the collision
* probability.
*/
void e_collisions::set_neutral_density(double density){
      neutral_density=density;
}

/*
* Method: set_vtn
* Usage:  set_vtn(neutral_thermal_velocity);
* -------------------------------------------
* Set the neutral thermal velocity.
*/
void e_collisions::set_vtn(double neutral_velocity){
      vtn=neutral_velocity;
}

/*
* Method: set_ion_mass
* Usage:  set_ion_mass(ION_MASS);
* --------------------------------
* Set the ion_mass of the species the electrons are interacting with.
*/
void e_collisions::set_ion_mass(double mass){
      ion_mass= mass;
}

/*
* Method: set_constant_elastic
* Usage:  set_constant_elastic(constant_cross_section_val);
* -----------------------------------------------------------  
* Sets the constant elastic cross-section value, if it is constant. This
* variable does not need to be set. If the cross-section is not constant
* do not set this value, otherwise, this must be set -in addition- to
* setting the cross-section function pointer to point to the
* sigma_const_elastic method defined in this file.
* 
*/
void e_collisions::set_constant_elastic(double val){
      constant_elastic = val;
}

/*
* Method: set_constant_ionization
* Usage:  set_constant_ionization(constant_cross_section_val);
* -----------------------------------------------------------  
* Sets the constant ionization cross-section value, if it is constant. This
* variable does not need to be set. If the cross-section is not constant
* do not set this value, otherwise, this must be set -in addition- to setting
* the cross-section function pointer to point to the 
* sigma_const_ionization method defined in this file.
*/
void e_collisions::set_constant_ionization(double val){
      constant_ionization = val;
}

 /*
* Method: set_crossSec_elastic
* Usage:  set_crossSec_elastic(&e_collisions::elastic_crossSec_function);
* ------------------------------------------------------------------------
* Sets the 'double (e_collisions::*sigma_elastic)(void)' function pointer to 
* point to the function passed in the arguement. User defined functions
* can be passed. This function returns energy-dependent elastic
* cross sections.
*/
void e_collisions::set_crossSec_elastic(double (e_collisions::*func1)(void)){
      sigma_elastic = func1;
}

/*
* Method: set_crossSec_ionization
* Usage: set_crossSec_ionization(&e_collisions::ionization_crossSec_func);
* ------------------------------------------------------------------------
* Sets the 'double (e_collisions::*sigma_ionization)(void)' function pointer 
* to point to the function passed in the arguement. User defined functions
* can be passed. This function returns energy-dependent ionization
* cross sections.
*/
void e_collisions::set_crossSec_ionization(double (e_collisions::*func1)(void)){
      sigma_ionization= func1;
}

/*
* Method: set_ionization_threshold
* Usage:  set_ionization_threshold(threshold_value);
* ----------------------------------------------------  
* Sets the ionization threshold in eV.
*/
void e_collisions::set_ionization_threshold(double val){
      threshold = val;
}

/*
* Method: calc_energy
* Usage: calc_energy(x_velocity,y_velocity,z_velocity)
* ------------------------------------------------------
* Calculates the particle energy for use in cross-section and
* collision calculations.
*/
void e_collisions::calc_energy(double vx,double vy,double vz){
     double u2;
     u2 = vx*vx + vy*vy +vz*vz;
     energy = (0.5/CHARGE)*ELECTRON_MASS*u2;
}

/*  
 * Implementation Notes: sigma_const_elastic
 * -------------------------------------------- 
 * For a constant value or no collision need to set constant_elastic to a value
 * or zero before passing this method to the set_crossSec_elastic method.
 *
 * e.g. 
 * object.set_constant_elastic(0.0);
 * object.set_crossSec_elastic(&e_collisions::sigma_const_elastic);
 */
double e_collisions::sigma_const_elastic(void){
        return constant_elastic;
}

/*  
 * Implementation Notes: sigma_const_ionization
 * --------------------------------------------- 
 * For a constant value or no collision need to set constant_ionization to a 
 * value or zero before passing this method to the set_crossSec_ionization
 * method.
 *
 * e.g. 
 * object.set_constant_ionization(1.0e-20);
 * object.set_crossSec_ionization(&e_collisions::sigma_const_ionization)
 */
double e_collisions::sigma_const_ionization(void){
        return constant_ionization;
}

/*
 * Implementation Notes: sigma_Ar_elastic
 * ---------------------------------------
 * Calculates the electron-Argon energy dependent cross-section for
 * energies in the range: threshold to 6 keV. The formulae are a collection
 * excel data fits to data presented in:
 * -Ashok Jain et al,Physica Scripta,vol 41,(1990)321
 */
double e_collisions::sigma_Ar_elastic(void){
   double sigma;

   if ( energy < 0.35 )
  
      sigma = 28019.0*pow(energy,6.0)- 32798.0*pow(energy,5.0)
              + 15048.0*pow(energy,4.0) - 34312.0*pow(energy,3.0)
              + 4084.0*pow(energy,2.0) - 249.0*energy + 7.576;

   else if( energy < 15.0 && energy > 0.35)

      sigma = - 0.013*pow(energy,3.0)
              + 0.232*pow(energy,2.0) + 1.066*energy + 0.085;

   else if( energy  > 15.0 ) //stictly only valid to 6k

      sigma = 126.2*pow(energy,-0.69);

   return sigma;
}

/*
 * Implementation Notes: sigma_Ar_ionize
 * ---------------------------------------
 * Calculates the electron-Argon energy dependent cross-section for
 * energies in the range: threshold to 6 keV. The formulae are a collection
 * excel data fits to data presented in:
 * -Ashok Jain et al,Physica Scripta,vol 41,(1990)321
 *
 */
double e_collisions::sigma_Ar_ionize(void){
   double sigma;
   double epsilon = 1.0e-5;  // use for equality of doubles

   if(energy < 15.755)
      sigma = 0.0;
   
   else if (energy < 100.0 && energy > 15.755) {

       ++countQualifying;
       sigma = (-3.0e-7) * pow(energy, 4.0) + (8.0e-5) * pow(energy, 3.0)
           - 0.008 * pow(energy, 2.0) + 0.416 * energy - 4.758;
   }   

   else if (energy > 100.0)  //only valid to 6k
   {
       ++countQualifying;
       sigma = 130.0 * pow(energy, -0.74);
   }

   return sigma;
}


/*
 * Implementation Notes: sigma_H_elastic
 * --------------------------------------
 * Calculates the electron-Hydrogen energy dependent cross-section for
 * energies in the range: threshold to 1 keV. The formulae are a collection
 * excel data fits to data found in the Aladin Database.
 *
 */
double e_collisions::sigma_H_elastic(void) {
    double sigma;

    if (energy < 1.0e-2) //not sure if extrapolation is correct

        sigma = 3.834 * pow(energy, 2.0) - 20.82 * energy + 40.98;

    else if (energy >= 1.0e-2 && energy < 20.0)

        sigma = -6.15 * log(energy) + 20.91;

    else if (energy >= 20.0)  //strictly only valid up to 1000 eV

        sigma = 102.2 * pow(energy, -1.19);

    return sigma; // (fabs(sigma * 1.0e-20));
}

/*
 * Implementation Notes: sigma_H_ionize
 * ---------------------------------------
 * Calculates the electron-Argon energy dependent cross-section for
 * energies in the range: threshold to 10 keV. The formulae are from:
 * -R K Janev et al, Atomic and Plasma-Material Interaction Data for Fusion,
 * vol 4,(1993)
 */
double e_collisions::sigma_H_ionize(void) {
    double sigma;
    double I = 13.6, I_E = energy / I, C = (1.0 - (1.0 / I_E));
    double A = 0.18450, B1 = -0.032226, B2 = -0.034539, B3 = 1.4003, B4 = -2.8115, B5 = 2.2986;

    if (energy < I)

        sigma = 0.0;

    else if (energy >= I) //valid up to 1e4

        sigma = (1.0 / (I * energy)) * (A * log(I_E) + B1 * C + B2 * pow(C, 2.0)
            + B3 * pow(C, 3.0) + B4 * pow(C, 4.0) + B5 * pow(C, 5.0));

    return sigma; // fabs(sigma * 1.0e-17);
}

/*
 * Implementation Notes: ionization 
 * ----------------------------------
 * Implements ionization collisions using the method outlined in :
 *  -Kenichi Nanbu et al,IEEE Trans. Plasma Sci.,vol 28,(2000) 971
 *   
 * In addition to calculating the post-collision velocities, this method
 * adds the new particles to the particle vectors.
 */   
void e_collisions::ionization(int i, vector<particleVec> &particleObj_e, vector<particleVec> &particleObj_i,double energy, double reduced_mass)
{
   int index,paws=0;
   int last_particle;
   double weight = particleObj_e[i].weight; 
   double eta0,eta1,ejected_energy,R,Wx,Wy,Wz;
   double vxe_after, vye_after, vze_after;
   double v2xe_after, v2ye_after, v2ze_after;
   double ve, vn, ve_after, v2e_after;
   double Rxe, Rye, Rze, R2xe, R2ye, R2ze;
   double excessEnergy;
   double vx_temp_ion, vy_temp_ion, vz_temp_ion, xi, yi; //FOR 2D
   double nvx = 0.0, nvy = 0.0, nvz = 0.0;
   double vxe=0.0,vye=0.0,vze=0.0;
   double xe=particleObj_e[i].x, x2e;
   double ye=particleObj_e[i].y, y2e; // FOR 2D


   //use Nanbu IEEE eqn 44b
   reduced_mass *= ELECTRON_MASS; 

   R = (double)rand()/((double)RAND_MAX);  //use for random splitting of excess energy

   // create a random unit vector for velocity after of impacting electron eqn 52a
   Rxe = 1.0 - 2.0*( (double)rand() / ((double)RAND_MAX));
   Rye = 1.0 - 2.0 * ((double)rand() / ((double)RAND_MAX));
   Rze = 1.0 - 2.0 * ((double)rand() / ((double)RAND_MAX));

   // create a random unit vector for velocity after of impacting electron eqn 52b
   R2xe = 1.0 - 2.0 * ((double)rand() / ((double)RAND_MAX));
   R2ye = 1.0 - 2.0 * ((double)rand() / ((double)RAND_MAX));
   R2ze = 1.0 - 2.0 * ((double)rand() / ((double)RAND_MAX));


  /*The neutral velocities are chosen randomly from an array that is 
   *populated before the simulation starts.
   */
   index= rand()%(INIT_NO_NEUTRALS - 1);

   nvx = neutral_vx[index];
   nvy = neutral_vy[index];
   nvz = neutral_vz[index];
  
   // electron velocity before collision
   vxe = particleObj_e[i].particleVec::vx;
   vye = particleObj_e[i].particleVec::vy;
   vze = particleObj_e[i].particleVec::vz;
  
  /* Center of mass velocity */  // eqn 49
   Wx = (ion_mass*nvx+ELECTRON_MASS*vxe)/(ion_mass+ELECTRON_MASS);
   Wy = (ion_mass*nvy+ELECTRON_MASS*vye)/(ion_mass+ELECTRON_MASS);
   Wz = (ion_mass*nvz+ELECTRON_MASS*vze)/(ion_mass+ELECTRON_MASS);

    //magnitude of the electron and neutral velocities 
   ve = sqrt(vxe * vxe + vye * vye + vze * vze);
   vn = sqrt(nvx * nvx + nvy * nvy + nvz * nvz);

   // equation 51 Nanbu right side aproximately, ve >> vn
   excessEnergy = (reduced_mass / 2.0)*(ve-vn)*(ve-vn) - threshold; 

   ve_after = sqrt(2.0 * R * excessEnergy / ELECTRON_MASS);  // eqn 52a

   v2e_after = sqrt(2.0 * (1.0- R )* excessEnergy / ELECTRON_MASS); // eqn 52b

   // equation 53a velocity of impinging electron after collision
   vxe_after = ve_after * Rxe;
   vye_after = ve_after * Rye;
   vze_after = ve_after * Rze;

   // equation 53b velocity of new ejected electron after collision
   v2xe_after = v2e_after * R2xe;
   v2ye_after = v2e_after * R2ye;
   v2ze_after = v2e_after * R2ze;

   //using Nanbu IEEE article exclusively
    // equation 50 - new ion's velocity after collision
   vx_temp_ion = Wx - (ELECTRON_MASS / ion_mass) * (vxe_after + v2xe_after);

   vy_temp_ion = Wy - (ELECTRON_MASS / ion_mass) * (vye_after + v2ye_after);

   vz_temp_ion = Wz - (ELECTRON_MASS / ion_mass) * (vze_after + v2ze_after);

    // advance the ion's positions
   xi = xe + vx_temp_ion * DTI;
   yi = ye + vy_temp_ion * DTI;

   // create a new ion
   particleObj_i.push_back(particleVec());
   last_particle = particleObj_i.size() - 1;
   particleObj_i[last_particle].x = xi;
   particleObj_i[last_particle].y = yi;
   particleObj_i[last_particle].vx = vx_temp_ion;
   particleObj_i[last_particle].vy = vy_temp_ion;
   particleObj_i[last_particle].vz = vz_temp_ion;
   particleObj_i[last_particle].weight = weight;

    // electron postcollision velocity taking into account the center of mass velocity
   vxe_after += Wx;
   vye_after += Wy;
   vze_after += Wz;

    // ejected electron's velocity
   v2xe_after += Wx;
   v2ye_after += Wy;
   v2ze_after += Wz;

   particleObj_e[i].particleVec::vx = vxe_after;
   particleObj_e[i].particleVec::vy = vye_after;
   particleObj_e[i].particleVec::vz = vze_after;

    // advance ejected electron's position
   x2e = xe + v2xe_after * DT;
   y2e = ye + v2ye_after * DT;
   
   // add ejected electron the electron vector
   particleObj_e.push_back(particleVec());
   last_particle = particleObj_e.size() - 1;
   particleObj_e[last_particle].x = x2e;
   particleObj_e[last_particle].y = y2e;
   particleObj_e[last_particle].vx = v2xe_after;
   particleObj_e[last_particle].vy = v2ye_after;
   particleObj_e[last_particle].vz = v2ze_after;
   particleObj_e[last_particle].weight = weight;
   
}

/*
 * Implementation Notes: collision
 * ---------------------------------
 * Determines if a collision occurs, and calls the collision methods if it 
 * does. This is done in a few steps.
 *
 * Originally followed the following method:
 *1)Determines if a collision occurs by comparing a random number to the 
 *  total collision probability as is described in the following articles:
 *  -V. Vahedi et al,Computer Physics Communications,vol 87,(1995) 179
 *  -C K Birdsall, IEEE Trans. Plasma Sci.,vol 19,(1991) 65
 *  -Kenichi Nanbu et al,IEEE Trans. Plasma Sci.,vol 28,(2000) 971
 * 
 *  Use Birdall, IEEE article's method for the probability
 *
 *2)If a collision occurs the type of collision is determined by comparing 
 *  each individual cross-section to random numbers as described in the
 *  articles listed in step one.
 *
 *3)If a collision occurs the relevant collision method to implement the 
 *  the collision is called with the coll_factor set to the negative
 *  of the reduced mass.
 *  See Kenichi Nanbu and Takizuka for more information about the coll_factor. 
 *  -Tomonori Takizuka et al,Journal of Computational Physics,vol 25,(1977)205
 *
 * 
  *The code which decides whether or not an electron collision occurs
  *is still a work in progress due to the need to limit the number of
  *ionization events. Currently an elastic collision occurs if a random
  *number between zero and one is less than the collision probability and
  *likewise for ionization. In the case of ionization a further condition
  *that the number of ionization events is less than a user defined limit
  *is also enforced.
 * 
 */
void e_collisions::collision(vector<particleVec> &particleObj_e, vector<particleVec> &particleObj_i)
{
   int count_escat=0,count_ionize=0;
   int paws=0;
   int i=0, particle_no, jp, kp;
   char whereabouts= 'n';
   double sigma_eEl, sigma_eIoniz; /*elastic and ionization cross-sections*/
   double sigma_total; /*sum of elastic and ionization cross-sections*/
   double R;
   double reduced_mass,u,coll_factor,ux;
   double offset = 1.0e-6;
   double probability = 0.0;  
   double epsilon = 1.0e-7; // use for = to zero or make > 0
   int countNumInRegion = 0;
   double maxElasticP = 0.0;
   double maxIonizeP = 0.0;
   

   countQualifying = 0;

   particle_no = particleObj_e.size();

   reduced_mass = ion_mass/(ELECTRON_MASS + ion_mass); 

   while(i< particle_no){

       // the grid position of the current electron
       jp = (int)(particleObj_e[i].x / DX);
       kp = (int)(particleObj_e[i].y / DY);

    // if the current electron is inside the domain
    if (jp >=0 && jp< X_NO_OF_CELLS && kp >=0 && kp< Y_NO_OF_CELLS) {

    ++countNumInRegion;

    // get current electron's energy
    calc_energy(particleObj_e[i].vx,particleObj_e[i].vy,particleObj_e[i].vz);

    // get current velocity magnitude
    u= sqrt( (2.0*CHARGE*this->energy)/ELECTRON_MASS);
    
    ux = particleObj_e[i].vx;

   /* Calculate collision cross-sections */
      sigma_eEl = (this->*sigma_elastic)();
      sigma_eIoniz = (this->*sigma_ionization)();
      sigma_total = sigma_eIoniz + sigma_eEl;
      sigma_total+=offset; 

    //probablity a scattering collision will occur
      R = ((double)rand() / ((double)RAND_MAX));  

    // eqn 3 Vahedi et al
    probability = 1.0 - exp(-(u * dt * 1.0e-20 * neutral_density * sigma_eEl));
    
    if (probability > maxElasticP)
        maxElasticP = probability;  // monitor peak magnitude of the probability
    
    // if too many elastic collisions kills the tail in the electron distribution
    // probabilities are a work in progress
    if (probability > 2.0*R)  
    {
        /*elastic collision*/
        coll_factor = -reduced_mass;
        scattered_velocity(i,particleObj_e,this->energy, reduced_mass,YES);
        count_escat += 1;   //used to monitor the number of collisions   
        whereabouts = 's';
        
    }

    // probability of ionizing collision
    probability = 1.0 - exp(-(u * dt * 1.0e-20 * neutral_density * sigma_eIoniz));
    if (probability > maxIonizeP)
        maxIonizeP = probability;
    

    if (probability > 0.0) // very few ionization events occur keep an eye on peak probability
    {
                /*ionizing collision*/
                whereabouts = 'I';
                ionization(i,particleObj_e,particleObj_i,this->energy, reduced_mass);
                count_ionize += 1;  //used to monitor the number of collisions
    } 


   /* Check for over-heating: freezes program and waits for integer input to move on */
    /*
      if(abs(particleObj_e[i].vx) > 7.9e7 && particleObj_e[i].x >=DX &&
           particleObj_e[i].x <= (XLENGTH-DX)){
               cout<<"in collisions, v>7.9e7,v before coll = "<<ux<<" v after= "<<particleObj_e[i].vx<<" x= "<<particleObj_e[i].x<<" weight= "<<particleObj_e[i].weight<<endl;
           cout<<"energy = "<<energy<<endl; 
           cout<<"i = "<<i<<endl;
           cout << " length = " << LENGTH << endl;
           cout << "whereabouts = " << whereabouts << endl;
           cout<<"This is typically due to overheating. Need to ensure the"
                 " electron time-step is less than the 1/plasma_frequency. The"
                 " electron density also needs to be less than the density" 
                 " used in the plasma frequency to help avoid over-heating." 
                " Things to try include reducing the macro-particle size and/or"
                 " the number of allowed ionization events."<<endl;   
           cout<<"Hit enter to move on "<<endl;
           cin>>paws;
        }
        */

     }/* if x within boundaries */
     i+=1;
   } /*goes with while loop*/

    // used to monitor the collision probability
    /*
    cout << " total scattered count = " << count_escat 
        << " total ionization events = "<< count_ionize 
        << " total ionize Qualify = " << countQualifying << endl;
    cout << "num of electrons inside region = " << countNumInRegion;
    cout << "\n max Elastic Prob = " << maxElasticP << " max Ionize Prob = " << maxIonizeP << endl;
    */
}
