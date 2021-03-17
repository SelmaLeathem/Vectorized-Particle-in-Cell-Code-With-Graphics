/*
 * File: boundary.cpp
 * -------------------
 * 
 * 
 * Implementation of the boundary class, which stores variables and 
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
 * Date: 2/28/2021
 * Description: Adjust functions to make calculations on a two-dimensional
 *              grid instead of a one-dimensional line. The original code
 *              above only accounted for two boundaries at x = 0 and x = length.
 *              Now there are two additional boundaries at y = 0 and y = length.
 *              Both boundaries undergo the same events.
 * 
 */

#include "boundary.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>   

boundary::boundary(double dt_in):dt(dt_in){

    //initialize values incase don't get set if don't use all parts of boundary
    wall_offset = 1.0e-6; 
    ion_proportion = 0.0;
    incident_mass = 0.0;
    energy = 0.0;
    angleX = 0.0;
    angleY = 0.0;
    max_yield =0.0;
    E0 = 0.0;
    Emax0 = 0.0;
    yield = 0.0;
    reflection_yield = 0.0;
    backscatter_yield = 0.0;
    threshold = 0.0;
    work_function = 0.0;

    // initialize function pointers to functions that return zero incase
    // not set by the user
    set_e_impact_yield(&boundary::set_e_impact_see_yield_to_zero);
    set_i_impact_yield(&boundary::set_i_impact_see_yield_Ar_to_zero);
    set_reflection_yield(&boundary::set_reflection_yield_to_zero);
    set_backscatter_yield(&boundary::set_backscatter_yield_to_zero);

}

/* Destructor: ~boundary */
boundary::~boundary(){
}

/*
* Method: set_max_yield
* Usage:  set_max_yield(max_yield_value)
* ---------------------------------------
* Sets the maximum secondary emmision coefficient at normal incidence.
*/
void boundary::set_max_yield(double val){
       max_yield=val;
}

/*
* Method: set_E0
* Usage:  set_E0(E0_value)
* ------------------------
* Sets the secondary emmision threshold in eV.
*/
void boundary::set_E0(double val){
       E0=val;
}

/*
* Method: set_Emax0
* Usage:  set_Emax0(Emax0_value)
* -------------------------------
* Sets the energy corresponding to the maximum secondary emmision
* coefficient at normal incidence max_yield.
*/
void boundary::set_Emax0(double val){
       Emax0=val;
}

/*
* Method: set_incident_mass
* Usage:  set_incident_mass(incident_mass_value)
* -----------------------------------------------
* Sets the mass of the particle impinging upon the wall in kg.
*/
void boundary::set_incident_mass(double val){
       incident_mass=val;
}

/*
* Method: set_vt
* Usage:  set_vt(thermal_velocity_value)
* -----------------------------------------------
* Sets the thermal velocity.
*/
void boundary::set_vt(double val){
       vt=val;
}

/*
* Method: set_threshold
* Usage:  set_threshold(threshold_value)
* ---------------------------------------
* Sets the ionization threshold of the atom species in eV. 
*/
void boundary::set_threshold(double val){
       threshold=val;
}

/*
* Method: set_work_function
* Usage:  set_work_function(work_function_value)
* -----------------------------------------------
* Sets the surface work function in eV.
*/
void boundary::set_work_function(double val){
       work_function=val;
}

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
void boundary::set_ion_proportion(double val){
       ion_proportion=val/100.0;
}

/*
* Method: set_implement_e_boundary
* Usage:  set_implement_e_boundary(&boundary::boundary_function);
* --------------------------------------------------------------
* Sets the 'void (boundary::*implement_e_boundary)(...)' function pointer
* to point to the function passed in the arguement.                  
*/ 
void boundary::set_implement_e_boundary(void (boundary::*func1)(vector<particleVec> &electronObj, vector<particleVec> &ionObj)){
      implement_e_boundary = func1;
}

/*
* Method: set_implement_i_boundary
* Usage:  set_implement_i_boundary(&boundary::boundary_function);
* --------------------------------------------------------------
* Sets the 'void (boundary::*implement_i_boundary)(...)' function pointer
* to point to the function passed in the arguement.                  
*/ 
void boundary::set_implement_i_boundary(void (boundary::*func1)(vector<particleVec> &electronObj, vector<particleVec> &ionObj)){
      implement_i_boundary = func1;
}

/* 
* Method: calc_energy
* Usage:  calc_energy(x_velocity,y_velocity,z_velocity);
* ------------------------------------------------------
* Calculates the particle energy in eV for use in yield
* calculations. Note, that user defined yield functions can access
* this value.
*/ 
void boundary::calc_energy(double vx,double vy,double vz){
     double u2;
     u2 = vx*vx + vy*vy +vz*vz;
     energy = (0.5/CHARGE)*incident_mass*u2; //convert to eV units
}

/* 
* Method: calc_angleX
* Usage:  calc_angleX(x_velocity);
* -------------------------------
* Calculates the electron impact angle for use in yield
* calculations. Note, that user defined yield functions can access
* this value.
*/ 
void boundary::calc_angleX(double vx){
    double offset=1.0e-6;
    double u2= 2.0*energy/incident_mass;
 
    angleX= acos(vx/(sqrt(u2)+offset));
}

/* 
* Method: calc_angleY
* Usage:  calc_angleY(y_velocity);
* -------------------------------
* Calculates the electron impact angle for use in yield
* calculations. Note, that user defined yield functions can access
* this value.
*/ 
void boundary::calc_angleY(double vy){
    double offset=1.0e-6;
    double u2= 2.0*energy/incident_mass;
 
    angleY= acos(vy/(sqrt(u2)+offset));
}

/*
* Method: set_e_impact_yield
* Usage:  set_e_impact_yield(&boundary::see_yield_function);
* ----------------------------------------------------------
* Sets the electron e_impact_see_yield function pointer to point to a
* function determined by the user, including a user defined one.
* This function calculates the electron-impact yield.
*
*/
void boundary::set_e_impact_yield(void (boundary::*func1)(void)){
               e_impact_see_yield = func1;
}

/*
* Method: set_i_impact_yield
* Usage:  set_i_impact_yield(&boundary::see_yield_function);
* ----------------------------------------------------------
* Sets the ion i_impact_see_yield function pointer to point to a
* function determined by the user, including a user defined one.
*
*/
void boundary::set_i_impact_yield(void (boundary::*func1)(void)){
                i_impact_see_yield = func1;
}

/*
* Method: set_reflection_yield
* Usage:  set_reflection_yield(&boundary::reflection_yield_function);
* ----------------------------------------------------------
* Sets the electron set_reflection_yield function pointer to point to a
* function determined by the user, including a user defined one.
*
*/
void boundary::set_reflection_yield(void (boundary::*func1)(void)){
                calc_reflection_yield=func1;
}

/*
* Method: set_backscatter_yield
* Usage:  set_backscatter_yield(&boundary::backscatter_yield_function);
* ----------------------------------------------------------
* Sets the electron set_backscatter_yield function pointer to point to a
* function determined by the user, including a user defined one.
*
*/
void boundary::set_backscatter_yield(void (boundary::*func1)(void)){
                calc_backscatter_yield=func1;
}

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
void boundary::e_impact_see_yield_v(void) {
    double w;
    double twoPi = 2.0 * PI;
    double v = (energy - E0) / (Emax0 - E0);
    double kconst = (v > 1) ? 0.25 : 0.62;  // equation 9 from Vaughan

    // equation 14 from Radmilovic

    if (onXboundary){  //if on the x-directed boundary
        w = (energy - E0) / (Emax0 * (1.0 + kconst * (angleX * angleX / twoPi)) - E0);
        yield = max_yield * (1.0 + kconst * (angleX * angleX / twoPi)) * pow((w * exp((1.0 - w))), kconst);
    }
    else  //if on the y-directed boundary
    {
        w = (energy - E0) / (Emax0 * (1.0 + kconst * (angleY * angleY / twoPi)) - E0);
        yield = max_yield * (1.0 + kconst * (angleY * angleY / twoPi)) * pow((w * exp((1.0 - w))), kconst);     
    }
    
   
}

/*
 * Implementation Notes: i_impact_see_yield_Ar
 * ---------------------------------------------
 * Ion impact secondary electron yield for Argon ions impinging a dirty surface.
 * See Also:
 *   -M Radmilovic-Radjenovic et al, Eur. Phys. J. D,Vol 54,(2009),445
 *
 */
void boundary::i_impact_see_yield_Ar(void) {
    yield = 0.0;
    //equation 16 from Radmilovic
    if (energy >= 80.0) {
        yield = 0.006 * energy / (1.0 + (energy / 10.0)) +
            1.05e-4 * (sqrt((energy - 80.0)) / (1.0 + pow((energy / 8000.0), 1.5)));
    }
}

/* Prevent excess at high energies: reflection and backscatter references:
* Ralf Krinke and H M Urbassek,  J. Phys. D: Appl. Phys., vol 29 (1996) 378
* M Radmilovic-Radjenovic and Z L Petrovic, Eur. Phys. J. D, vol 54 (2009) 445
*/
void boundary::calc_reflection_yield_v(void) {
    reflection_yield = 0.0;
    if (energy > 50.0)
        reflection_yield = 0.1 * yield;
}

/* Prevent excess at high energies: reflection and backscatter references:
* Ralf Krinke and H M Urbassek,  J. Phys. D: Appl. Phys., vol 29 (1996) 378
* M Radmilovic-Radjenovic and Z L Petrovic, Eur. Phys. J. D, vol 54 (2009) 445
*/
void boundary::calc_backscatter_yield_v(void) {
    backscatter_yield = 0.0;
    if (energy > 50.0)
        backscatter_yield = 0.1 * yield;
}

/* Set yield functions to return zero in constructor. The user and program then
* only sets the yields they want to use.
*/

// backscatter
void boundary::set_backscatter_yield_to_zero(void)
{
    backscatter_yield = 0.0;
}

// reflection
void boundary::set_reflection_yield_to_zero(void)
{
    reflection_yield = 0.0;
}

//ion impact secondary electron emmision
void boundary::set_i_impact_see_yield_Ar_to_zero(void)
{
    yield = 0.0;
}

//electron impact secondary electron emmision
void boundary::set_e_impact_see_yield_to_zero(void)
{
    yield = 0.0;
}

/*
* Method: reflect_electron
* Usage:  reflect_electron(particle_number,particle_vector)
* --------------------------------------------------------------
* Reflects the particle_vector element identified by the particle_number.
*
*/
void boundary::reflect_electron(int current_particle,vector<particleVec> &particleObj){
    int jp, kp;
    double R;
    double lengthX = DX * X_NO_OF_CELLS; // length in x direction
    double lengthY = DY * Y_NO_OF_CELLS; // length in y direction
    R = (double)(rand())/((double)RAND_MAX);

   // get the current grid location
   jp = (int)(particleObj[current_particle].x / DX);
   kp = (int)(particleObj[current_particle].y / DY);

    // right boundary
    if (jp >= X_NO_OF_CELLS) {
        // reverse velocity
        particleObj[current_particle].vx = -particleObj[current_particle].vx;
        // advance position from surface by a random amount
        particleObj[current_particle].x = lengthX -
        R*dt*abs(particleObj[current_particle].vx);
    }
    else if (jp < 0) // left boundary
    {
        // reverse velocity
        particleObj[current_particle].vx = -particleObj[current_particle].vx;
        // advance position from surface by a random amount
        particleObj[current_particle].x = R*dt*abs(particleObj[current_particle].vx);
    }
        
    
    if (kp >= Y_NO_OF_CELLS) { // top boundary
        // reverse velocity
        particleObj[current_particle].vy = -particleObj[current_particle].vy;
        // advance position from surface by a random amount
        particleObj[current_particle].y = lengthY -
        R*dt*abs(particleObj[current_particle].vy);
    }
    else if (kp < 0) // bottom boundary
    {
        // reverse velocity
        particleObj[current_particle].vy = -particleObj[current_particle].vy;
        // advance position from surface by a random amount
        particleObj[current_particle].y = R*dt*abs(particleObj[current_particle].vy);
    }
      
}

/*
 * Implementation Notes: secondary_vel
 * --------------------------------------
 * Calculates the position and velocity of backscattered and true secondary
 * electrons using the formulae presented in:
 *   -Ralf Krimke et al, J. Phys. D: Appl. Phys.,vol 29,(1996),378
 *
 */
void boundary::secondary_vel(int current_particle,double velocity,vector<particleVec> &particleObj){
   int paws;
   int jp, kp;
   double R,Rphi,Rchi,phi,chi; // emmision angles from surface
   double lengthX = DX * X_NO_OF_CELLS; // length in x direction
   double lengthY = DY * Y_NO_OF_CELLS; // length in y direction

   // get the current grid location
   jp = (int)(particleObj[current_particle].x / DX);
   kp = (int)(particleObj[current_particle].y / DY);

    // left-right boundaries
    if (jp < 0 || jp >= X_NO_OF_CELLS)
    {
        R = (double)(rand())/((double)RAND_MAX); 
        Rphi = (double)(rand())/((double)RAND_MAX);
        Rchi = (double)(rand())/((double)RAND_MAX);
        chi=asin(sqrt(Rchi));  //equation 5 from Krimke
        phi= 2.0*PI*Rphi; //equation 6 from Krimke

        //equations 8-11 from Krimke
        particleObj[current_particle].vx = abs(velocity*cos(chi));
        particleObj[current_particle].vy = velocity*sin(chi)*cos(phi);
        particleObj[current_particle].vz = velocity*sin(chi)*sin(phi);
        particleObj[current_particle].x = R*dt*abs(particleObj[current_particle].vx);

        if (jp >= X_NO_OF_CELLS)
        {
            // account for coming in in -x direction
            particleObj[current_particle].x = lengthX - particleObj[current_particle].x;
            particleObj[current_particle].vx = -abs(velocity*cos(chi));
        }
            
        
    }

    // top-bottom boundaries
    if (kp < 0 || kp >= Y_NO_OF_CELLS) 
    {
        R = (double)(rand())/((double)RAND_MAX);
        Rphi = (double)(rand())/((double)RAND_MAX);
        Rchi = (double)(rand())/((double)RAND_MAX);
        chi=asin(sqrt(Rchi));   //equation 5 from Krimke
        phi= 2.0*PI*Rphi;   //equation 6 from Krimke

        //equations 8-11 from Krimke
        particleObj[current_particle].vy = abs(velocity*cos(chi));
        particleObj[current_particle].vx = velocity*sin(chi)*cos(phi);
        particleObj[current_particle].vz = velocity*sin(chi)*sin(phi);
        particleObj[current_particle].y = R*dt*abs(particleObj[current_particle].vy);

        if (kp >= Y_NO_OF_CELLS)
        {
            // account for coming in in -y direction
            particleObj[current_particle].y = lengthY - particleObj[current_particle].y;
            particleObj[current_particle].vy = -abs(velocity*cos(chi));
        }           
    }

}

/*
 * Implementation Notes: secondary_electrons
 * ---------------------------------------------
 * Electron secondary velocity, 'velocity' corresponds to 'v' in:
 *   -Ralf Krimke et al, J. Phys. D: Appl. Phys.,vol 29,(1996),378
 */
void boundary::secondary_electrons(int current_particle,vector<particleVec> &particleObj){
   int last_particle;
   int jp;
   double R; 
   double velocity;

    // get the current grid location
   jp = (int)(particleObj[current_particle].x / DX);

   particleObj.push_back(particleVec());
   last_particle  = particleObj.size()-1;
   //initialize values
   particleObj[last_particle].x = 0.0;
   particleObj[last_particle].y = 0.0;
   particleObj[last_particle].vx = 0.0;
   particleObj[last_particle].vy = 0.0;
   particleObj[last_particle].vz = 0.0;
   particleObj[last_particle].weight = 1.0;

    //use the weight of the incident particle
    particleObj[last_particle].weight = particleObj[current_particle].weight;

    // get x or y position depending upon which boundary are at, if reach this
    // function at the particle is not at the x boundary then it must be at
    // a y boundary
    if (jp < 0 || jp >= X_NO_OF_CELLS)
        particleObj[last_particle].x = particleObj[current_particle].x;
    else
        particleObj[last_particle].y = particleObj[current_particle].y;

    R = (double)(rand())/((double)RAND_MAX);

    //equation 7 from Krimke
    velocity= sqrt((2.0*R*CHARGE*(threshold-2.0*work_function)/ELECTRON_MASS));
    secondary_vel(last_particle,velocity,particleObj);
}

/*
* Method: absorb
* Usage:  absorb(particle_vector)
* ---------------------------------
* Implements an absorbing boundary, which is quicker than setting the
* yields to zero. Currently only used for debugging.
* 
*/ 
void boundary::absorb(vector<particleVec> &particleObj){
     int i=0;
     int n = particleObj.size();
     int jp, kp;

    // go through each particle if it is at a boundary delete it from the
    // list of particles
     while( i < n){

         // get the current grid location
         jp = (int)(particleObj[i].x / DX);
         kp = (int)(particleObj[i].y / DY);
         
        // if at left-right boundary
         if (jp < 0 ||jp >= X_NO_OF_CELLS) {

            // swap before delete to prevent vector shuffle
            swap(particleObj[i],particleObj[(n-1)]);
            particleObj.pop_back();
              
         } // if at top-bottom boundary
         else if (kp < 0 ||kp >= Y_NO_OF_CELLS)
         {
            swap(particleObj[i],particleObj[(n-1)]);
            particleObj.pop_back();
         }
         else
           i++;
        n = particleObj.size();
     
     }
}


 /*
* Method: e_boundary
* Usage:  e_boundary(e_particle_vector,ion_particle_vector)
* ------------------------------------------------------
* Implements the electron particle boundary. Electrons are either
* reflected, backscattered, absorbed, or create secondary electrons.
*
*/
void boundary::e_boundary(vector<particleVec>& particleObj, vector<particleVec>& ionsObj)
{
    int sign = 1, paws = 1, no_inject = 0; // number of electrons to inject
    // count the number of each type of surface interaction
    int i = 0, j, count_ref = 0, count_scat = 0, count_sec = 0, count_absorb = 0;
    int jp, kp;
    int n = particleObj.size();
    double R = 0.0, velocity;
    double yield_max = 0.0;
    double offset = 1.0e-6;

    while (i < n) {       // go through all electrons to see if at boundary 

         // get the current grid location
        jp = (int)(particleObj[i].x / DX);
        kp = (int)(particleObj[i].y / DY);

        //if at a grid boundary
        if (jp < 0 || jp >= X_NO_OF_CELLS || kp < 0 || kp >= Y_NO_OF_CELLS) {

            if (jp < 0 || jp >= X_NO_OF_CELLS)
                onXboundary = true;     // if at left-right boundary
            else
                onXboundary = false;    // if at top-bottom boundary

            // calculate the incident energy
            calc_energy(particleObj[i].vx, particleObj[i].vy, particleObj[i].vz);

            // calculate the incident angle
            calc_angleX(particleObj[i].vx);
            calc_angleY(particleObj[i].vy);

            // get the electron impact secondary yield
            (this->*e_impact_see_yield)();
            R = (double)rand() / ((double)RAND_MAX);

            if (this->yield > yield_max) yield_max = this->yield;

            if (R < this->yield) { // if a random number < the yield then this event occurs

                (this->*calc_reflection_yield)();  // get the reflection yield
                (this->*calc_backscatter_yield)();  // get the backscatter yield

                // estimate number of electrons to inject from surface based on these yields
                no_inject = (int)floor(abs((this->yield - this->reflection_yield - this->backscatter_yield))) + 1;

                R = (double)rand() / ((double)RAND_MAX);

                // does reflection occur?
                if (R < this->reflection_yield / (this->yield + offset)) {
                    reflect_electron(i, particleObj); // implement reflection
                    count_ref++;
                    i++;
                }

                // does backscattering occur?
                else if (R > this->reflection_yield / this->yield && R < (this->reflection_yield + this->backscatter_yield) / (this->yield + offset)) {
                    R = (double)(rand() % 10000) / 10000.0;

                    //implement backscattering
                    velocity = R * sqrt((2.0 * CHARGE * this->energy) / ELECTRON_MASS);
                    secondary_vel(i, velocity, particleObj);
                    count_scat++;
                    i++;
                }
                // does secondary emmission occur?
                else if (R > (this->reflection_yield + this->backscatter_yield) / (this->yield + offset)) {
                    for (j = 0; j < no_inject; j++) {

                        //implement secondary emmision
                        secondary_electrons(i, particleObj);
                        n = particleObj.size();
                    }
                    swap(particleObj[i], particleObj[(n - 1)]);
                    particleObj.pop_back();  // remove incident particle
                    //i--;
                    count_sec++;
                }
            } //if R < yield 

            else {
                // particle gets absorbed by the surface
                swap(particleObj[i], particleObj[(n - 1)]);
                particleObj.pop_back();  // remove particle
                //i--;
                count_absorb++;
            }

        } // end if particle outside boundary 
        else {
            i++;  // particle is within the boundary so move on
        }

        n = particleObj.size();

    } //while loop
}

/*
 * Implementation Notes: i_seconondary_vel
 * -----------------------------------------
 * Calculates the position and velocity of backscattered and true secondary
 * electrons using the formulae presented in:
 *   -Ralf Krimke et al, J. Phys. D: Appl. Phys.,vol 29,(1996),378
 *
 */
void boundary::i_secondary_vel(int last_particle,vector<particleVec> &electrons,double x, double y, double weight){
   int no_inject=0, jp, kp;
   double R=0.0;
   double velocity;
   double Rphi,Rchi,phi,chi;    // emmision angles from surface
   double lengthX = DX * X_NO_OF_CELLS; // length in x direction
   double lengthY = DY * Y_NO_OF_CELLS; // length in y direction

    // get the current grid location
    jp = (int)(x / DX);
    kp = (int)(y / DY);

    if (jp < 0 || jp >= X_NO_OF_CELLS) // left-right boundaries
    {
        R = (double)(rand())/((double)RAND_MAX);

        //equation 7 from Krimke
        velocity= sqrt((2.0*R*CHARGE*(threshold-2.0*work_function)/ELECTRON_MASS));

        Rphi = (double)(rand())/((double)RAND_MAX);
        Rchi = (double)(rand())/((double)RAND_MAX);
        
        chi=asin(sqrt(Rchi)); //equation 5 from Krimke
        phi= 2.0*PI*Rphi; //equation 6 from Krimke

        //equations 8-11 from Krimke
        electrons[last_particle].vx = abs(velocity*cos(chi));
        electrons[last_particle].vy = velocity*sin(chi)*cos(phi);
        electrons[last_particle].vz = velocity*sin(chi)*sin(phi);
        electrons[last_particle].weight= weight;
        electrons[last_particle].x = R*dt*abs(electrons[last_particle].vx);

        if (jp >= X_NO_OF_CELLS)
        {
            // account for coming in in -x direction
            electrons[last_particle].x = lengthX -electrons[last_particle].x;
            electrons[last_particle].vx = -abs(velocity*cos(chi));
        }
            
   }

   if (kp < 0 || kp >= Y_NO_OF_CELLS)  // top-bottom boundaries
   {
        R = (double)(rand())/((double)RAND_MAX);

        //equation 7 from Krimke
        velocity= sqrt((2.0*R*CHARGE*(threshold-2.0*work_function)/ELECTRON_MASS));

        Rphi = (double)(rand())/((double)RAND_MAX);
        Rchi = (double)(rand())/((double)RAND_MAX);
        
        chi=asin(sqrt(Rchi));  //equation 5 from Krimke
        phi= 2.0*PI*Rphi;  //equation 6 from Krimke

        //equations 8-11 from Krimke
        electrons[last_particle].vy = abs(velocity*cos(chi));
        electrons[last_particle].vx = velocity*sin(chi)*cos(phi);
        electrons[last_particle].vz = velocity*sin(chi)*sin(phi);
        electrons[last_particle].weight= weight;
        electrons[last_particle].y = R*dt*abs(electrons[last_particle].vy);

        if (kp >= Y_NO_OF_CELLS)
        {
            // account for coming in in -y direction
            electrons[last_particle].y = lengthY -electrons[last_particle].y;
            electrons[last_particle].vy = -abs(velocity*cos(chi));
        }           
   }
    
}

/*
* Method: ion_boundary
* Usage:  ion_boundary(e_particle_vector,ion_particle_vector)
* --------------------------------------------------------------
* Implements the ion particle boundary. All ions are absorbed at the
* surface, and some create secondary electrons prior to absorption.
*/
void boundary::ion_boundary(vector<particleVec> &electrons,vector<particleVec> &ions){
   int no_inject=0,paws=1;
   int i=0,j,last_particle, jp, kp;
   int n = ions.size();
   double R;

   while( i < n){     // go through all ions to see if at boundary   

        // get the current grid location
        jp = (int)(ions[i].x / DX);
        kp = (int)(ions[i].y / DY);

        //if at a grid boundary
        if(jp < 0 || jp >= X_NO_OF_CELLS || kp < 0 || kp >= Y_NO_OF_CELLS){
            
            // calculate the incident energy
            calc_energy(ions[i].vx,ions[i].vy,ions[i].vz);

            // get the ion impact secondary yield
            (this->*i_impact_see_yield)();
            
            R = (double)rand() / ((double)RAND_MAX);
            
            // if a random number < the yield then this event occurs
            if ( R < this->yield ){
            
            // estimate number of electrons to inject from surface based on this yield
             no_inject = (int)floor(this->yield)+1;

             for (j=0;j<no_inject; j++){

                 // add newly injected electrons
                electrons.push_back(particleVec());
                last_particle  = electrons.size()-1;

                // initialize new particle 
                electrons[last_particle].x = 0.0;
                electrons[last_particle].y = 0.0;
                electrons[last_particle].vx = 0.0;
                electrons[last_particle].vy = 0.0;
                electrons[last_particle].vz = 0.0;
                electrons[last_particle].weight = 1.0;

                // calculate electron's properties
                i_secondary_vel(last_particle, electrons, ions[i].x, ions[i].y, ions[i].weight);
             }
          }  //if R<yield
                 
          swap(ions[i],ions[(n-1)]); // remove the incident ion
          ions.pop_back();

      } //if x outside boundary 
      else{
         i++;
      }
      n=ions.size();
   
   } //while i<n

}

/*
 * Implementation Notes: get_particle_velocity
 * ---------------------------------------------
 * This method is a temporary place holder for a more accurate way to
 * calculate the particle velocity taken randomly off a particle velocity
 * distribution. In this case Maxwellian.
 */
double boundary::get_particle_velocity(double thermal_velocity)
{
    double R1, R2;
    double v = 0.0;

    do{
       R1 =(double)rand()/((double)RAND_MAX);
       R2 =(double)rand()/((double)RAND_MAX);
         v = thermal_velocity*cos(2.0*PI*R2)*sqrt(abs(log(R1)));
    }while( abs(v) > thermal_velocity*10000);

    return v;
}

/*
* Method: plop_particles
* Usage:  plop_particles(elecElement,ionElement,vte,vti);
* ------------------------------------------------
* Plops an ion-electron pair into the simulation domain.
*
*/
void boundary::plop_particles(particleVec &elecElement, particleVec &ionElement,double vte, double vti){
    double random_numberX;
    double random_numberY;

    // Select a random location to plop each particle. Each particle should
    // be plopped at the same location.
    random_numberX = (double)rand() / ((double)RAND_MAX);
    random_numberY = (double)rand() / ((double)RAND_MAX);

    elecElement.x = (DX + random_numberX * (XLENGTH - 2.0 * DX));
    ionElement.x = (DX + random_numberX * (XLENGTH - 2.0 * DX));
    elecElement.y = (DY + random_numberY * (YLENGTH - 2.0 * DY));
    ionElement.y = (DY + random_numberY * (YLENGTH - 2.0 * DY));

    // get the particle's velocity components based on their thermal velocity

    elecElement.vx = get_particle_velocity(vte);
    ionElement.vx = get_particle_velocity(vti);

    elecElement.vy = get_particle_velocity(vte);
    ionElement.vy = get_particle_velocity(vti);

    elecElement.vz = get_particle_velocity(vte);
    ionElement.vz = get_particle_velocity(vti);
}

/*
 * Implementation Notes: steady_e_boundary
 * ---------------------------------------------
 * Implements a type of boundary that results in a rapid steady state, by
 * inserting an electron-ion pair in the simulation domain everytime an
 * electron arrives at the surface. Ionization should be turned off and
 * the ion time-step should equall the electron value. The particle number
 * is constant so the desired number of simulation particles should be given
 * the particle::init_particle_no parameter.
 */  
void boundary::steadyBoundaryX(vector<particleVec> &electronObj, vector<particleVec> &ionObj) {
         int i, j, k, n = electronObj.size();
         int n_ion = ionObj.size();
         int count_LHS = 0, i_before;
         int count_RHS = 0;
         int amount_to_plopLHS;  // number of electrons available to plop on the left side
         int amount_to_plopRHS;   // number of electrons available to plop on the right side
         int klast_lhs = 0, klast_rhs = 0;
         int no_plopped = 0;
         int jp, jpi;

        // This was originally coded with multiple ion species in mind where it is
        // possible to have more ions at the boundaries than free electrons from the
        // original vector. Note: that the vector lengths are constant for this boundary
        // type

         for (i = 0; i < n; i++) {
             jp = (int)(electronObj[i].x / DX); // get the current grid location

             if (electronObj[i].x < 0.0) // if at left boundary
                 count_LHS++;

             else if (jp >= X_NO_OF_CELLS) // if at right boundary
                 count_RHS++;
         }

        // number of available electrons to plop at each boundary
         amount_to_plopLHS = (int)((ion_proportion) * ((double)count_LHS));
         amount_to_plopRHS = (int)((ion_proportion) * ((double)count_RHS));
         i = 0;

         while (i < n) {  // go through all the particles

             jp = (int)(electronObj[i].x / DX);

             if (electronObj[i].x < 0.0) {  // if an electron is outside of the left boundary

                 k = klast_lhs;
                 i_before = i;
                 while (amount_to_plopLHS > 0 && k < n_ion) { // while there are enough electrons to plop

                     jpi = (int)(ionObj[k].x / DX);

                     if (ionObj[k].x < 0.0) // if an ion is at the left boundary
                     {
                         //plop the electron-ion pair into the domain
                         plop_particles(electronObj[i], ionObj[k], VTE, VTI);
                         amount_to_plopLHS--;
                         i++;
                         no_plopped++;
                         break; //break out of while
                     }
                     k++;
                 } //while k>

                 klast_lhs = k + 1;
                 if (i == i_before) {
                     i++;
                 }
             }
             else if (jp >= X_NO_OF_CELLS) { // if an electron is outside of the right boundary

                 k = klast_rhs;
                 i_before = i;
                 while (amount_to_plopRHS > 0 && k < n_ion) { // while there are enough electrons to plop

                     jpi = (int)(ionObj[k].x / DX);

                     if (jpi >= X_NO_OF_CELLS)  // if an ion is at the right boundary
                     {
                         //plop the electron-ion pair into the domain
                         plop_particles(electronObj[i], ionObj[k], VTE, VTI);
                         amount_to_plopRHS--;
                         i++;
                         no_plopped++;
                         break; //break out of while
                     }
                     k++;
                 } //while k> 

                 klast_rhs = k + 1;
                 if (i == i_before) {
                     i++;
                 }
             }
             else {

                 i++;
             }

         } //while i<n

     cout<<" no plopped in x = "<<no_plopped<<endl;
}

/*
 * Implementation Notes: steady_e_boundary
 * ---------------------------------------------
 * Implements a type of boundary that results in a rapid steady state, by
 * inserting an electron-ion pair in the simulation domain everytime an
 * electron arrives at the surface. Ionization should be turned off and
 * the ion time-step should equall the electron value. The particle number
 * is constant so the desired number of simulation particles should be given
 * the particle::init_particle_no parameter.
 */  
void boundary::steadyBoundaryY(vector<particleVec> &electronObj,vector<particleVec> &ionObj){
     int i, j, k, n = electronObj.size();
     int n_ion = ionObj.size();
     int count_BOTTOM = 0, i_before;
     int count_TOP = 0;
     int amount_to_plopBOTTOM; // number of electrons available to plop at the bottom
     int amount_to_plopTOP;  // number of electrons available to plop at the top
     int klast_lhs = 0, klast_rhs = 0;
     int no_plopped = 0;
     int jp, jpi;

     for (i = 0; i < n; i++) {
         jp = (int)(electronObj[i].y / DY);  // get the current grid location

         if (electronObj[i].y < 0.0)    // if at bottom boundary
             count_BOTTOM++;

         else if (jp >= Y_NO_OF_CELLS) // if at top boundary
             count_TOP++;
     }

     amount_to_plopBOTTOM = (int)((ion_proportion) * ((double)count_BOTTOM));
     amount_to_plopTOP = (int)((ion_proportion) * ((double)count_TOP));
     i = 0;

     while (i < n) {    // go through all the particles

         jp = (int)(electronObj[i].y / DY);

         if (electronObj[i].y < 0.0) {  // if an electron is outside of the bottom boundary

             k = klast_lhs;
             i_before = i;
             while (amount_to_plopBOTTOM > 0 && k < n_ion) { // while there are enough electrons to plop

                 jpi = (int)(ionObj[k].y / DY);

                 if (ionObj[k].y < 0.0)  // if an ion is at the bottom boundary
                 {
                     //plop the electron-ion pair into the domain
                     plop_particles(electronObj[i], ionObj[k], VTE, VTI);
                     amount_to_plopBOTTOM--;
                     i++;
                     no_plopped++;
                     break; //break out of while
                 }
                 k++;
             } //while k>

             klast_lhs = k + 1;
             if (i == i_before) {
                 i++;
             }
         }
         else if (jp >= Y_NO_OF_CELLS) { // if an electron is outside of the top boundary

             k = klast_rhs;
             i_before = i;
             while (amount_to_plopTOP > 0 && k < n_ion) { // while there are enough electrons to plop

                 jpi = (int)(ionObj[k].y / DY);

                 if (jpi >= Y_NO_OF_CELLS) // if an ion is at the top boundary
                 {
                     //plop the electron-ion pair into the domain
                     plop_particles(electronObj[i], ionObj[k], VTE, VTI);
                     amount_to_plopTOP--;
                     i++;
                     no_plopped++;
                     break; //break out of while
                 }
                 k++;
             } //while k> 

             klast_rhs = k + 1;
             if (i == i_before) {
                 i++;
             }
         }
         else {

             i++;
         }

     } //while i<n

     cout<<" no plopped in y = "<<no_plopped<<endl;
}

/*
 * Implementation Notes: steady_e_boundary
 * ---------------------------------------------
 * Implements a type of boundary that results in a rapid steady state, by
 * inserting an electron-ion pair in the simulation domain everytime an
 * electron arrives at the surface. Ionization should be turned off and
 * the ion time-step should equall the electron value. The particle number
 * is constant so the desired number of simulation particles should be given
 * the particle::init_particle_no parameter.
 */  
void boundary::steady_e_boundary(vector<particleVec> &electronObj,vector<particleVec> &ionObj){
        steadyBoundaryX(electronObj,ionObj);  // left-right boundaries
        steadyBoundaryY(electronObj,ionObj);  // top-bottom boundaries
}

void boundary::steady_i_boundary(vector<particleVec> &electronObj,vector<particleVec> &ionObj){
    /* This is intentionally blank,everything is done in the electron boundary
     * method.
     */ 
}


