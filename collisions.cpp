/*
 * File: collisions.cpp
 * ---------------------
 * 
 * 
 * This file implements the collisions.h interface.
 *
 * 6/15
 * 
 *  Date: 2/28/2021
 * 
 * The scattering function defined in this class currently only applies to electrons
 * although is expandable to include ions. This class is kept as a parent to 
 * e_collisions and ion_collisions since at some future date it might provide abstract
 * function headers should ion scattering additionaly be accounted for.
 */

#include "collisions.h"

/*
* Constructor: collisions
* Usage: collisions coll_object;
*        collisions coll_object(is_an_electron,dt);  
*/ 
collisions::collisions(int if_electron_in,double dt_in):if_electron(if_electron_in),dt(dt_in)
{

   neutral_vx = new double[INIT_NO_NEUTRALS]();
   if (!neutral_vx){
      printf("Error allocating memory for neutral_vx\n");
      exit(1);
   }
   
   neutral_vy = new double[INIT_NO_NEUTRALS]();
   if (!neutral_vy){
      printf("Error allocating memory for neutral_vy\n");
      exit(1);
   }

   neutral_vz = new double[INIT_NO_NEUTRALS]();
   if (!neutral_vz){
      printf("Error allocating memory for neutral_vz\n");
      exit(1);
   }

}

/* Destructor: ~collisions */
collisions::~collisions()
{
    delete [] neutral_vx;
    delete [] neutral_vy;
    delete [] neutral_vz;
}

/* 
 * Implementation Notes: calc_neutral_velocity
 * --------------------------------------------
 * This is a place-holder for a more accurate way of taking
 * velocities off a Maxwellian. 
 *
 * Arrays of neutral velocities are calculated before the simulation starts.
 * Neutral velocities are randomly chosen from these arrays during the 
 * simulation in collision calculations.
 */
void collisions::calc_neutral_velocity(void)
{
   int i;
   double R1,R2;

   for (i=0; i< INIT_NO_NEUTRALS; i++)
   {
    do{
       R1 = (double)rand()/((double)RAND_MAX);
       R2 = (double)rand()/((double)RAND_MAX);

       neutral_vx[i] = vtn*cos(PI*R2)*sqrt(abs(log(R1)));
    }while(abs(neutral_vx[i]) > 10000.0*vtn );

   do{

      R1 = (double)rand()/((double)RAND_MAX);
      R2 = (double)rand()/((double)RAND_MAX);
      neutral_vy[i] = vtn*cos(PI*R2)*sqrt(abs(log(R1)));
      neutral_vz[i] = vtn*sin(PI*R2)*sqrt(abs(log(R1)));
    }while(abs(neutral_vy[i]) > 10000.0*vtn || abs(neutral_vz[i]) > 10000.0*vtn);
   }

}

/*
 * Implementation Notes: scattered_velocity
 * ------------------------------------------
 * Calculates the post-collision velocites after a scattering event, as 
 * presented in:
 *  -Kenichi Nanbu,IEEE Trans. Plasma Sci.,vol 28,(2000)971
 *  -Tomonori Takizuka et al,Journal of Computational Physics,vol 25,(1977)205
 * 
 * This method is called for elastic scattering postcollision velocities, and
 * post ionization scattering.
 *
 *  - V. Vahedi and M. Surendra in Computer Physics Communications, vol 87 (1995) 179
 *  -Kenichi Nanbu,IEEE Trans. Plasma Sci.,vol 28,(2000)971
 *  for cos/sinChi and cos/sinPhi
 * Note: Energy is in eV.
 */
void collisions::scattered_velocity(int i,vector<particleVec> &particleObj, double energy, double coll_factor,int yes_or_no)
{
   int index;
   double R,cosChi,sinChi,cosPhi,sinPhi;  //scattering angles
   double nvx=0.0,nvy=0.0,nvz=0.0;
   double ux,uy,uz,u,uperp, uyp, uzp;  //intermediatory variables
   double two_pi = 2.0*PI;
   double min_vel = 1.0;
   double minEnergy=1.0e-2;  //to avoid division by zero
   

   R = (double)rand()/((double)RAND_MAX);
   
   // eqn 55 in Nanbu IEEE
   if (energy < minEnergy)
       cosChi = sqrt(1.0 - R);
   else
        cosChi = 2.0 * (1.0 - pow((1.0 + energy), R)) / energy + 1.0;

   if (cosChi > 1.0)
       cout << "cosChi > 1\n";

   sinChi = sqrt(1.0 - cosChi*cosChi);

   //eqn 42 in Nanbu IEEE
   cosPhi = cos(two_pi*R);
   sinPhi = sin(two_pi*R);  

   index= rand()%(INIT_NO_NEUTRALS - 1);
  
   nvx = neutral_vx[index]*((double)yes_or_no);
   nvy = neutral_vy[index]*((double)yes_or_no);
   nvz = neutral_vz[index]*((double)yes_or_no);

   // g = v - v(neutral) in eqn 56 Nanbu IEEE
   ux = particleObj[i].particleVec::vx - nvx;
   uy = particleObj[i].particleVec::vy - nvy;
   uz = particleObj[i].particleVec::vz - nvz;

   // perpendicular velocity
   uperp = sqrt(uy*uy + uz*uz);

   // total velocity
   u = sqrt(ux * ux + uperp * uperp);

   if (uperp < min_vel)  //avoid division by zero
   {
         uyp = 1.0;
         uzp = 1.0;
   }
   else
   {
         uyp = uy / uperp;
         uzp = uz / uperp;
   }

   // from Nanbu IEEE eqns 56

   particleObj[i].particleVec::vz -= coll_factor * (uz*(1.0 - cosChi) - (ux*uzp*cosPhi - u*uyp*sinPhi)*sinChi );

   particleObj[i].particleVec::vy -= coll_factor * (uy * (1.0 - cosChi) - (ux*uyp*cosPhi  + u*uzp*sinPhi)*sinChi);

   particleObj[i].particleVec::vx -= coll_factor * (ux * (1.0 - cosChi) + uperp*cosPhi*sinChi );
  
}
