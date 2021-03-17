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
 *       Multigrid solves Poission equation for the potential on a two-dimensional grid 
 *       in O(n) as compared to O(n^3) for straight LU decomposition.
 */

#include "fields.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

/*
* Constructor: fields
* Usage: fields field_object;
*
*/
fields::fields(void)
{
   int i, j;
   double zero = 0.0;

   // charge density
   rho = new double*[Y_NO_OF_CELLS]();
   if (!rho){
      printf("Error allocating memory for rho\n");
      exit(1);
   }

   for(i=0;i<(Y_NO_OF_CELLS);i++){
      rho[i] = new double[X_NO_OF_CELLS]();
      if(!rho[i]){
         printf("Error allocating memory for rho\n");
         exit(1);
      }
      for(j=0;j<(X_NO_OF_CELLS);j++){
         rho[i][j] = 0.0;
      }
   }

   /* Electrostatic potential is initiated using routines
    * from nr3.h within mglin.h (Numerical Recipes Code)
    */
   phi.assign(Y_NO_OF_CELLS,X_NO_OF_CELLS,zero);

   /* electric field in the x direction is an OpenCL type double */
   electric_fieldX = new cl_double[(X_NO_OF_CELLS + 1)* (Y_NO_OF_CELLS + 1)]();
   if (!electric_fieldX) {
       printf("Error allocating memory for electric_fieldX\n");
       exit(1);
   }

    /* electric field in the y direction is an OpenCL type double */
   electric_fieldY = new cl_double[(X_NO_OF_CELLS + 1) * (Y_NO_OF_CELLS + 1)]();
   if (!electric_fieldY) {
       printf("Error allocating memory for electric_fieldY\n");
       exit(1);
   }

}

/* Destructor: ~fields */
fields::~fields(){

   for(int i=0;i<Y_NO_OF_CELLS;i++)
   {
      delete [] rho[i];
   }

   delete [] rho;
   delete [] electric_fieldX;
   delete [] electric_fieldY;
}

/*
* Method: grid_rho
* Usage:  grid_rho(ion_objects,e_grid_charge);
* ----------------------------------------------
* Calculates the net charge density on each grid point.
*
*/
void fields::grid_rho(particle *ion, cl_double *qe)
{
   int i,j;

   /*
     Rho = electron charge denisty + ion charge density
     where it is understood that the electron charge denisty
     has a negative value since the electron charge is negative
     while the ion charge density is positive

   */

  // Initialize rho with the ion charge density
   for (j = 0; j < Y_NO_OF_CELLS; j++)
   {
      for (i = 0; i < X_NO_OF_CELLS; i++)
      {
         rho[j][i] = ion[0].get_q(j,i);
      }
   }

   // If there are other ion species then add their charge density
   // to rho
   for(int k =1; k<NO_ION_SPECIES;k++)
   {
      for (j = 0; j < Y_NO_OF_CELLS; j++)
      {
         for (i = 0; i < X_NO_OF_CELLS; i++)
         {
            rho[j][i] += ion[k].get_q(j,i);
         }
      }
   }

   // Add the electron charge density to rho
   for (j = 0; j < Y_NO_OF_CELLS; j++)
   {
      for (i = 0; i < X_NO_OF_CELLS; i++)
      {
         rho[j][i]  += (double)qe[j * X_NO_OF_CELLS + i];
      }
   }

} 

/* 
 * Implementation Notes: grid_potential
 * ---------------------------------------
 * Calculate the electrostatic potential at each grid point from the finite-
 * difference Poisson equation using the method outlined in:
 * C K Birdsall and A B Langdon, Plasma Physics via Computer Simulation,
 * Chapter 2 and Appendix D.
 *
 * The resulting matrix for phi is solved in O(N) using the multigrid
 * method. The code for the multigrid method is imported from Numerical
 * Recipes below. 
 * 
 * The grid must be sqare with (n+1) x (n+1) dimension size, where n
 * is a power of two. 
 * 
 * Use Numerical Recipes: The art of scientific computing by William Press
 * 
 * Boundary Conditions:
 * phi0X = phi value of left surface
 * phiLX = phi value of right surface
 * phi0Y = phi value of bottom surface
 * phiLY = phi value of top surface
 * 
 */
void fields::grid_potential(double phi0X, double phiLX, double phi0Y, double phiLY, double timestep, double phi_period)
{  
   int i,j;
   double coef = 1.0/EPSILON;  
   // See pg 1072 of Numerical Recipes typically only need one or two cycles per level
   const int nCycle = 2; 
   // the right boundary can be set to a constant or sinusoindally varing voltage value
   // in the paramters.h file
   double phi_right = phiLX*cos(phi_period*timestep);
   double divdxdy = 1.0 / (DX * DY);

    /* The Numerical Recipes function for calculating phi (Mglin) takes as
     * input a reference to phi with the values of the charge denisty stored
     * and replaces these values with the values for phi
    */
   for (j = 0; j < Y_NO_OF_CELLS; j++)
   {
      for (i = 0; i < X_NO_OF_CELLS; i++)
      {
         phi[j][i] = -coef*rho[j][i]; 
      }
   }

   // Implement the boundary voltages in the x direction
   for (j = 0; j < Y_NO_OF_CELLS; j++)
   {
      phi[j][1] = phi[j][1] - phi0X* divdxdy;
      phi[j][X_NO_OF_CELLS-2] = phi[j][X_NO_OF_CELLS-2] - phi_right* divdxdy;
   }

   // Implement the boundary voltages in the y direction
   for (i = 0; i < X_NO_OF_CELLS; i++)
   {
      phi[1][i] = phi[1][i] - phi0Y* divdxdy;
      phi[Y_NO_OF_CELLS -2][i] = phi[Y_NO_OF_CELLS -2][i] - phiLY* divdxdy;
   }
   
   // Calculate the potential using Mglin, a multi-grid solver from Numerical Recipes
   Mglin(phi,nCycle, XLENGTH);  

   for (j = 0; j < Y_NO_OF_CELLS; j++)
   {
       phi[j][0] = phi0X;
       phi[j][X_NO_OF_CELLS - 1] = phi_right; // - phi_right/coef2;
   }

   for (i = 0; i < X_NO_OF_CELLS; i++)
   {
       phi[0][i] = phi0Y;
       phi[Y_NO_OF_CELLS - 1][i] = phiLY;
   }
}

/* 
 * Implementation Notes: grid_electric_field
 * ------------------------------------------
 * The electric field at each grid point is calculated using the formula in
 * C K Birdsall and A B Langdon, Plasma Physics via Computer Simulation,
 * Chapter 2.
 *
 */
void fields::grid_electric_field(void) 
{
   int i,j;
   double twoDX=2.0*DX; //dx = dy since on square grid

   /* 
      In two dimensions have Ex and Ey given by the potential difference
      over dx (=dy)
   */

   for (j = 1; j< Y_NO_OF_CELLS-1; j++) // Note: stops before boundary (see below)
   {
      for (i =1; i < X_NO_OF_CELLS-1; i++)
      {
         electric_fieldX[j* X_NO_OF_CELLS +  i] = (phi[j][i-1] - phi[j][i+1])/twoDX;
         electric_fieldY[j * X_NO_OF_CELLS + i] = (phi[j-1][i] - phi[j+1][i])/twoDX;
      }
   }

   // At x boundary use Ex[boundary-1] = (Ex[boundary] + Ex[boundary-2])/2
   // and solve for Ex[boundary]
   for(j = 0; j < Y_NO_OF_CELLS; j++)
   {
      electric_fieldX[j * X_NO_OF_CELLS] = 2.0*electric_fieldX[j * X_NO_OF_CELLS + 1] - electric_fieldX[j * X_NO_OF_CELLS + 2];
      electric_fieldX[j * X_NO_OF_CELLS + X_NO_OF_CELLS-1] =
          2.0*electric_fieldX[j * X_NO_OF_CELLS + X_NO_OF_CELLS -2] - electric_fieldX[j * X_NO_OF_CELLS + X_NO_OF_CELLS -3];
      electric_fieldX[j * X_NO_OF_CELLS + X_NO_OF_CELLS] = 0.0;
   }

   // At y boundary use Ey[boundary-1] = (Ey[boundary] + Ey[boundary-2])/2
   // and solve for Ey[boundary]
   for(i = 0; i < X_NO_OF_CELLS; i++)
   {
      electric_fieldY[i] = 2.0*electric_fieldY[X_NO_OF_CELLS + i] - electric_fieldY[2*X_NO_OF_CELLS + i];
      electric_fieldY[(Y_NO_OF_CELLS-1)* X_NO_OF_CELLS + i] =
          2.0*electric_fieldY[(Y_NO_OF_CELLS - 2) * X_NO_OF_CELLS + i] - electric_fieldY[(Y_NO_OF_CELLS - 3) * X_NO_OF_CELLS + i];
      electric_fieldY[Y_NO_OF_CELLS*X_NO_OF_CELLS + i] = 0.0;
   }

}
