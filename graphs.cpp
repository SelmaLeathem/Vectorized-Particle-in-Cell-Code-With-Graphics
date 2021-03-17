/*
 * File: graphs.cpp
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

#include "graphs.h"

/*
* Constructor: graphs
* Usage: graphs graph_object;
*        graphs graph_object(ION_TYPE,SPREAD,VT,QCI);
*/
graphs::graphs(char *type_inp,int spread_in,double vt_in,double qc_in):type_in(type_inp),spread(spread_in),vt(vt_in),qc(qc_in)
{
   int i, j;
   int no_graph_steps = 1;

   // for plotting the net charge density 
   rho_plot = new double*[Y_NO_OF_CELLS]();
      if (!rho_plot){
         printf("Error allocating memory for rho_plot\n");
         exit(1);
      }

      for(i=0;i<(Y_NO_OF_CELLS);i++){
         rho_plot[i] = new double[X_NO_OF_CELLS]();
         if(!rho_plot[i]){
            printf("Error allocating memory for rho_plot\n");
            exit(1);
         }
         for(j=0;j<(X_NO_OF_CELLS);j++){
            rho_plot[i][j] = 0.0;
         }
      }

   // for plotting the electron or ion charge density
   q_plot = new double*[Y_NO_OF_CELLS]();
   if (!q_plot){
      printf("Error allocating memory for q_plot\n");
      exit(1);
   }

   for(i=0;i<(Y_NO_OF_CELLS);i++){
      q_plot[i] = new double[X_NO_OF_CELLS]();
      if(!q_plot[i]){
         printf("Error allocating memory for q_plot\n");
         exit(1);
      }
      for(j=0;j<(X_NO_OF_CELLS);j++){
         q_plot[i][j] = 0.0;
      }
   }

   // for plotting the potential
   phi_plot = new double*[Y_NO_OF_CELLS]();
      if (!phi_plot){
         printf("Error allocating memory for phi_plot\n");
         exit(1);
      }

      for(i=0;i<(Y_NO_OF_CELLS);i++){
         phi_plot[i] = new double[X_NO_OF_CELLS]();
         if(!phi_plot[i]){
            printf("Error allocating memory for phi_plot\n");
            exit(1);
         }
         for(j=0;j<(X_NO_OF_CELLS);j++){
            phi_plot[i][j] = 0.0;
         }
      }

   // for plotting the electric field in the x direction
   Egrid_plotX = new double*[Y_NO_OF_CELLS]();
      if (!Egrid_plotX){
         printf("Error allocating memory for Egrid_plotX\n");
         exit(1);
      }

      for(i=0;i<(Y_NO_OF_CELLS);i++){
         Egrid_plotX[i] = new double[X_NO_OF_CELLS]();
         if(!Egrid_plotX[i]){
            printf("Error allocating memory for Egrid_plotX\n");
            exit(1);
         }
         for(j=0;j<(X_NO_OF_CELLS);j++){
            Egrid_plotX[i][j] = 0.0;
         }
      }

   // for plotting the electric field in the y direction
   Egrid_plotY = new double*[Y_NO_OF_CELLS]();
   if (!Egrid_plotY){
      printf("Error allocating memory for Egrid_plotY\n");
      exit(1);
   }

   for(i=0;i<(Y_NO_OF_CELLS);i++){
      Egrid_plotY[i] = new double[X_NO_OF_CELLS]();
      if(!Egrid_plotY[i]){
         printf("Error allocating memory for Egrid_plotY\n");
         exit(1);
      }
      for(j=0;j<(X_NO_OF_CELLS);j++){
         Egrid_plotY[i][j] = 0.0;
      }
   }

   // for plotting the electron or ion particle x-velocity distribution at location 1
   d_plot1X = new double[NO_BINS+1]();
   if(!d_plot1X){
     printf("Error allocating memory for d_plot1X\n");
     exit(1);
   }

   // for plotting the electron or ion particle x-velocity distribution at location 2
   d_plot2X = new double[NO_BINS+1]();
   if(!d_plot2X){
     printf("Error allocating memory for d_plot2X\n");
     exit(1);
   }

   // for plotting the electron or ion particle y-velocity distribution at location 1
   d_plot1Y = new double[NO_BINS+1]();
   if(!d_plot1Y){
     printf("Error allocating memory for d_plot1Y\n");
     exit(1);
   }

   // for plotting the electron or ion particle y-velocity distribution at location 2
   d_plot2Y = new double[NO_BINS+1]();
   if(!d_plot2Y){
     printf("Error allocating memory for d_plot2Y\n");
     exit(1);
   }

   for(int j=0; j < NO_BINS+1; j++)
   {
      d_plot1X[j] = 0.0;
      d_plot2X[j] = 0.0;
      d_plot1Y[j] = 0.0;
      d_plot2Y[j] = 0.0;
   }

   // for holding temporary data when calculating the above distribution in x
   temp_distX = new cl_double[NO_BINS+1]();
   if(!temp_distX){
     printf("Error allocating memory for temp_distX\n");
     exit(1);
   }

   // for holding temporary data when calculating the above distribution in y
   temp_distY = new cl_double[NO_BINS+1]();
   if(!temp_distY){
     printf("Error allocating memory for temp_distY\n");
     exit(1);
   }

   // for holding temporary data when calculating the above distribution in x
   temp_distX2 = new cl_double[NO_BINS+1]();
   if(!temp_distX2){
     printf("Error allocating memory for temp_distX2\n");
     exit(1);
   }

   // for holding temporary data when calculating the above distribution in y
   temp_distY2 = new cl_double[NO_BINS+1]();
   if(!temp_distY2){
     printf("Error allocating memory for temp_distY2\n");
     exit(1);
   }

   // for plotting the average electron or ion x-velocity per cell along a slice of 
   // the grid in y
   v_plotX = new double[X_NO_OF_CELLS]();
   if(!v_plotX){
     printf("Error allocating memory for v_plotX\n");
     exit(1);
   }

   // for plotting the average electron or ion y-velocity per cell along a slice of 
   // the grid in x
   v_plotY = new double[Y_NO_OF_CELLS]();
   if(!v_plotY){
     printf("Error allocating memory for v_plotY\n");
     exit(1);
   }

   for (int i = 0; i< X_NO_OF_CELLS; i++)
   {
      v_plotX[i] = 0.0;
      v_plotY[i] = 0.0;
   }

   // for summing the average electron or ion x-velocity per cell along a slice of 
   // the grid in x
   v_grid_tempX = new cl_double[X_NO_OF_CELLS]();
   if(!v_grid_tempX){
     printf("Error allocating memory for v_grid_tempX\n");
     exit(1);
   }

   // for summing the average electron or ion y-velocity per cell along a slice of 
   // the grid in x
   v_grid_tempY = new cl_double[X_NO_OF_CELLS]();
   if(!v_grid_tempY){
     printf("Error allocating memory for v_grid_tempY\n");
     exit(1);
   }
   
   // holds the number of particles per cell for the above average x-velocity plots
   count_particlesX = new cl_double[X_NO_OF_CELLS]();
   if(!count_particlesX){
     printf("Error allocating memory for count_particlesX\n");
     exit(1);
   }

   // holds the number of particles per cell for the above average x-velocity plots
   count_particlesY = new cl_double[Y_NO_OF_CELLS]();
   if(!count_particlesY){
     printf("Error allocating memory for count_particlesY\n");
     exit(1);
   }

}

/* Destructor: ~graphs */ 
graphs::~graphs(){

  for(int i=0;i<Y_NO_OF_CELLS;i++)
   {
      delete [] rho_plot[i]; 
      delete [] phi_plot[i]; 
      delete [] Egrid_plotX[i]; 
      delete [] Egrid_plotY[i];
      delete [] q_plot[i];
   }
   delete [] rho_plot;
   delete [] phi_plot;
   delete [] Egrid_plotX;
   delete [] Egrid_plotY;
   delete [] q_plot;

   delete [] v_grid_tempX;
   delete [] v_grid_tempY;

   delete [] d_plot1X;
   delete [] d_plot2X;
   delete [] d_plot1Y;
   delete [] d_plot2Y;
   delete [] v_plotX;
   delete [] v_plotY;

   delete [] count_particlesX;
   delete [] count_particlesY;
   
   delete [] temp_distX;
   delete [] temp_distY;
   delete [] temp_distX2;
   delete [] temp_distY2;
}




/*
* Method: set_type_in
* Usage:  set_type_in("i1");
* ----------------------------
* Define a string that is used to distinguish which species each diagnostic
* text file is reffering to.
*/
void graphs::set_type_in(char* str){
     type_in=str;
} 

/*
* Method: set_spread
* Usage:  set_spread(number_of_cells);
* -------------------------------------
* Set the number of cells around the central one from which to pick
* particle velocities from. These are used to calculate the particle
* velocity distribution.
*/
void graphs::set_spread(int val){
     spread=val;
} 

/*
* Method: set_vt
* Usage:  set_vt(thermal_velocity);
* ----------------------------------
* Set the value of the thermal velocity used to normalize the particle
* velocity distribution.
*/
void graphs::set_vt(double vt_in){
     vt=vt_in;
}

/*
* Method: set_qc
* Usage:  set_qc(charge_of_macro_particle);
* -------------------------------------------
* Use to set the the value of the charge of each macro-particle that is
* multiplied by the particle number density of each cell. To plot the number
* density set this value to 1.0.
*/
void graphs::set_qc(double qc_in){
     qc=qc_in;
}

/* 
* Method: set_dist_location
* Usage:  set_dist_location(cell_numberY, cell_numberX);
* ----------------------------------------   
* Pick the central cell from which particle velocities are picked to
* determine the particle velocity distribution.
*/
void graphs::set_dist_location(int locationX, int locationY){
     dist_locationX=locationX;
     dist_locationY=locationY;
}

/* 
* Method: set_dist_location2
* Usage:  set_dist_location2(cell_number);
* ----------------------------------------   
* Pick a second central cell from which particle velocities are picked to
* determine the particle velocity distribution.
*/
void graphs::set_dist_location2(int locationX, int locationY){
     dist_location2X=locationX;
     dist_location2Y=locationY;
}

/* Sets parameters needed for openCL. The parameters are initialized
* and destroyed in the PIC class and are sent here for use.
*/
void graphs::init_openCL(cl_context context, cl_command_queue cmdQueue,
               cl_kernel kernelVg, cl_kernel kernelVd,cl_mem dCountpX,
               cl_mem dCountpY, cl_mem dv_grid_tempX, cl_mem dv_grid_tempY,
               cl_mem dtemp_distX1, cl_mem dtemp_distX2,
               cl_mem dtemp_distY1, cl_mem dtemp_distY2)
{
      this->context = context;
      this->cmdQueue = cmdQueue;
      this->kernelVg = kernelVg;
      this->kernelVd = kernelVd;
      this->dCountpX = dCountpX;
      this->dCountpY = dCountpY;
      this->dv_grid_tempX = dv_grid_tempX;
      this->dv_grid_tempY = dv_grid_tempY;
      this->dtemp_distX1 = dtemp_distX1;
      this->dtemp_distX2 = dtemp_distX2;
      this->dtemp_distY1 = dtemp_distY1;
      this->dtemp_distY2 = dtemp_distY2;
}

/* Sets parameters needed for openCL for the case of a steadystate
* boundary. In this instance the particle arrays don't change size
* therefore memory need only be reserved during initialization.
*/
void graphs::init_openCL_steady(cl_mem dP)
{
   this->dP = dP;
}

/*
* Method: sum_fields
* Usage:  sum_fields(net_charge_density,electric_potential,electric_field);
* -------------------------------------------------------------------------
* Sums field quantities, including the net_charge density, electric 
* potential and the electric field.
*/
void graphs::sum_fields(double **rho, MatDoub_IO phi,double *EgridX, double *EgridY )
{
   int i,j;

   for (j=0; j < Y_NO_OF_CELLS; j++)
   {
     for (i=0; i < X_NO_OF_CELLS; i++)
     {
       rho_plot[j][i] += rho[j][i];
       phi_plot[j][i] += phi[j][i];
       Egrid_plotX[j][i] += EgridX[j*X_NO_OF_CELLS +i];
       Egrid_plotY[j][i] += EgridY[j * X_NO_OF_CELLS + i];
     }
   }
  
}

/*
* Method: sum_density
* Usage:  sume_density(particle_charge_density);
* -----------------------------------------------
* Sums the charge density of each species.
*/
void graphs::sum_density(cl_double *q)
{
   int i,j;
   
  for (j=0; j < Y_NO_OF_CELLS; j++)
   {
     for (i=0; i < X_NO_OF_CELLS; i++)
     {
        q_plot[j][i] += ((double)q[j * X_NO_OF_CELLS + i] /qc);
     
     }
   }
}


/* 
* Method: sum_av_velocities
* Usage:  sum_av_velocities(particle_vector);
* --------------------------------------------
* Sums the average velocity of each grid cell for the species passed to 
* the arguement. Assumes a square grid. The grid cells are either the
* x directed velocity along x at yPosition in y or the y directed velocity
* directed along the y axis at xPosition in x.
*/
void graphs::sum_av_velocities(vector<particleVec> &particleObj){
   int i,j,k, yPosition, xPosition; //position to take slice of grid at
   double dx = DX, dy = DY;  // cell sizes
   int no_x_cells = X_NO_OF_CELLS, no_y_cells = Y_NO_OF_CELLS;
   int particle_no; // number of particles

   particle_no = particleObj.size();

   // initialize arrays to hold cell average data
   for (j = 0; j < X_NO_OF_CELLS; j++)
   {
      v_grid_tempX[j] = (cl_double)0.0;
      v_grid_tempY[j] = (cl_double)0.0;
      count_particlesX[j] = (cl_double)0.0;
      count_particlesY[j] = (cl_double)0.0;
   }
   
  
   // size of memory areas to reserve on the GPU
   size_t particleSizeVec = particle_no * sizeof(particleVec);
   size_t gridSize = X_NO_OF_CELLS * sizeof(cl_double);

   // 5. allocate the device memory buffers:

   if (STEADYSTATE == "NO")  // then the particle array size and hence dP size varies
   {
       // holds the particle data
       dP = clCreateBuffer(context, CL_MEM_READ_ONLY, particleSizeVec, NULL, &status);
       if (status != CL_SUCCESS)
       {
           fprintf(stderr, "clCreateBuffer failed (1)\n");
           cout << "Status = " << status << endl;
       }
   }
  

   // 6. enqueue the 2 commands to write the data from the host buffers to the device buffers:

   // transfer the local particle data to GPU memory
   status = clEnqueueWriteBuffer( cmdQueue, dP, CL_FALSE, 0, particleSizeVec, &particleObj[0], 0, NULL, NULL );
   if (status != CL_SUCCESS)
   {
         fprintf(stderr, "clEnqueueWriteBuffer failed (1)\n");
         cout << "Status = " << status << endl;
   } 

   // transfer dCountpX
   status = clEnqueueWriteBuffer( cmdQueue, dCountpX, CL_FALSE, 0, gridSize, count_particlesX, 0, NULL, NULL );
   if (status != CL_SUCCESS)
   {
      fprintf(stderr, "clEnqueueWriteBuffer failed (3)\n");
      cout << "Status = " << status << endl;
   }

   // transfer dCountpY
     status = clEnqueueWriteBuffer( cmdQueue, dCountpY, CL_FALSE, 0, gridSize, count_particlesY, 0, NULL, NULL );
   if (status != CL_SUCCESS)
   {
      fprintf(stderr, "clEnqueueWriteBuffer failed (4)\n");
      cout << "Status = " << status << endl;
   }

   // transfer dv_grid_tempX
   status = clEnqueueWriteBuffer( cmdQueue, dv_grid_tempX, CL_FALSE, 0, gridSize, v_grid_tempX, 0, NULL, NULL );
   if (status != CL_SUCCESS)
   {
      fprintf(stderr, "clEnqueueWriteBuffer failed (5)\n");
      cout << "Status = " << status << endl;
   }

   // transfer dv_grid_tempY
   status = clEnqueueWriteBuffer( cmdQueue, dv_grid_tempY, CL_FALSE, 0, gridSize, v_grid_tempY, 0, NULL, NULL );
   if (status != CL_SUCCESS)
   {
      fprintf(stderr, "clEnqueueWriteBuffer failed (6)\n");
      cout << "Status = " << status << endl;
   }

   // ensure the above steps have completed
    Wait(cmdQueue);
		

    // 10. setup the arguments to the kernel object:

   // particle data array
   status = clSetKernelArg(kernelVg, 0, sizeof(cl_mem), &dP);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (0)\n");
       cout << "Status = " << status << endl;
   }

   // number of particles per cell where y = length/2
   status = clSetKernelArg(kernelVg, 1, sizeof(cl_mem), &dCountpX);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (1)\n");
       cout << "Status = " << status << endl;
   }

   // number of particles per cell where x = length/2
   status = clSetKernelArg(kernelVg, 2, sizeof(cl_mem), &dCountpY);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (2)\n");
       cout << "Status = " << status << endl;
   }

   //holds temp values for summing average velocities
   status = clSetKernelArg(kernelVg, 3, sizeof(cl_mem), &dv_grid_tempX);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (3)\n");
       cout << "Status = " << status << endl;
   }

   //holds temp values for summing average velocities
   status = clSetKernelArg(kernelVg, 4, sizeof(cl_mem), &dv_grid_tempY);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (4)\n");
       cout << "Status = " << status << endl;
   }

   // DX cell size in x
   status = clSetKernelArg(kernelVg, 5, sizeof(cl_double), &dx);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (5)\n");
       cout << "Status = " << status << endl;
   }

   // DY cell size in y
   status = clSetKernelArg(kernelVg, 6, sizeof(cl_double), &dy);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (6)\n");
       cout << "Status = " << status << endl;
   }

   // X_NO_CELLS
   status = clSetKernelArg(kernelVg, 7, sizeof(cl_int), &no_x_cells);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (7)\n");
       cout << "Status = " << status << endl;
   }

   // Y_NO_CELLS
   status = clSetKernelArg(kernelVg, 8, sizeof(cl_int), &no_y_cells);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (8)\n");
       cout << "Status = " << status << endl;
   }

   // 11. enqueue the kernel object for execution:

   size_t globalWorkSize[3] = { particle_no , 1, 1 };
   size_t localWorkSize[3] = { LOCAL_SIZE,   1, 1 };

   Wait(cmdQueue);

   status = clEnqueueNDRangeKernel(cmdQueue, kernelVg, 1, NULL, globalWorkSize, NULL, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueNDRangeKernel failed: %d\n", status);
       cout << "Status = " << status << endl;
   }

   Wait(cmdQueue); 

    // 12. read the results buffer back from the device to the host:

   status = clEnqueueReadBuffer(cmdQueue, dCountpX, CL_TRUE, 0, gridSize, count_particlesX, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueReadBuffer failed (1)\n");
       cout << "Status = " << status << endl;
   }

   status = clEnqueueReadBuffer(cmdQueue, dCountpY, CL_TRUE, 0, gridSize, count_particlesY, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueReadBuffer failed (2)\n");
       cout << "Status = " << status << endl;
   }

   status = clEnqueueReadBuffer(cmdQueue, dv_grid_tempX, CL_TRUE, 0, gridSize, v_grid_tempX, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueReadBuffer failed (3)\n");
       cout << "Status = " << status << endl;
   }

   status = clEnqueueReadBuffer(cmdQueue, dv_grid_tempY, CL_TRUE, 0, gridSize, v_grid_tempY, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueReadBuffer failed (4)\n");
       cout << "Status = " << status << endl;
   }

   // At the end of each function call flush the queue
   // reference: https://github.com/conanwu777/particle_system/blob/master/main.cpp
   clFlush(cmdQueue);

   if (STEADYSTATE == "NO")  // release reserved memory
   {
       clReleaseMemObject(dP);
   }
 

   // divide the total sum of the particle velocities for each cell by the a count  
   // of the number of particles for that cell
   for(j=0; j< X_NO_OF_CELLS ; j++){
      if(count_particlesX[j] > 0.5)
      {
         v_plotX[j] += (double)(v_grid_tempX[j]/count_particlesX[j]);  
      }
      if(count_particlesY[j] > 0.5)
      {
         v_plotY[j] += (double)(v_grid_tempY[j]/count_particlesY[j]);
      }
   }
} 

/*
* Methods: sum_av_dist
* Usage:   sum_av_dist(particle_vector), sum_av_dist2(particle_vector);
* -----------------------------------------------------------------------
* Sums the particle velocity distribution at cell locations defined by the
* dist_location and dist_location2 variables. Particle velocities are 
* included from several cells around each location indicated by the 
* spread variable.
*/ 
void graphs::sum_av_dist(vector<particleVec> &particleObj){
   int i,j,k, particle_no; // number of particles
   double dx = DX, dy = DY; // cell sizes
   int no_x_cells = X_NO_OF_CELLS, no_y_cells = Y_NO_OF_CELLS;
   int no_bins = NO_BINS; //number of velocity bins 
   double dspread = spread; //number of cells about the positions to take velocities
   // positions on the grid to calculated distributions
   double dist_locationX1 = dist_locationX; // for vx
   double dist_locationX2 = dist_location2X; // for vx
   double dist_locationY1 = dist_locationY;  // for vy
   double dist_locationY2 = dist_location2Y; // for vy
   double dvt = vt; // thermal velocity
   double divide_bins_by = DIVIDE_BINS_BY; // for scaling 

   particle_no = particleObj.size();

    // initialize arrays to hold cell average data
   for(i=0; i< NO_BINS + 1; i++)
   {
      temp_distX[i] = (cl_double)0.0;
      temp_distY[i] = (cl_double)0.0;
      temp_distX2[i] = (cl_double)0.0;
      temp_distY2[i] = (cl_double)0.0;
   }

   // size of memory areas to reserve on the GPU
   size_t particleSizeVec = particle_no * sizeof(particleVec);
   size_t binSize = (NO_BINS + 1) * sizeof(cl_double);

   // 5. allocate the device memory buffers:

   if (STEADYSTATE == "NO")
   {
       // holds the particle data
       dP = clCreateBuffer(context, CL_MEM_READ_ONLY, particleSizeVec, NULL, &status);
       if (status != CL_SUCCESS)
       {
           fprintf(stderr, "clCreateBuffer failed (1)\n");
           cout << "Status = " << status << endl;
       }
   }

   // 6. enqueue the 2 commands to write the data from the host buffers to the device buffers:

   // transfer the local particle data to GPU memory
   status = clEnqueueWriteBuffer( cmdQueue, dP, CL_FALSE, 0, particleSizeVec, &particleObj[0], 0, NULL, NULL );
   if (status != CL_SUCCESS)
   {
         fprintf(stderr, "clEnqueueWriteBuffer failed (1)\n");
         cout << "Status = " << status << endl;
   } 

   // transfer dtemp_distX1
   status = clEnqueueWriteBuffer( cmdQueue, dtemp_distX1, CL_FALSE, 0, binSize, temp_distX, 0, NULL, NULL );
   if (status != CL_SUCCESS)
   {
      fprintf(stderr, "clEnqueueWriteBuffer failed (3)\n");
      cout << "Status = " << status << endl;
   }   

   // transfer dtemp_distX2
   status = clEnqueueWriteBuffer( cmdQueue, dtemp_distX2, CL_FALSE, 0, binSize, temp_distX2, 0, NULL, NULL );
   if (status != CL_SUCCESS)
   {
      fprintf(stderr, "clEnqueueWriteBuffer failed (3)\n");
      cout << "Status = " << status << endl;
   } 

   //transfer dtemp_distY1
   status = clEnqueueWriteBuffer( cmdQueue, dtemp_distY1, CL_FALSE, 0, binSize, temp_distY, 0, NULL, NULL );
   if (status != CL_SUCCESS)
   {
      fprintf(stderr, "clEnqueueWriteBuffer failed (3)\n");
      cout << "Status = " << status << endl;
   } 

   //transfer dtemp_distY2
   status = clEnqueueWriteBuffer( cmdQueue, dtemp_distY2, CL_FALSE, 0, binSize, temp_distY2, 0, NULL, NULL );
   if (status != CL_SUCCESS)
   {
      fprintf(stderr, "clEnqueueWriteBuffer failed (3)\n");
      cout << "Status = " << status << endl;
   } 

   // ensure the above steps have completed
    Wait(cmdQueue);
		

    // 10. setup the arguments to the kernel object:

   // particle data array
   status = clSetKernelArg(kernelVd, 0, sizeof(cl_mem), &dP);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (0)\n");
       cout << "Status = " << status << endl;
   }

   //holds temp values for summing average vx distribution
   status = clSetKernelArg(kernelVd, 1, sizeof(cl_mem), &dtemp_distX1);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (1)\n");
       cout << "Status = " << status << endl;
   }

   //holds temp values for summing average vx distribution
   status = clSetKernelArg(kernelVd, 2, sizeof(cl_mem), &dtemp_distX2);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (2)\n");
       cout << "Status = " << status << endl;
   }

   //holds temp values for summing average vy distribution
   status = clSetKernelArg(kernelVd, 3, sizeof(cl_mem), &dtemp_distY1);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (3)\n");
       cout << "Status = " << status << endl;
   }

   //holds temp values for summing average vy distribution
   status = clSetKernelArg(kernelVd, 4, sizeof(cl_mem), &dtemp_distY2);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (4)\n");
       cout << "Status = " << status << endl;
   }

    // DX cell size in x
   status = clSetKernelArg(kernelVd, 5, sizeof(cl_double), &dx);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (5)\n");
       cout << "Status = " << status << endl;
   }

   // DY cell size in y
   status = clSetKernelArg(kernelVd, 6, sizeof(cl_double), &dy);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (6)\n");
       cout << "Status = " << status << endl;
   }

   // X_NO_CELLS
   status = clSetKernelArg(kernelVd, 7, sizeof(cl_int), &no_x_cells);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (7)\n");
       cout << "Status = " << status << endl;
   }

   // Y_NO_CELLS
   status = clSetKernelArg(kernelVd, 8, sizeof(cl_int), &no_y_cells);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (8)\n");
       cout << "Status = " << status << endl;
   }

   //number of velocity bins
   status = clSetKernelArg(kernelVd, 9, sizeof(cl_int), &no_bins);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (9)\n");
       cout << "Status = " << status << endl;
   }

   //number of cells about the positions to take velocities
   status = clSetKernelArg(kernelVd, 10, sizeof(cl_double), &spread);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (10)\n");
       cout << "Status = " << status << endl;
   }

   // positions on the grid to calculated distributions

   status = clSetKernelArg(kernelVd, 11, sizeof(cl_double), &dist_locationX1);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (11)\n");
       cout << "Status = " << status << endl;
   }

   status = clSetKernelArg(kernelVd, 12, sizeof(cl_double), &dist_locationX2);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (12)\n");
       cout << "Status = " << status << endl;
   }

      status = clSetKernelArg(kernelVd, 13, sizeof(cl_double), &dist_locationY1);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (13)\n");
       cout << "Status = " << status << endl;
   }

      status = clSetKernelArg(kernelVd, 14, sizeof(cl_double), &dist_locationY2);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (14)\n");
       cout << "Status = " << status << endl;
   }

   // thermal velocity
   status = clSetKernelArg(kernelVd, 15, sizeof(cl_double), &dvt);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (15)\n");
       cout << "Status = " << status << endl;
   }

   // for scaling
   status = clSetKernelArg(kernelVd, 16, sizeof(cl_double), &divide_bins_by);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (16)\n");
       cout << "Status = " << status << endl;
   }

   // 11. enqueue the kernel object for execution:

   size_t globalWorkSize[3] = { particle_no , 1, 1 };
   size_t localWorkSize[3] = { LOCAL_SIZE,   1, 1 };

   Wait(cmdQueue);

   status = clEnqueueNDRangeKernel(cmdQueue, kernelVd, 1, NULL, globalWorkSize, NULL, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueNDRangeKernel failed: Vd %d\n", status);
       cout << "Status = " << status << endl;
   }

   Wait(cmdQueue); 

   // 12. read the results buffer back from the device to the host:

   status = clEnqueueReadBuffer(cmdQueue, dtemp_distX1, CL_TRUE, 0, binSize, temp_distX, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueReadBuffer failed (1)\n");
       cout << "Status = " << status << endl;
   }  

      status = clEnqueueReadBuffer(cmdQueue, dtemp_distX2, CL_TRUE, 0, binSize, temp_distX2, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueReadBuffer failed (2)\n");
       cout << "Status = " << status << endl;
   } 

      status = clEnqueueReadBuffer(cmdQueue, dtemp_distY1, CL_TRUE, 0, binSize, temp_distY, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueReadBuffer failed (3)\n");
       cout << "Status = " << status << endl;
   } 

      status = clEnqueueReadBuffer(cmdQueue, dtemp_distY2, CL_TRUE, 0, binSize, temp_distY2, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueReadBuffer failed (4)\n");
       cout << "Status = " << status << endl;
   } 

   // At the end of each function call flush the queue
   // reference: https://github.com/conanwu777/particle_system/blob/master/main.cpp
   clFlush(cmdQueue);

   if (STEADYSTATE == "NO") // release the gpu memory
   {
       clReleaseMemObject(dP);
   }

   // sum the distributions for each iteration to later divide by 
   // the number of iterations to get an average 
   for(j=0; j< NO_BINS; j++)
   {
      d_plot1X[j] += (double)temp_distX[j];
      d_plot2X[j] += (double)temp_distX2[j];
      d_plot1Y[j] += (double)temp_distY[j];
      d_plot2Y[j] += (double)temp_distY2[j];
   }
 
}

/* 
* Set the number of steps over which to take averages for graph
* output.
*/
void graphs::set_no_graph_steps(int numSteps)
{
    no_graph_steps = numSteps;
}

/*
* The following function takes an averages of the indicated quantity by
* dividing it over the number of steps each item was summed up over. 
* In addition it calculates the maximum average value of each field
* quantity.
*/
void graphs::average_fields(void){
  int i, j;

  // initialize maximum values
  phiMax = phi_plot[0][0] / ((double)no_graph_steps);
  ExMax = Egrid_plotX[0][0] / ((double)no_graph_steps);
  EyMax = Egrid_plotY[0][0] / ((double)no_graph_steps);

   for (j=0; j < Y_NO_OF_CELLS; j++)
   {
     for (i=0; i < X_NO_OF_CELLS; i++)
     {
       rho_plot[j][i] /= ((double)no_graph_steps);

       phi_plot[j][i] /= ((double)no_graph_steps);

       if (phiMax < phi_plot[j][i])  
           phiMax = phi_plot[j][i];  // get maximum

       Egrid_plotX[j][i] /= ((double)no_graph_steps);

       if (ExMax < Egrid_plotX[j][i]) // get maximum
           ExMax = Egrid_plotX[j][i];

       Egrid_plotY[j][i] /= ((double)no_graph_steps);

       if (EyMax < Egrid_plotY[j][i])  // get maximum
           EyMax = Egrid_plotY[j][i];
     }
   }
}

/*
* Used to retrieve averaged graph class data for output plots rendered  
* in the graphics output functions.
*/

// electrostatic potential
double** graphs::get_phiPlot()
{
    return phi_plot;
}

// x component of the electric field
double** graphs::getExPlot(void)
{
    return Egrid_plotX;
}

// y component of the electric field
double** graphs::getEyPlot(void)
{
    return Egrid_plotY;
}

// max value of the electrostatic potential
double graphs::getPhiMax(void)
{
    return phiMax;
}

// max value of the x component of the electric field
double graphs::getExMax(void)
{
    return ExMax;
}

// max value of the y component of the electric field
double graphs::getEyMax(void)
{
    return EyMax;
}

// density
double** graphs::getQPlot(void)
{
    return q_plot;
}

// maximum value of the density
double graphs::getQMax(void)
{
    return qMax;
}

/*
* The following function takes an average of the indicated quantity by
* dividing it over the number of steps each item was summed up over. 
* In addition it calculates the maximum average value of each field
* quantity.
*/
void graphs::average_density(void){  
  int i,j;
  // initialize maximum value
  qMax = q_plot[1][1] / ((double)no_graph_steps);
  
   for (j=0; j < Y_NO_OF_CELLS; j++)
   {
     for (i=0; i < X_NO_OF_CELLS; i++)
     {
        q_plot[j][i] /= ((double)no_graph_steps);

        // get maximum
        if (qMax < q_plot[j][i])
            qMax = q_plot[j][i];
     }
   }

}

/*
* The following function takes an average of the indicated quantity by
* dividing it over the number of steps each item was summed up over. 
*/
void graphs::average_velocities(void){
  int j;

  for(j=0;j < X_NO_OF_CELLS; j++){
     *(v_plotX+j) = *(v_plotX+j)/(VTE*((double)no_graph_steps));
     *(v_plotY+j) = *(v_plotY+j)/(VTE*((double)no_graph_steps));
  }
  
}

/*
* The following function takes an average of the indicated quantity by
* dividing it over the number of steps each item was summed up over. 
*/
void graphs::average_dist1(void){
  int j;

  for(j=0; j< NO_BINS; j++)
  {
     *(d_plot1X + j) = *(d_plot1X + j)/((double)no_graph_steps);
     *(d_plot1Y + j) = *(d_plot1Y + j)/((double)no_graph_steps);
  }
    
}

/*
* The following function takes an average of the indicated quantity by
* dividing it over the number of steps each item was summed up over. 
*/
void graphs::average_dist2(void){
  int j;

  for(j=0; j< NO_BINS; j++)
  {
     *(d_plot2X + j) = *(d_plot2X + j)/((double)no_graph_steps);
     *(d_plot2Y + j) = *(d_plot2Y + j)/((double)no_graph_steps);
  }

}

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

// rho, phi, Ex, Ey slice 
void graphs::write_fields_to_file(void){
  int i, yPosition, xPosition;
  FILE *filePtr;
  FILE *filePtr2;

  yPosition = Y_NO_OF_CELLS/2;
  xPosition = X_NO_OF_CELLS/2;

  filePtr = fopen("rho_phi_EgridX.txt", "w");
  filePtr2 = fopen("rho_phi_EgridY.txt", "w");
  
  // just plot along the center line in y
  for(i=0; i< X_NO_OF_CELLS; i++){
     fprintf(filePtr, "%e %e %e %e\n",(double)i, rho_plot[yPosition][i], phi_plot[yPosition][i], Egrid_plotX[yPosition][i]);
     fprintf(filePtr2, "%e %e %e %e\n",(double)i, rho_plot[i][xPosition], phi_plot[i][xPosition], Egrid_plotY[i][xPosition]);
  }
     
  fclose(filePtr);
  fclose(filePtr2);
}

// Rho on the two-dimensional grid
void graphs::write_rho_to_file(void)
{
   int i,j;
   FILE *filePtr;

   filePtr = fopen("rho2d.txt", "w");

   for (j=0; j < Y_NO_OF_CELLS; j++)
   {

     for (i=0; i < X_NO_OF_CELLS; i++)
     {
       fprintf(filePtr, "%e ",rho_plot[j][i]);
     }

     fprintf(filePtr, "\n");
   }

    fclose(filePtr);
}

// Phi on the two-dimensional grid
void graphs::write_phi_to_file(void)
{
   int i,j;
   FILE *filePtr;

   filePtr = fopen("phi2d.txt", "w");

   for (j=0; j < Y_NO_OF_CELLS; j++)
   {
      
     for (i=0; i < X_NO_OF_CELLS; i++)
     {
       fprintf(filePtr, "%e ",phi_plot[j][i]);
     }

     fprintf(filePtr, "\n");
   }

    fclose(filePtr);
}

// EgridX on the two-dimensional grid
void graphs::write_EgridX_to_file(void)
{
   int i,j;
   FILE *filePtr;

   filePtr = fopen("egrid2dX.txt", "w");

   for (j=0; j < Y_NO_OF_CELLS; j++)
   {
      
     for (i=0; i < X_NO_OF_CELLS; i++)
     {
       fprintf(filePtr, "%e ",Egrid_plotX[j][i]);
     }

     fprintf(filePtr, "\n");
   }

    fclose(filePtr);
}

// EgridY on the two-dimensional grid
void graphs::write_EgridY_to_file(void)
{
   int i,j;
   FILE *filePtr;

   filePtr = fopen("egrid2dY.txt", "w");

   for (j=0; j < Y_NO_OF_CELLS; j++)
   {
      
     for (i=0; i < X_NO_OF_CELLS; i++)
     {
       fprintf(filePtr, "%e ",Egrid_plotY[j][i]);
     }

     fprintf(filePtr, "\n");
   }

    fclose(filePtr);
}

// Density on the two-dimensional grid
void graphs::write_density_to_file(void){
  int j, i;
  FILE *filePtr;

  char filename[80] = "density2d_"; 

  strcat(filename, type_in);
  strcat(filename, ".txt");

  filePtr = fopen(filename, "w");

   for (j=0; j < Y_NO_OF_CELLS; j++)
   {
      
     for (i=0; i < X_NO_OF_CELLS; i++)
     {
       fprintf(filePtr, "%e ",q_plot[j][i]);
     }

     fprintf(filePtr, "\n");
   }

  fclose(filePtr);
}

// average velocity per cell in a slice
void graphs::write_v_to_file(void){
  int j;
  FILE *filePtr;

  char filename[80] = "av_velocity_";
  strcat(filename, type_in);
  strcat(filename, ".txt");

  filePtr = fopen(filename, "w");

  for(j=0; j< X_NO_OF_CELLS; j++)
    fprintf(filePtr,"%e %e %e \n",(double)j,*(v_plotX+j),*(v_plotY+j) );

  fclose(filePtr);
}

// average velocity distribution at location 1
void graphs::write_dist1_to_file(void){
  int j;
  FILE *filePtrX;
  FILE *filePtrY;
  char filenameX[80] = "first_v_distX_";
  char filenameY[80] = "first_v_distY_";

  strcat(filenameX, type_in);
  strcat(filenameX, ".txt");
  filePtrX = fopen(filenameX, "w");

  strcat(filenameY, type_in);
  strcat(filenameY, ".txt");
  filePtrY = fopen(filenameY, "w");
   
  for(j=0; j< NO_BINS; j++)
  {
      fprintf(filePtrX,"%e %e \n",(double)(j-(NO_BINS-1)/2)/((double)DIVIDE_BINS_BY), *(d_plot1X+j));
      fprintf(filePtrY,"%e %e \n",(double)(j-(NO_BINS-1)/2)/((double)DIVIDE_BINS_BY), *(d_plot1Y+j));
  }
   
  fclose(filePtrX);
  fclose(filePtrY);
}

// average velocity distribution at location 2
void graphs::write_dist2_to_file(void){
  int j;
   FILE *filePtrX;
  FILE *filePtrY;
  char filenameX[80] = "second_v_distX_";
  char filenameY[80] = "second_v_distY_";

  strcat(filenameX, type_in);
  strcat(filenameX, ".txt");
  filePtrX = fopen(filenameX, "w");

  strcat(filenameY, type_in);
  strcat(filenameY, ".txt");
  filePtrY = fopen(filenameY, "w");
   
  for(j=0; j< NO_BINS; j++)
  {
      fprintf(filePtrX,"%e %e \n",(double)(j-(NO_BINS-1)/2)/((double)DIVIDE_BINS_BY), *(d_plot2X+j));
      fprintf(filePtrY,"%e %e \n",(double)(j-(NO_BINS-1)/2)/((double)DIVIDE_BINS_BY), *(d_plot2Y+j));
  }

  fclose(filePtrX);
  fclose(filePtrY);
}

