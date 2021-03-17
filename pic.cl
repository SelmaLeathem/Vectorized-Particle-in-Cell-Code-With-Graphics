/************************************************************************************************
 * Date:        2/25/2021
 * Description: Function definitions for a PIC simulation. Functions replaced in the
                code include:
				     collect_charge
					 get_particle_electric_field
					 move_particles
					 sum_av_velocities
					 sum_av_dist
 ***********************************************************************************************/
typedef struct particleVec
{
	double x;        /*position along x axis*/
	double y;        /*position along y axis*/
	double vx;        /*velocity along x axis*/
	double vy;       /*y axis velocity*/
	double vz;       /*z axis velocity*/
	double weight;   /*wieght of particle*/
} particle;

 /* replaces collect charge function */
kernel
void
collectCharge(global const particle* dP, global double* dQ, const double qc_dxdy, const double dx, const double dy,
	const int no_x_cells, const int no_y_cells)
{
	int gid = get_global_id(0);
	int jp, kp;  // the nearest grid number to the particle's current position
	double wj, wk, a1, a2, a3, a4; //used to wieght to grid

	jp = (int)(dP[gid].x / dx);
	kp = (int)(dP[gid].y / dy);

	// check the particle is within the boundary
	if (jp < no_x_cells - 1 && dP[gid].x > 0.0 && kp < no_y_cells - 1 && dP[gid].y > 0.0) {
		wj = dP[gid].x / dx - (double)jp;
		wk = dP[gid].y / dy - (double)kp;

		//Birdsall and Langdon pg 308
		a4 = (1.0 - wj) * (1.0 - wk);
		a3 = wj * (1.0 - wk);
		a2 = wk * (1.0 - wj);
		a1 = wj * wk;

		// charge density at the grid points
		dQ[kp * no_x_cells + jp] += qc_dxdy * a4 * dP[gid].weight;
		dQ[kp * no_x_cells + jp + 1] += qc_dxdy * a3 * dP[gid].weight;
		dQ[(kp + 1) * no_x_cells + jp] += qc_dxdy * a2 * dP[gid].weight;
		dQ[(kp + 1) * no_x_cells + jp + 1] += qc_dxdy * a1 * dP[gid].weight;

	}
}


 /* replaces get particle electric field function */
kernel
void
get_particle_electric_field(global const particle* dP, global double* dExp, global double* dEyp,
	global const double* dEx, global const double* dEy, const double dx,
	const double dy, const int no_x_cells, const int no_y_cells)
{
	int gid = get_global_id(0);
	int jp, kp; // the nearest grid number to the particle's current position
	double wj, wk, a1, a2, a3, a4; //used to wieght to grid

	jp = (int)(dP[gid].x / dx);
	kp = (int)(dP[gid].y / dy);

	//initialize fields to zero
	dExp[gid] = 0.0;
	dEyp[gid] = 0.0;

	// check the particle is within the boundary
	if (jp < no_x_cells && dP[gid].x > 0.0 && kp < no_y_cells && dP[gid].y > 0.0) {
		wj = dP[gid].x / dx - (double)jp;
		wk = dP[gid].y / dy - (double)kp;

		//Birdsall and Langdon pg 308
		a4 = (1.0 - wj) * (1.0 - wk);
		a3 = wj * (1.0 - wk);
		a2 = wk * (1.0 - wj);
		a1 = wj * wk;

		// https://people.nscl.msu.edu/~lund/uspas/sbp_2018/lec_intro/04.pic.pdf
		dExp[gid] = a4 * dEx[kp * no_x_cells + jp] + a3 * dEx[kp * no_x_cells + jp + 1] +
			a2 * dEx[(kp + 1) * no_x_cells + jp] + a1 * dEx[(kp + 1) * no_x_cells + jp + 1];
		dEyp[gid] = a4 * dEy[kp * no_x_cells + jp] + a3 * dEy[kp * no_x_cells + jp + 1] +
			a2 * dEy[(kp + 1) * no_x_cells + jp] + a1 * dEy[(kp + 1) * no_x_cells + jp + 1];

	}

}

/* replaces move_particles function */
kernel
void
move_particles(global  particle* dP,  global const double* dExp, global const double* dEyp,
	                         const double dx, const double dy, const int no_x_cells, const int no_y_cells,
							 const double dvdt_coef, const double omega_dt, const double cos_angle, const double sin_angle,
							 const double dt, const double cos_omega, const double sin_omega, const double magnetic_field
	                        )
{
	int gid = get_global_id(0);
	int jp, kp; // the nearest grid number to the particle's current position
	//particle velocities perpendicular and parallel to B
	double vperp,vpar,vperp_after; 

	jp = (int)(dP[gid].x / dx); 
	kp = (int)(dP[gid].y / dy);

	// check the particle is within the boundary
	if (jp < no_x_cells && dP[gid].x > 0.0 && kp < no_y_cells && dP[gid].y > 0.0) {

		 /***  - 1/2 accelerate due to the electric field- ***/

        dP[gid].vx += dvdt_coef* dExp[gid];
        dP[gid].vy += dvdt_coef* dEyp[gid];

		  /*** -rotate around magnetic field- ***/

        if(magnetic_field > 0.0 ){

          /*perpendicular and parallel velocities*/

          vpar = cos_angle*dP[gid].vx  + sin_angle*dP[gid].vy;
          vperp = -sin_angle*dP[gid].vx + cos_angle*dP[gid].vy;

         /* rotate perpendicular components */
      
          vperp_after = vperp*cos_omega +sin_omega*dP[gid].vz;
          dP[gid].vz = -vperp*sin_omega + cos_omega*dP[gid].vz;

         /* get vx and vy from vperp_after and vpar */
   
          dP[gid].vx = cos_angle*vpar - sin_angle*vperp_after;
          dP[gid].vy = sin_angle*vpar + cos_angle*vperp_after;
        }

		/***  - 1/2 accelerate due to the electric field - ***/

		dP[gid].vx += dvdt_coef * dExp[gid];
		dP[gid].vy += dvdt_coef * dEyp[gid];

		 /* Move particle along the x-axis according to standard equations of motion*/
        dP[gid].x  += dt*dP[gid].vx;
        dP[gid].y  += dt*dP[gid].vy;
	 
	}
}

/* replaces sum_av_velocities function */
kernel
void
sum_av_velocities(global const particle* dP, global double* dCountpX,
	global double* dCountpY, global double* dv_grid_tempX, global double* dv_grid_tempY,
	const double dx, const double dy, const int no_x_cells, const int no_y_cells
)
{
	int gid = get_global_id(0);
	int jp, kp;	// the nearest grid number to the particle's current position
	int xPosition = no_x_cells / 2; // position in x where to take y slice
	int yPosition = no_y_cells / 2; // position in y where to take x slice

	jp = (int)(dP[gid].x / dx); 
	kp = (int)(dP[gid].y / dy);


	// if x is located in a cell along x at yPosition in y then 
	// use it in the average velocity calculation

	if (jp < no_x_cells && dP[gid].x > 0.0 && kp == yPosition) {

		dv_grid_tempX[jp] += dP[gid].vx;
		dCountpX[jp] += 1.0;

	}

	// if  y is located in a cell along y at yPosition in x then 
	// use it in the average velocity calculation

	if (jp == xPosition && kp < no_y_cells && dP[gid].y > 0.0) {
		dv_grid_tempY[kp] += dP[gid].vy;
		dCountpY[kp] += 1.0;

	}
}

/* replaces sum_av_dist functions */
kernel
void
sum_av_dist(global const particle* dP, global double* dtemp_distX1,
	global double* dtemp_distX2, global double* dtemp_distY1, global double* dtemp_distY2,
	const double dx, const double dy, const int no_x_cells, const int no_y_cells,
	const int no_bins, const double spread, const double dist_locationX1,
	const double dist_locationX2, const double dist_locationY1, const double dist_locationY2,
	const double vt, const double divide_bins_by
)
{
	int gid = get_global_id(0);
	int jp, kp; // the nearest grid number to the particle's current position
	int vbin; // the selected velocity bin
	int xPosition = no_y_cells / 2; // position in x to take distribution
	int yPosition = no_y_cells / 2; // position in y to take distribution

	jp = (int)(dP[gid].x / dx); 
	kp = (int)(dP[gid].y / dy);


	/*
		If x or y is located within spread distance from the above positions then
		plop that particle in the appropriate velocity bin:
	*/


	if (jp >= (dist_locationX1 - spread) && dP[gid].x > 0.0 && jp <= (dist_locationX1 + spread)
		&& jp < no_x_cells && kp == yPosition) {

		vbin = (int)((dP[gid].vx / vt) * divide_bins_by);

		if (abs(vbin) <= (no_bins - 1) / 2)
			dtemp_distX1[vbin + (no_bins - 1) / 2] += 1.0;  //shift about zero
	}

	if (jp >= (dist_locationX2 - spread) && dP[gid].x > 0.0 && jp <= (dist_locationX2 + spread)
		&& jp < no_x_cells && kp == yPosition) {

		vbin = (int)((dP[gid].vx / vt) * divide_bins_by);

		if (abs(vbin) <= (no_bins - 1) / 2)
			dtemp_distX2[vbin + (no_bins - 1) / 2] += 1.0; //shift about zero
	}

	if (kp >= (dist_locationY1 - spread) && dP[gid].y > 0.0 && kp <= (dist_locationY1 + spread)
		&& kp < no_y_cells && jp == xPosition) {

		vbin = (int)((dP[gid].vy / vt) * divide_bins_by);

		if (abs(vbin) <= (no_bins - 1) / 2)
			dtemp_distY1[vbin + (no_bins - 1) / 2] += 1.0; //shift about zero
	}

	if (kp >= (dist_locationY2 - spread) && dP[gid].y > 0.0 && kp <= (dist_locationY2 + spread)
		&& kp < no_y_cells && jp == xPosition) {

		vbin = (int)((dP[gid].vy / vt) * divide_bins_by);

		if (abs(vbin) <= (no_bins - 1) / 2)
			dtemp_distY2[vbin + (no_bins - 1) / 2] += 1.0; //shift about zero
	}

}







