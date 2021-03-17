/*
 * File: particle.cpp
 * --------------------
 * 
 * 
 * This file implements the particle class which stores particle vectors and
 * associated functions that operate on particles, such as acceleration by
 * the grid electric field, and the collection of the total particle charge    
 * in each cell.
 *
 * 4/15
 * 
 * Date: 2/2/2021
 * 
 * Using the formulae from "Plasma Physics Via Computer Simulation" by
 * Birdsall and Langdon in chapter 14 the code was adapted to account for
 * a two-dimensional grid. This required adding an additional function that
 * weights the x and y components of the electric field at the grid to the
 * individual particles.
 * 
 * Date: 2/8/2021
 * 
 * Vectorizing onto the GPU was applied to the following functions:
 * 
 *                            collect_charge
 *                            get_particle_electric_field
 *                            move_particles
 *
 */

#include "particle.h"


/*
* Constructor: particle
* Usage: particle;
*        particle(dt,vt,q_over_m,qc,mass,dxp);
*
*/
particle::particle(double dt_p, double vt_p, double q_over_m_p, double qc_p,double mass_p,double dxp_p, double dyp_p):dt(dt_p),
                     vt(vt_p), q_over_m(q_over_m_p), qc(qc_p),mass(mass_p),dxp(dxp_p),dyp(dyp_p)
{
   count_reduce = 0;
   int i,j;

   //charge density on the grid
   q = new cl_double[X_NO_OF_CELLS*Y_NO_OF_CELLS]();
   if (!q) {
       printf("Error allocating memory for q\n");
       exit(1);
   }

  // the upper velocity limit of a reduction group  
   v_limit = new double[NO_GROUPS]();
   if(!v_limit){
     printf("Error allocating memory for v_limit\n");
     exit(1);
   }

  // the upper size limit of reduction group
   group_limit = new int[NO_GROUPS]();
   if(!group_limit){
     printf("Error allocating memory for group_limit\n");
     exit(1);
   }

  // the size of a velocity group after a reduction
   new_group_size = new int[NO_GROUPS]();
   if(!new_group_size){
     printf("Error allocating memory for new_group_size\n");
     exit(1);
   }

  // the new wieghting of the group
   group_weight = new double[NO_GROUPS]();
   if(!group_weight){
     printf("Error allocating memory for group_weight\n");
     exit(1);
   }

  // used to identify group membership
   j_group = new int[NO_GROUPS]();
   if(!j_group){
     printf("Error allocating memory for j_group\n");
     exit(1);
   }

  // used to identify group sizes
   Ngroup = new int[NO_GROUPS]();
   if(!Ngroup){
     printf("Error allocating memory for Ngroup\n");
     exit(1);
   }

  // used to identify how many particles of each group to remove
   Nremove = new int[NO_GROUPS]();
   if(!Nremove){
     printf("Error allocating memory for Nremove\n");
     exit(1);
   }

}

/* Destructor: ~particle */
particle::~particle()
{
   delete [] q;
   delete [] v_limit;
   delete [] group_weight;
   delete [] group_limit;
   delete [] new_group_size;
   delete [] j_group;
   delete [] Ngroup;
   delete [] Nremove;

}


/*
* Method: set_vt
* Usage: set_vt(thermal_velocity)
* ---------------------------------    
* Set the particle thermal velocity, sqrt( (kT/e)*e/m).
*/
void particle::set_vt(double vt){
     this->vt=vt;
}

/*
* Method: set_q_over_m
* Usage:  set_q_over_m(mass_over_charge);
* ----------------------------------------
* Set the particle charge to mass ratio.                 
*/
void particle::set_q_over_m(double q_over_m){
     this->q_over_m=q_over_m;
}

/*
* Method: set_mass
* Usage:  set_mass(particle_mass);
* ---------------------------------
* Set the particle mass.
*/
void particle::set_mass(double mass){
     this->mass=mass;
}

/*
* Method: set_dxp
* Usage:  set_dxp(distance_between_particles);
* ---------------------------------------------
* Particles are loaded uniformly, at a distance dxp apart,
* at the start of the simulation if pic::restart=false.
*/
void particle::set_dxp(double dxp){
     this->dxp=dxp;
}

/*
* Method: set_dyp
* Usage:  set_dyp(distance_between_particles);
* ---------------------------------------------
* Particles are loaded uniformly, at a distance dxp apart,
* at the start of the simulation if pic::restart=false.
*/
void particle::set_dyp(double dyp){
     this->dyp=dyp;
}

/*
* Method: set_init_no
* Usage:  set_init_no(initial_no_of_particles)
* ---------------------------------------------
* Set the initial number of particles to load into the simulation if
* pic::restart = false.
*/ 
void particle::set_init_no(int init_particle_no){
     this->init_particle_no=init_particle_no;
}

/*
* Method: set_no_inject
* Usage:  set_no_inject(no_particles_to_inject);
* -----------------------------------------------
* Set the number of particles to inject each time-step. Currently,
* particles are inserted randomly along the x-axis. This variable is
* currently mandatory. For no particles, set it to zero.
*/
void particle::set_no_inject(int no_to_inject){
     this->no_inject = no_to_inject;
}

 /*
* Method: set_dt
* Usage:  set_dt(time_step);
* ----------------------------
* Set the time-step. The ion time-step can be optionally larger than
* the electron.
*/
void particle::set_dt(double dt){
     this->dt=dt;
}

/*
* Method: set_qc
* Usage:  set_qc(charge_of_each_macroParticle); 
* ----------------------------------------------
* Set the charge of each macro-particle. Typically this is equall to
* the size of each macro-particle multiplied by the charge.
*/
void particle::set_qc(double qc){
     this->qc=qc;
}

/*
* Method: set_no_groups
* Usage:  set_no_groups(no_groups);
* -----------------------------------------------
* Set the number of velocity groups to be used when decreasing
* the particle number.
*/
void particle::set_no_groups(int val){
     this->no_groups = val;
}

/*
* Method: set_v_limit
* Usage:  set_v_limit(v_limit,group_element);
* -----------------------------------------------
* Set the upper velocity of a velocity group. Used to increase particle
* weighting of particles that are first placed in user defined
* velocity groups.
*/
void particle::set_v_limit(double val,int element){
     v_limit[element]=val;
}

/*
* Method: set_group_limit
* Usage:  set_group_limit(group_limit,group_element);
* -----------------------------------------------
* Set the particle number limit of a velocity group. Used to increase 
* particle weighting of particles that are first placed in user defined
* velocity groups.
*/
void particle::set_group_limit(int val,int element){
     group_limit[element]=val;
}

/*
* Method: set_new_group_size
* Usage:  set_new_group_size(group_limit,group_element);
* -----------------------------------------------
* Set the particle desired particle number of a velocity group. Used to  
* increase particle weighting of particles that are first placed in user 
* defined velocity groups.
*/
void particle::set_new_group_size(int val,int element){
     new_group_size[element]=val;
}

/*
* Method: set_size_to_check
* Usage:  set_size_to_check(val);
* -----------------------------------------------
* Putting particles into groups to check each group size is time
* consuming so this is only done when the number of particles of a given
* species reaches the value given by size_to_check.
*/ 
void particle::set_size_to_check(int val){
     this->size_to_check=val;
}

/* Sets parameters needed for openCL. The parameters are initialized
* and destroyed in the PIC class and are sent here for use.
*/
void particle::init_openCL(cl_context context, cl_command_queue cmdQueue,
                     cl_kernel kernelCC, cl_kernel kernelE, cl_kernel kernelV,
                      cl_mem dQ, cl_mem dEx, cl_mem dEy )
{
    this->context = context;
    this->cmdQueue = cmdQueue;
    this->kernelCC = kernelCC;
    this->kernelE = kernelE;
    this->kernelV = kernelV;
    this->dQ = dQ;
    this->dEx = dEx;
    this->dEy = dEy;
}

/* Sets parameters needed for openCL for the case of a steadystate
* boundary. In this instance the particle arrays don't change size
* therefore memory need only be reserved during initialization.
*/
void particle::init_openCL_steady(cl_mem dP, cl_mem dExp, cl_mem dEyp)
{
    this->dP = dP;
    this->dExp = dExp;
    this->dEyp = dEyp;
}

/*
* Method: get_init_no
* Usage:  int init_particle_no = get_init_no();
* ------------------------
* Returns the initial number of particles to be loaded if not loading  
* particles from a previous run.
*/
int particle::get_init_no(void){
     return this->init_particle_no;
}

/*
* Method: create_particles
* Usage:  create_particles();
* -----------------------------
* Adds init_particle_no vector elements to the particle vectors before
* the simulation starts, if pic::restart=false.
*/
void particle::create_particles(void){
   int i;
   for(i=0;i<init_particle_no;i++){
       particleObj.push_back(particleVec());
       particleObj[i].x = 0.0;
       particleObj[i].y = 0.0;
       particleObj[i].vx = 0.0;
       particleObj[i].vy = 0.0;
       particleObj[i].vz = 0.0;
       particleObj[i].weight = 1.0;
   }
}

 /*
  * Method: is_equall_to
  * Usage:  is_equall_to(vector_to_be_changed,vector_to_copy,
  *                        first_vec_element_no,second_vec_element_no)
  * --------------------------------------------------------------------
  * Sets an element of the first vector equall to an element in the second.
  * This might not be a necessary method.
  */
void particle::is_equall_to(vector<particleVec> &obj1, vector<particleVec> &obj2,int i, int j){

    obj1[i].x = obj2[j].x;
    obj1[i].y = obj2[j].y;
    obj1[i].vx = obj2[j].vx;
    obj1[i].vy = obj2[j].vy;
    obj1[i].vz = obj2[j].vz;
    obj1[i].weight = obj2[j].weight;
}

/*
* Method: load_particles
* Usage:  load_particles();
* --------------------------
* Starts the simulation by loading particle data dumped from a previous
* run if the pic::restart is initialized to true via the pic constructor.
*/
void particle::load_particles(void){  
     int i=0;
     int last_particle;
     double x,y, vx,vy,vz,weight;
     ifstream in;

     if(!in){
       cout<<"Unable to open file"<<endl;
       exit(1);
     }

    // the filename ends indicating if it holds electron or ion data
     char filename[80]="dataDump_";
     strcat(filename,type);
     strcat(filename,".txt"); 
  
     in.open(filename);

     if (in.fail())
     {
         cout << "failed to open file\n";
         exit(1);
     }

    // read in the particle data
     
     while (in >> x >> y >> vx >> vy >> vz >> weight) {
      
       particleObj.push_back(particleVec());
       last_particle = particleObj.size() - 1;
       particleObj[last_particle].x = x;
       particleObj[last_particle].y = y;
       particleObj[last_particle].vx = vx;
       particleObj[last_particle].vy = vy;
       particleObj[last_particle].vz = vz;
       particleObj[last_particle].weight = weight;
    }
     
     cout << "number of particles read in is " << particleObj.size() << endl;
     set_init_no(particleObj.size());
    
     in.close(); 
}

/*
 * Implementation Notes: initiate_positions
 * -------------------------------------------  
 * Particles are placed evenly in the grid domain at a distance dxp apart
 * from its neighbour in x and dyp apart from its neighbour in y. 
 * other.
 */
void particle::initiate_positions(void)
{
  int i = 0;
  double x_pos = 0.0, y_pos = 0.0;
  double epsilon_x = 0.0; // dxp / 2.0;
  double epsilon_y = 0.0; // dyp / 2.0;

  // Note: the user decides what values to make dxp and dyp from the input file

  while(i < init_particle_no) // exit after have placed all particles
  {
      /*
      * start back at the beginning again and place a
      * second particle where the previous ones were laid
      * if there are more particles than what dxp and
      * dyp account for
      */
      
      y_pos = dyp;
     while(y_pos < YLENGTH + epsilon_y && i < init_particle_no)
     {
        
        x_pos = dxp;
        while(x_pos < XLENGTH + epsilon_x && i < init_particle_no)
        {
           particleObj[i].x = x_pos;
           particleObj[i].y = y_pos;
           x_pos += dxp;
           ++i;
        }
        y_pos += dyp;
     }
  }
}

/*
 * Implementation Notes: get_particle_velocity
 * ---------------------------------------------
 * This method is a temporary place holder for a more accurate way to 
 * calculate the particle velocity taken randomly off a particle velocity
 * distribution. In this case Maxwellian.
 */
double particle::get_particle_velocity(void)
{
    double R1, R2;
    double v = 0.0;

    do{
       R1 =(double)rand()/((double)RAND_MAX);
       R2 =(double)rand()/((double)RAND_MAX);
         v = vt*cos(2.0*PI*R2)*sqrt(abs(log(R1)));
    }while( abs(v) > vt*10000);

    return v;
}

/*
* Method: initiate_velocities
* Usage:  initiate_velocities();
* --------------------------------
* Initiates particle velocites from a Maxwellian distribution before the  
* simulation starts if pic::restart=false.
*/
void particle::initiate_velocities(void)
{
   int i;

   for(i=0; i< init_particle_no; i++){
      particleObj[i].vx = get_particle_velocity();
      particleObj[i].vy = get_particle_velocity();
      particleObj[i].vz = get_particle_velocity();
   } 
}

/*
* Method: inject_more_particles
* Usage:  inject_more_particles();
* ----------------------------------
* Inject no_inject number of particles each time-step randomly in the
* grid.
*/
void particle::inject_more_particles(void){
   int i;
   double random_numberX;
   double random_numberY;

   for(i=0;i<no_inject;i++){
       particleObj.push_back(particleVec());

       // put particles at a random location
       random_numberX = (double)rand() / ((double)RAND_MAX);
       random_numberY = (double)rand() / ((double)RAND_MAX);
       particleObj[(particleObj.size() - 1)].x = (DX + random_numberX * (XLENGTH - 2.0 * DX));
       particleObj[(particleObj.size() - 1)].y = (DY + random_numberY * (YLENGTH - 2.0 * DY));

       // get the velocities
       particleObj[(particleObj.size()-1)].vx = get_particle_velocity();
       particleObj[(particleObj.size()-1)].vy = get_particle_velocity();
       particleObj[(particleObj.size()-1)].vz = get_particle_velocity();
   }
}

 /*
* Method: initialize_map
* Usage:  initialize_map();
* ---------------------------
* Initializes values used to reduce the particle population and increase
* the particle wieghting. This function is called by 
* contain_particle_growth.
*/
void particle::initialize_map(void){
   int i;

   for(i=0;i<no_groups;i++){
      j_group[i] = 0;
      Ngroup[i] = 0;
      Nremove[i] = 0;
      group_weight[i]=1.0;
   }
}


/* Implementation Notes: reduce_population
 * -----------------------------------------
 * See the implementation notes for the contain_particle_growth method for
 * more information. To summarize particles are separated into velocity groups
 * according to the magnitude of their velocity and an amount of particles
 * determined by the user is taken away from each group and then the wieghting
 * of the particles left over in the group is increased. In this way the tails
 * of the velocity distribution can be maintained.
 */
void particle::reduce_population(void){
   int i,g, n = particleObj.size();
   int particles_to_remove=0;
   int R;
   int paws=0;
   double vxy;

  /* Calculate the number of particles in each velocity group */
   for(i=0; i< n; i++){

      vxy = sqrt((particleObj[i].vx)*(particleObj[i].vx) + (particleObj[i].vy)*(particleObj[i].vy));

   /* Ngroup holds the number of particles in each velocity group */
       g=0; 
       
      while( vxy > v_limit[g] ) {
           g++;
       }
 
       Ngroup[g]+=1; 
   }

   /* Calculate Nremove, the number of particles to delete from each velocity
    * group.
    */
   for(g=0;g<no_groups;g++){
      if(Ngroup[g] > group_limit[g])
            Nremove[g]= Ngroup[g] - new_group_size[g];
   }

   for(g=0;g<no_groups;g++){
   /* The total number of particles to remove from the particle vector */
      particles_to_remove+=Nremove[g];
   }
    
    i=0;
   
    while( i < particles_to_remove){

     /* Randomly select a particle from the particle vector */
      R =rand()%(particleObj.size()-1);

      vxy = sqrt((particleObj[R].vx)*(particleObj[R].vx) + (particleObj[R].vy)*(particleObj[R].vy));  

      /* Determine which velocity group the particle vector element at R 
         belongs to. */ 
      g=0; 

      //while( abs(particleObj[R].v) > v_limit[g] ){
      while( vxy > v_limit[g] ){
            g++;
      }

      /* j_group is a counter that keeps track of the number of particles
       * deleted from that group 
       */
      /* Remove vector particle element R if j_group < Nremove */
      if( j_group[g] < Nremove[g] ){
            j_group[g]++;
            swap(particleObj[R],particleObj[(particleObj.size()-1)]);
            particleObj.pop_back();
            i++;
      }
   }//while i<particles_to_remove

   /* Calculate the particle weight factor of each velocity group */

   if(particles_to_remove >0){
      for(g=0;g<no_groups;g++){
         if(j_group[g]>0){
          group_weight[g]= (double)Ngroup[g]/( (double)(Ngroup[g] - j_group[g]));
          }
      }
 
      n=particleObj.size();

   /* Increase particle wieght of all remaining particles */

      for(i=0;i<n;i++){
          g=0;

          vxy = sqrt((particleObj[i].vx)*(particleObj[i].vx) + (particleObj[i].vy)*(particleObj[i].vy));

          //while(abs(particleObj[i].v) > v_limit[g]){
         while(vxy> v_limit[g]){
            g++;
          }
        particleObj[i].weight = (particleObj[i].weight)*group_weight[g];
      }

  } //if particles_to_remove >0
}

/*
 * Implementation Notes: contain_particle_growth
 * -----------------------------------------------
 * Limit the size of the vector particle arrays by reducing the particle number
 * according to the scheme set out by E E Kunhardt et al:
 * -Journal of Computational Physics,vol 67,(1986) 279
 * -----------------------------------------------------
 * This scheme involves dividing the particles up into velocity groups,
 * and decreasing the size of each group once the upper group size is reached.
 * Each remaining particle is then rewieghted according to the amount that the 
 * group is reduced by. The number of particle velocity groups and the upper 
 * size limit of each group is determined by the user. 
 */
void particle::contain_particle_growth(void){
   if( particleObj.size() > size_to_check ){
        initialize_map();
        reduce_population();
        count_reduce++;
   }
}

/*
 * Implementation Notes: collect_charge
 * --------------------------------------
 * Calculates charge on the grid using the first order scheme presented in
 * C K Birdsall and A B Langdon, Plasma Physics via Computer Simulation, 
 * Institute of Physics Publishing, chapter 14, page 308
 * 
 */
void particle::collect_charge(void)  
{
   int i,j,k; 

   //Birdsall pg 308 use q/area
   double qc_dxdy = qc/(DX*DY); // macrocharge per area
   double dx = DX, dy = DY; // cell size
   int no_x_cells = X_NO_OF_CELLS, no_y_cells = Y_NO_OF_CELLS;
   int no_particles = particleObj.size(); // number of particles

    // size of the memory to reserve on the GPU for:
    // the charge density q array
   size_t qSize = Y_NO_OF_CELLS * X_NO_OF_CELLS * sizeof(cl_double);
   // for the particle vector
   size_t particleSizeVec = no_particles * sizeof(particleVec);
   

   // initialize q before send to openCL function
   for(k=0; k< Y_NO_OF_CELLS; k++){
      for(j=0; j< X_NO_OF_CELLS; j++){
         q[k*X_NO_OF_CELLS + j] = (cl_double)0.0;
      }    
   }

  /* Implement OpenCL code to vectorize this function onto the GPU. The
   * following steps have been taken from Mike Bailey's course on parallel
   * programming and he has given permission to replicate them in this
   * project. 
   * 
   * The following steps allocate enough memory to the GPU memory and write
   * the particle data arrays to this memory. The GPU then pulls each element
   * from the particle data arrays and applies the original function to collect
   * the charge at the grid points due to nearby particles. After the GPU has
   * finished collecting charge at the grid points it returns the charge density
   * array which is utilized elsewhere in the simulation.
   * 
   * To know what OpenCL commands to include in this function as compared to
   * else where the following resources were used:
   * https://www.olcf.ornl.gov/tutorials/opencl-vector-addition/
   * https://github.com/conanwu777/particle_system/blob/master/main.cpp
   * 
   */

   // 5. allocate the device memory buffers:

    if (STEADYSTATE == "NO") // then memory requirements change per iteration
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
    status = clEnqueueWriteBuffer(cmdQueue, dP, CL_FALSE, 0, particleSizeVec, &particleObj[0], 0, NULL, NULL);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clEnqueueWriteBuffer failed (1)\n");
        cout << "Status = " << status << endl;
    }

    // transfer a zero initialized grid charge array to GPU memory
    status = clEnqueueWriteBuffer(cmdQueue, dQ, CL_FALSE, 0, qSize, q, 0, NULL, NULL);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clEnqueueWriteBuffer failed (3)\n");
        cout << "Status = " << status << endl;
    }

    // ensure the above steps have completed
    Wait(cmdQueue);

    // 10. setup the arguments to the kernel object:

    // particle data array
    status = clSetKernelArg(kernelCC, 0, sizeof(cl_mem), &dP);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clSetKernelArg failed (1)\n");
        cout << "Status = " << status << endl;
    }

    // charge density array
    status = clSetKernelArg(kernelCC, 1, sizeof(cl_mem), &dQ);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clSetKernelArg failed (2)\n");
        cout << "Status = " << status << endl;
    }

    // macro-charge/area
    status = clSetKernelArg(kernelCC, 2, sizeof(cl_double), &qc_dxdy);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clSetKernelArg failed (3)\n");
        cout << "Status = " << status << endl;
    }

    // DX cell size in x
    status = clSetKernelArg(kernelCC, 3, sizeof(cl_double), &dx);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clSetKernelArg failed (4)\n");
        cout << "Status = " << status << endl;
    }

    // DY cell size in y
    status = clSetKernelArg(kernelCC, 4, sizeof(cl_double), &dy);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clSetKernelArg failed (5)\n");
        cout << "Status = " << status << endl;
    }

    // X_NO_CELLS
    status = clSetKernelArg(kernelCC, 5, sizeof(cl_int), &no_x_cells);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clSetKernelArg failed (6)\n");
        cout << "Status = " << status << endl;
    }

    // Y_NO_CELLS
    status = clSetKernelArg(kernelCC, 6, sizeof(cl_int), &no_y_cells);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clSetKernelArg failed (7)\n");
        cout << "Status = " << status << endl;
    }


    // 11. enqueue the kernel object for execution:

    size_t globalWorkSize[3] = { no_particles , 1, 1 };
    size_t localWorkSize[3] = { LOCAL_SIZE,   1, 1 };

    Wait(cmdQueue);

    status = clEnqueueNDRangeKernel(cmdQueue, kernelCC, 1, NULL, globalWorkSize, NULL, 0, NULL, NULL);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clEnqueueNDRangeKernel CC failed: %d\n", status);
        cout << "Status = " << status << endl;
    }


    Wait(cmdQueue); 

    // 12. read the results buffer back from the device to the host:

    status = clEnqueueReadBuffer(cmdQueue, dQ, CL_TRUE, 0, qSize, q, 0, NULL, NULL);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clEnqueueReadBuffer failed\n");
        cout << "Status = " << status << endl;
    }


    // do below outside of vectorization

    // double the value of q along the sides
    for (int j = 0; j < Y_NO_OF_CELLS; j++)
    {
        q[j * X_NO_OF_CELLS] *= 2.0;
        q[j * X_NO_OF_CELLS + X_NO_OF_CELLS - 1] *= 2.0;
    }
    for (int i = 0; i < X_NO_OF_CELLS; i++)
    {
        q[i] *= 2.0;
        q[(Y_NO_OF_CELLS - 1) * X_NO_OF_CELLS + i] *= 2.0;
    }

    //double the value again in the corners
    q[0] *= 2.0;
    q[X_NO_OF_CELLS - 1] *= 2.0;
    q[(Y_NO_OF_CELLS - 1) * X_NO_OF_CELLS] *= 2.0;
    q[(Y_NO_OF_CELLS - 1) * X_NO_OF_CELLS + X_NO_OF_CELLS - 1] *= 2.0;


    // At the end of each function call flush the queue
    // reference: https://github.com/conanwu777/particle_system/blob/master/main.cpp
    clFlush(cmdQueue);

    if (STEADYSTATE == "NO") // release variable sized reserved memory
    {
        clReleaseMemObject(dP);

    }

}



/*
* 2d grid charge-in-cloud weighting of grid electric field to particle electric field
* https://people.nscl.msu.edu/~lund/uspas/sbp_2018/lec_intro/04.pic.pdf and
* Birdsal and Langdon in Plasma Physics via Computer Simulation Chapter 14
*
* This calculate what the nearby values of the electric field on the grid feels like
* at each specific particle position.
*/
void particle::get_particle_electric_field(double *EgridX, double *EgridY)
{
   int no_particles = particleObj.size(); // number of particles
   double dx = DX, dy = DY; // cell size
   int no_x_cells = X_NO_OF_CELLS, no_y_cells = Y_NO_OF_CELLS;

    /*
     * The electric field at the particle needs to vary with the size of the particle
     * array. Since this value is not used outside of the particle class it is declared
     * as an array of openCL doubles before being passed to openCL.
     */
    EfieldXp = new cl_double[no_particles];
    if (!EfieldXp) {
        printf("Error allocating memory for EfieldXp\n");
    }

    EfieldYp = new cl_double[no_particles];
    if (!EfieldYp) {
        printf("Error allocating memory for EfieldYp\n");
    }

    // Initialize the electric field each iteration
    for (int i = 0; i < no_particles; i++)
    {
        EfieldXp[i] = (cl_double)0.0;
        EfieldYp[i] = (cl_double)0.0;
    }

    // Size of the memory area to reserve on the GPU for:
    // the electric field on the grid
    size_t EgridSize = (X_NO_OF_CELLS + 1) * (Y_NO_OF_CELLS + 1) *sizeof(cl_double);
    // particles array
    size_t particleSizeVec = no_particles * sizeof(particleVec);
    // the field at the particles
    size_t fieldParticleSize = no_particles * sizeof(cl_double);


   // 5. allocate the device memory buffers:

   if (STEADYSTATE == "NO") // then memory requirements change per iteration
   {
       // holds the particle data
       dP = clCreateBuffer(context, CL_MEM_READ_ONLY, particleSizeVec, NULL, &status);
       if (status != CL_SUCCESS)
       {
           fprintf(stderr, "clCreateBuffer failed (1)\n");
           cout << "Status = " << status << endl;
       }

        // will store the resulting x component of the electric field
       dExp = clCreateBuffer(context, CL_MEM_WRITE_ONLY, fieldParticleSize, NULL, &status);
       if (status != CL_SUCCESS)
       {
           fprintf(stderr, "clCreateBuffer failed (2)\n");
           cout << "Status = " << status << endl;
       }

        // will store the resulting y component of the electric field
       dEyp = clCreateBuffer(context, CL_MEM_WRITE_ONLY, fieldParticleSize, NULL, &status);
       if (status != CL_SUCCESS)
       {
           fprintf(stderr, "clCreateBuffer failed (3)\n");
           cout << "Status = " << status << endl;
       }
   }

   // 6. enqueue the 2 commands to write the data from the host buffers to the device buffers:

    // transfer the local particle data to GPU memory
   status = clEnqueueWriteBuffer(cmdQueue, dP, CL_FALSE, 0, particleSizeVec, &particleObj[0], 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueWriteBuffer failed (1)\n");
       cout << "Status = " << status << endl;
   }

   // transfer a zero initialized grid charge array to GPU memory
   status = clEnqueueWriteBuffer(cmdQueue, dEx, CL_FALSE, 0, EgridSize, EgridX, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueWriteBuffer failed (2)\n");
       cout << "Status = " << status << endl;
   }

   // transfer a zero initialized grid charge array to GPU memory
   status = clEnqueueWriteBuffer(cmdQueue, dEy, CL_FALSE, 0, EgridSize, EgridY, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueWriteBuffer failed (3)\n");
       cout << "Status = " << status << endl;
   }

   // ensure the above steps have completed
   Wait(cmdQueue);

   // 10. setup the arguments to the kernel object:

   // particle data array
    status = clSetKernelArg(kernelE, 0, sizeof(cl_mem), &dP);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clSetKernelArg failed (1)\n");
        cout << "Status = " << status << endl;
    }

    // the x component of the electric field at the particle  
    status = clSetKernelArg(kernelE, 1, sizeof(cl_mem), &dExp);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clSetKernelArg failed (2)\n");
        cout << "Status = " << status << endl;
    }

    // the y component of the electric field at the particle
    status = clSetKernelArg(kernelE, 2, sizeof(cl_mem), &dEyp);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clSetKernelArg failed (3)\n");
        cout << "Status = " << status << endl;
    }

    // the x component of the electric field at the grid
    status = clSetKernelArg(kernelE, 3, sizeof(cl_mem), &dEx);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clSetKernelArg failed (4)\n");
        cout << "Status = " << status << endl;
    }

   // the y component of the electric field at the grid
   status = clSetKernelArg(kernelE, 4, sizeof(cl_mem), &dEy);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (5)\n");
       cout << "Status = " << status << endl;
   }

   // DX cell size in x
   status = clSetKernelArg(kernelE, 5, sizeof(cl_double), &dx);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (6)\n");
       cout << "Status = " << status << endl;
   }

   // DY cell size in y
   status = clSetKernelArg(kernelE, 6, sizeof(cl_double), &dy);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (7)\n");
       cout << "Status = " << status << endl;
   }

   // X_NO_CELLS
   status = clSetKernelArg(kernelE, 7, sizeof(cl_int), &no_x_cells);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (8)\n");
       cout << "Status = " << status << endl;
   }

   // Y_NO_CELLS
   status = clSetKernelArg(kernelE, 8, sizeof(cl_int), &no_y_cells);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (9)\n");
       cout << "Status = " << status << endl;
   }


   // 11. enqueue the kernel object for execution:

   size_t globalWorkSize[3] = { no_particles , 1, 1 };
   size_t localWorkSize[3] = { LOCAL_SIZE,   1, 1 };

   Wait(cmdQueue);

   status = clEnqueueNDRangeKernel(cmdQueue, kernelE, 1, NULL, globalWorkSize, NULL, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueNDRangeKernel E failed: %d\n", status);
       cout << "Status = " << status << endl;
       cout << "local size = " << LOCAL_SIZE << endl;
   }


   Wait(cmdQueue); 

 // 12. read the results buffer back from the device to the host:

   status = clEnqueueReadBuffer(cmdQueue, dExp, CL_TRUE, 0, fieldParticleSize, EfieldXp, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueReadBuffer failed (1)\n");
       cout << "Status = " << status << endl;
   }

   status = clEnqueueReadBuffer(cmdQueue, dEyp, CL_TRUE, 0, fieldParticleSize, EfieldYp, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueReadBuffer failed (2)\n");
       cout << "Status = " << status << endl;
   }

   // At the end of each function call flush the queue
 // reference: https://github.com/conanwu777/particle_system/blob/master/main.cpp
   clFlush(cmdQueue);

   if (STEADYSTATE == "NO")  // release variable sized reserved memory
   {
       clReleaseMemObject(dP);
       clReleaseMemObject(dExp);
       clReleaseMemObject(dEyp);
   }
  
}

/*   
 *
 * Implementation Notes: move_particles
 * ----------------------------------------
 * Calculates the new particle velocities and positions due to the electric
 * and optional magnetic fields using the leap-frog and  
 * accelerate-rotate-accelerate methods presented in:
 * C K Birdsall and A B Langdon, Plasma Physics via Computer Simulation, 
 * Institute of Physics Publishing, chapter 4.
 * 
 * 2d grid charge-in-cloud weighting of grid electric field to particle electric field
 * https://people.nscl.msu.edu/~lund/uspas/sbp_2018/lec_intro/04.pic.pdf and
 * Birdsall and Langdon chapter 14
 */                                                                                                                                               
void particle::move_particles(){
   int no_particles = particleObj.size(); // number of particles
   double coef = 0.5*dt*q_over_m; // coeficient for the equation of motion due to an electric field
   double omega_dt = dt*q_over_m*MAGNETIC_FIELD; // coeficient due to movement in a magnetic field
   double magnetic_field = MAGNETIC_FIELD; // value of the magnetic field strength
   // angles involved in rotating particles about a magnetic field line B
   double cos_angle = cos(ANGLE);
   double sin_angle = sin(ANGLE);
   double cos_omega = cos(omega_dt);
   double sin_omega = sin(omega_dt);
   double dx = DX, dy = DY;  // cell sizes
   int no_x_cells = X_NO_OF_CELLS, no_y_cells = Y_NO_OF_CELLS;
   double vperp,vpar,vperp_after,ux, uy; //particle velocities perpendicular and parallel to B

    // size of the memory areas to reserve on the gpu
   size_t particleSizeVec = no_particles * sizeof(particleVec);
   size_t fieldParticleSize = no_particles * sizeof(cl_double);

   if (STEADYSTATE == "NO") // then memory requirements change per iteration
   {
       // holds the particle data
       dP = clCreateBuffer(context, CL_MEM_READ_WRITE, particleSizeVec, NULL, &status);
       if (status != CL_SUCCESS)
       {
           fprintf(stderr, "clCreateBuffer failed (1)\n");
           cout << "Status = " << status << endl;
       }

        // will store the resulting x component of the electric field
       dExp = clCreateBuffer(context, CL_MEM_READ_ONLY, fieldParticleSize, NULL, &status);
       if (status != CL_SUCCESS)
       {
           fprintf(stderr, "clCreateBuffer failed (3)\n");
           cout << "Status = " << status << endl;
       }

        // will store the resulting y component of the electric field
       dEyp = clCreateBuffer(context, CL_MEM_READ_ONLY, fieldParticleSize, NULL, &status);
       if (status != CL_SUCCESS)
       {
           fprintf(stderr, "clCreateBuffer failed (4)\n");
           cout << "Status = " << status << endl;
       }

   }

   // 6. enqueue the 2 commands to write the data from the host buffers to the device buffers:

   // transfer the local particle data to GPU memory

   status = clEnqueueWriteBuffer(cmdQueue, dP, CL_FALSE, 0, particleSizeVec, &particleObj[0], 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueWriteBuffer failed (1)\n");
       cout << "Status = " << status << endl;
   }

   status = clEnqueueWriteBuffer(cmdQueue, dExp, CL_FALSE, 0, fieldParticleSize, EfieldXp, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueReadBuffer failed (3)\n");
       cout << "Status = " << status << endl;
   }

   status = clEnqueueWriteBuffer(cmdQueue, dEyp, CL_FALSE, 0, fieldParticleSize, EfieldYp, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueReadBuffer failed (4)\n");
       cout << "Status = " << status << endl;
   }

   // ensure the above steps have completed
   Wait(cmdQueue);


   // 10. setup the arguments to the kernel object:

   // particle data array
   status = clSetKernelArg(kernelV, 0, sizeof(cl_mem), &dP);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (0)\n");
       cout << "Status = " << status << endl;
   }

   // the x component of the electric field at the particle 
   status = clSetKernelArg(kernelV, 1, sizeof(cl_mem), &dExp);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (1)\n");
       cout << "Status = " << status << endl;
   }

    // the y component of the electric field at the particle
   status = clSetKernelArg(kernelV, 2, sizeof(cl_mem), &dEyp);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (2)\n");
       cout << "Status = " << status << endl;
   }

   // DX cell size in x
   status = clSetKernelArg(kernelV, 3, sizeof(cl_double), &dx);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (3)\n");
       cout << "Status = " << status << endl;
   }

   // DY cell size in y
   status = clSetKernelArg(kernelV, 4, sizeof(cl_double), &dy);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (4)\n");
       cout << "Status = " << status << endl;
   }

   // X_NO_CELLS
   status = clSetKernelArg(kernelV, 5, sizeof(cl_int), &no_x_cells);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (5)\n");
       cout << "Status = " << status << endl;
   }

   // Y_NO_CELLS
   status = clSetKernelArg(kernelV, 6, sizeof(cl_int), &no_y_cells);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (6)\n");
       cout << "Status = " << status << endl;
   }

   // coeficient for the equation of motion due to an electric field 
   status = clSetKernelArg(kernelV, 7, sizeof(cl_double), &coef);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (7)\n");
       cout << "Status = " << status << endl;
   }

    // coeficient due to movement in a magnetic field
   status = clSetKernelArg(kernelV, 8, sizeof(cl_double), &omega_dt);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (8)\n");
       cout << "Status = " << status << endl;
   }

   // angles involved in rotating particles about a magnetic field line B

   status = clSetKernelArg(kernelV, 9, sizeof(cl_double), &cos_angle);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (9)\n");
       cout << "Status = " << status << endl;
   }

   status = clSetKernelArg(kernelV, 10, sizeof(cl_double), &sin_angle);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (10)\n");
       cout << "Status = " << status << endl;
   }

    //time step
   status = clSetKernelArg(kernelV, 11, sizeof(cl_double), &dt);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (11)\n");
       cout << "Status = " << status << endl;
   }

    // angles involved in rotating particles about a magnetic field line B
   status = clSetKernelArg(kernelV, 12, sizeof(cl_double), &cos_omega);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (12)\n");
       cout << "Status = " << status << endl;
   }

   status = clSetKernelArg(kernelV, 13, sizeof(cl_double), &sin_omega);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (13)\n");
       cout << "Status = " << status << endl;
   }

   // value of the magnetic field strength 
   status = clSetKernelArg(kernelV, 14, sizeof(cl_double), &magnetic_field);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clSetKernelArg failed (14)\n");
       cout << "Status = " << status << endl;
   }

   // 11. enqueue the kernel object for execution:

   size_t globalWorkSize[3] = { no_particles , 1, 1 };
   size_t localWorkSize[3] = { LOCAL_SIZE,   1, 1 };

   Wait(cmdQueue);

   status = clEnqueueNDRangeKernel(cmdQueue, kernelV, 1, NULL, globalWorkSize, NULL, 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueNDRangeKernel V failed: %d\n", status);
       cout << "Status = " << status << endl;
   }


   Wait(cmdQueue); 

 // 12. read the results buffer back from the device to the host:

   status = clEnqueueReadBuffer(cmdQueue, dP, CL_TRUE, 0, particleSizeVec, &particleObj[0], 0, NULL, NULL);
   if (status != CL_SUCCESS)
   {
       fprintf(stderr, "clEnqueueReadBuffer failed (1)\n");
       cout << "Status = " << status << endl;
   }

   // At the end of each function call flush the queue
    // reference: https://github.com/conanwu777/particle_system/blob/master/main.cpp
   clFlush(cmdQueue);

   if (STEADYSTATE == "NO") // release variable sized reserved memory
   {
       clReleaseMemObject(dP);
       clReleaseMemObject(dExp);
       clReleaseMemObject(dEyp);
   }

    // Delete the electric field at the particles arrays. These arrays
    // change size with the particle size.
    delete[] EfieldXp;
    EfieldXp = NULL;

    delete[] EfieldYp;
    EfieldYp = NULL;
}

/*
* Method: set_type
* Usage:  set_type("i1");
* -------------------------
* Sets a string expression that is used to identify each species type
* being dumped in a text file or read in for particle loading.
*/
void particle::set_type(char *type_in){
     type = type_in;
}

/*
* Method: dump_data
* Usage:  dump_data();
* ---------------------
* Dumps particle data to a text file for each species at the end of a run
* if pic::dump == true. Note, that this is set to true by default.
*/
void particle::dump_data(void){
     int i;
     int n = particleObj.size();
     ofstream out;

     if(!out){
       cout<<"Unable to open file"<<endl;
       exit(1);
     }

     char filename[80]="dataDump_";
     strcat(filename,type);  // the end of the file name indicates the particle type
     strcat(filename,".txt"); 
  
     out.open(filename);

     for(i=0;i< n;i++){
       out<<particleObj[i].x<<" "<<particleObj[i].y<<" "<<particleObj[i].vx<<" "<<particleObj[i].vy<<" "<<particleObj[i].vz<<" "<<particleObj[i].weight<<endl;
    }

     out.close();
}

