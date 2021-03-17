/*
 * File: pic.cpp
 * ------------
 * 
 * 
 * This class is the engine of the simulation. It is responsible for
 * initializing the data, running the simulation in a time loop and doing
 * closing up duties such as data dumps. The pic class initializes
 *  and runs a particle-in-cell simulation by declaring objects
 * and calling the methods defined in the boundary, particle, fields,  
 * collisions, e_collisions, ion_collisions, and graphs classes.
 *
 * 10/15
 * 
 * Date: 3/1/2021
 * 
 * Add functions to retreive graph averaged plot data and particle data
 * for use in graphics output displays.
 * 
 * Date: 1/28/2021
 * 
 * Description: Update file to account for a two-dimensional simulation. This 
 * file is pretty insular from what it calls and so the only change to 
 * the simulation runtime was the addition of a function call to generate
 * particle fields from grid fields and a few function calls to dump more
 * graph data to output files.
 * 
 * In addition the initialization of simulation data from user input was
 * added to the pic constructor.
 * 
 * The infrastructure is in place to account for more than one ion species
 * but this is not yet implemented in all classes.
 * 
 * Initialize and destruct openCL variables in this file for other classes to use.
 *
 */

#include "pic.h"




/*
 * Constructor: pic
 * ------------------    
 * Initializes a particle-in-cell simulation as described in:
 *
 * -C K Birdsall and A B Langdon, Plasma Physics Via Computer Simulation
 * and the software on www.eecs.berkeley.edu
 *
 * Note: The variable restart can only be initialized via the constructor, and
 * is false by default.
 *
 */
pic::pic(bool restart_p,bool dump_p,int run_time_p,int start_graph_p,int end_graph_p):restart(restart_p),dump(dump_p),run_time(run_time_p),start_graph(start_graph_p),end_graph(end_graph_p){

   int i,ion_particle_no=0;
   writeGraphDataToFiles = true;  // is set to true by default

    // an array of particle objects; one for each ion species
   ions = new particle[NO_ION_SPECIES]();
   if (!ions){
      printf("Error allocating memory for ions\n");
      exit(1);
   }

    // an electron collision object
   e_collision = new e_collisions[NO_ION_SPECIES]();
   if (!e_collision){
      printf("Error allocating memory for e_collision\n");
      exit(1);
   }

    // an ion collision object
   ion_collision = new ion_collisions[NO_ION_SPECIES]();
   if (!ion_collision){
      printf("Error allocating memory for ion_collision\n");
      exit(1);
   }

    // an array of ion boundary objects; one for each species
   sec_ions = new boundary[NO_ION_SPECIES]();
   if (!sec_ions){
      printf("Error allocating memory for sec_ions\n");
      exit(1);
   }
   
    // an array of graph objects; one for each species
   graphs_ions = new graphs[NO_ION_SPECIES]();
   if (!graphs_ions){
      printf("Error allocating memory for graphs_ions\n");
      exit(1);
   }
   
   /* Initialize parameters from user input for the different classes:
   * particle, collision, boundary and graph
   */
   initialize_ion_particles();
   initialize_electron_particles();
   initialize_ion_collision();
   initialize_electron_collision();
   initialize_boundary();
   initialize_graph_parameters();
   initialize_OpenCL();

    // initialize with user input 
   haveCollisions = HAVE_COLLISIONS;

/* load particles into simulation region */
   if(restart==true){   // if user is starting simulation from a previous dump
       cout << "loading particles....\n";

      electrons.particle::load_particles();
      cout << "loaded electrons\n";

      for(i=0;i<NO_ION_SPECIES;i++)
         ions[i].particle::load_particles();
      cout << "loaded ions\n";
   }
   else{  // user is starting simulation from time zero
      for(i=0;i< NO_ION_SPECIES;i++){
         ions[i].particle::create_particles();
         ions[i].particle::initiate_positions();
         ions[i].particle::initiate_velocities();
      }
   electrons.particle::create_particles();
   electrons.particle::initiate_positions();
   electrons.particle::initiate_velocities();
   }

  


/* initialize neutral velocity arrays */
   for(i=0;i<NO_ION_SPECIES;i++){
      ion_collision[i].ion_collisions::calc_neutral_velocity();
      e_collision[i].e_collisions::calc_neutral_velocity();
   }


/* Variables used to store data variations in time */

   timedata.open("timeData.txt");
   if(!timedata){
      cout<<"Cannot open file"<<endl;
      exit(1);
   }

   time_plot = new int[SIZE_TIME_DATA]();
   if (!time_plot){
       printf("Error allocating memory for time_plot\n");
       exit(1);
   }

   for (int i = 0; i < SIZE_TIME_DATA; i++)
       time_plot[i] = 0.0;

    // number of electrons at one position
   no_e = new int[SIZE_TIME_DATA]();  
   if (!no_e){
       printf("Error allocating memory for no_e\n");
       exit(1);
   }

    // number of ions at one position
   no_i = new int*[NO_ION_SPECIES]();
   if (!no_i){
       printf("Error allocating memory for no_i\n");
       exit(1);
   }
   for(i=0;i<NO_ION_SPECIES;i++){
        no_i[i]= new int[SIZE_TIME_DATA];
        if (!no_i[i]){
             printf("Error allocating memory for no_i[i]\n");
             exit(1);
        }
   }

    // potential to left of grid
   phi_left = new double[SIZE_TIME_DATA]();
   if (!phi_left){
       printf("Error allocating memory for phi_left\n");
       exit(1);
   }

    // potential in middle of grid
   phi_mid = new double[SIZE_TIME_DATA]();
   if (!phi_mid){
       printf("Error allocating memory for phi_mid\n");
       exit(1);
   }

    // potential to right of grid
  phi_right = new double[SIZE_TIME_DATA]();
  if (!phi_right){
       printf("Error allocating memory for phi_right\n");
       exit(1);
   }

     // electron denisty to right of grid
   eden_right = new double[SIZE_TIME_DATA]();
   if (!eden_right){
       printf("Error allocating memory for eden_right\n");
       exit(1);
   }

    // ion denisty to right of grid
   iden_right = new double*[NO_ION_SPECIES]();
   if (!iden_right){
       printf("Error allocating memory for iden_right\n");
       exit(1);
   }
   for(i=0;i<NO_ION_SPECIES;i++){
       iden_right[i]= new double[SIZE_TIME_DATA];
       if (!iden_right[i]){
           printf("Error allocating memory for iden_right[i]\n");
           exit(1);
       }
   }

}

/* Initialize electron particle data with user defined input. 
*  All the data is defined in the user input file.
*/
void pic::initialize_electron_particles(void)
{
    //general particle data
    electrons.particle::set_init_no(INIT_PARTICLE_NO); // starting particle number
    electrons.particle::set_dt(DT);  // time step size
    electrons.particle::set_vt(VTE); // thermal velocity
    electrons.particle::set_q_over_m(Q_OVER_ME); // charge over mass
    electrons.particle::set_qc(QCE);  // macrocharge value
    electrons.particle::set_mass(ELECTRON_MASS); // particle mass
    electrons.particle::set_dxp(DXPE); // distance between particles for initial loading
    electrons.particle::set_no_inject(NUMBER_TO_INJECT); // number of particles to inject each timestep

    /* 
    * parameters used for reducing the particle number when the particle vectors
    * get too large:
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

    /* Number of velocity groups to account for when reducing the particle number.  
    * Choose more than one velocity group to preserve high velocity tails in the
    * random reduction of the number of particles. This is a mandatory constant 
    * that is also used to construct group arrays.
    */     
    electrons.particle::set_no_groups(NO_GROUPS);

    //Set for three velocity groups to use in reducing the number
    //of particles in the simulation
    electrons.particle::set_v_limit(VE1, 0);
    electrons.particle::set_group_limit(UPPER_GROUP_LIMIT, 0);
    electrons.particle::set_new_group_size(NEW_GROUP_SIZE, 0);
    electrons.particle::set_size_to_check(SIZE_TO_CHECK_E);

    electrons.particle::set_v_limit(VE2, 1);
    electrons.particle::set_group_limit(UPPER_GROUP_LIMIT, 1);
    electrons.particle::set_new_group_size(NEW_GROUP_SIZE, 1);

    electrons.particle::set_v_limit(UPPER_V_LIMIT, 2);
    electrons.particle::set_group_limit(UPPER_GROUP_LIMIT, 2);
    electrons.particle::set_new_group_size(NEW_GROUP_SIZE, 2);

    //Used to name data dump at end of simulation and output data files 
    electrons.particle::set_type(ELECTRON_TYPE);
}

/* Initialize ion particle data with user defined input. 
*  All the data is defined in the user input file.
*/
void pic::initialize_ion_particles(void)
{
    //general particle data
    ions[0].particle::set_vt(VTI); // thermal velocity
    ions[0].particle::set_q_over_m(Q_OVER_MI); // charge over mass
    ions[0].particle::set_mass(ION_MASS); // particle mass
    ions[0].particle::set_dxp(DXPI); // distance between particles for initial loading
    ions[0].particle::set_init_no(INIT_PARTICLE_NO); // starting particle number
    ions[0].particle::set_dt(DTI);  // time step size
    ions[0].particle::set_qc(QCI); // macrocharge value
    ions[0].particle::set_no_inject(NUMBER_TO_INJECT); // number of particles to inject each timestep

    /* 
    * parameters used for reducing the particle number when the particle vectors
    * get too large:
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

    /* Number of velocity groups to account for when reducing the particle number.  
    * Choose more than one velocity group to preserve high velocity tails in the
    * random reduction of the number of particles. This is a mandatory constant 
    * that is also used to construct group arrays.
    */    
    ions[0].particle::set_no_groups(NO_GROUPS);

    //Set for three velocity groups to use in reducing the number
    //of particles in the simulation
    ions[0].particle::set_v_limit(VI1, 0);
    ions[0].particle::set_group_limit(UPPER_GROUP_LIMIT, 0);
    ions[0].particle::set_new_group_size(NEW_GROUP_SIZE, 0);
    ions[0].particle::set_size_to_check(SIZE_TO_CHECK);

    ions[0].particle::set_v_limit(VI2, 1);
    ions[0].particle::set_group_limit(UPPER_GROUP_LIMIT2, 1);
    ions[0].particle::set_new_group_size(NEW_GROUP_SIZE2, 1);

    ions[0].particle::set_v_limit(UPPER_V_LIMIT_I, 2);
    ions[0].particle::set_group_limit(UPPER_GROUP_LIMIT2, 2);
    ions[0].particle::set_new_group_size(NEW_GROUP_SIZE2, 2);

    //Used to name data dump at end of simulation and output data files 
    ions[0].particle::set_type(ION_TYPE);
    
}

/*
* Initialize ion collision data using user defined input.
*  All the data is defined in the user input file.
*/
void pic::initialize_ion_collision(void)
{

    // general parameters
    ion_collision[0].ion_collisions::set_neutral_density(NEUTRAL_DENSITY);
    ion_collision[0].ion_collisions::set_vtn(VTN); // neutral thermal velocity
    ion_collision[0].ion_collisions::set_ion_mass(ION_MASS); //species mass

    //set the following to true if the collision cross-section calculation
    //requires the center of mass energy
    ion_collision[0].set_CMtrue_CX(false);

    //define the collision cross sections to be used. 

    /* ion collisions:
       predefined options currently include:
         sigma_Ar_cx();
         sigma_H_cx();
         sigma_Ar_elastic();
         sigma_H_elastic();
         For a constant value or no collisions need to set constant_elastic or
         constant_cx = value or zero.
    */

    //set ion elastic collision cross-section to zero for no collision
    if (SIGMA == "CONSTANT")
    {
        // const sigma or no collision
        ion_collision[0].ion_collisions::set_constant_CX(SIGMA_VALUE);
        ion_collision[0].ion_collisions::set_crossSec_CX(&ion_collisions::sigma_const_cx);
    }
    else if (SPECIES == "ARGON")
    {
        // energy dependent cross section
        ion_collision[0].ion_collisions::set_crossSec_CX(&ion_collisions::sigma_Ar_cx);
    }
    else
    {
        // energy dependent cross section
        ion_collision[0].ion_collisions::set_crossSec_CX(&ion_collisions::sigma_H_cx);
    }
      
}

/*
* Initialize electron collision data using user defined input.
*  All the data is defined in the user input file.
*/
void pic::initialize_electron_collision(void)
{
    // general parameters
    e_collision[0].e_collisions::set_neutral_density(NEUTRAL_DENSITY);
    e_collision[0].e_collisions::set_vtn(VTN); // neutral thermal velocity
    e_collision[0].e_collisions::set_ion_mass(ION_MASS); //species mass

    // threshold electron energy for ionization to occur
    e_collision[0].e_collisions::set_ionization_threshold(THRESHOLD);

    /* electron collisions
    predefined options currently include:
    sigma_Ar_elastic()
    sigma_Ar_ionize()
    sigma_H_elastic()
    sigma_H_ionize()
    sigma_const_elastic()
    sigma_const_ionization()
    user defined functions
     For a constant value or no collisions need to set constant_elastic or
     constant_ionization = value or zero. 
    */

   // set ELASTIC_SIGMA_VALUE to zero for no collisions 
    if (SIGMA_ELASTIC == "CONSTANT")
    {
        e_collision[0].e_collisions::set_constant_elastic(ELASTIC_SIGMA_VALUE);
        e_collision[0].e_collisions::set_crossSec_elastic(&e_collisions::sigma_const_elastic);
    }
    else if (SPECIES == "ARGON") //implement collisions with energy dependent cross-section
        e_collision[0].e_collisions::set_crossSec_elastic(&e_collisions::sigma_Ar_elastic);
    else // hydrogen, again with energy dependent cross-section
        e_collision[0].e_collisions::set_crossSec_elastic(&e_collisions::sigma_H_elastic);
    
    // set IONIZATION_SIGMA_VALUE to zero for no collisions
    if (SIGMA_IONIZATION == "CONSTANT")
    {
        e_collision[0].e_collisions::set_constant_ionization(IONIZATION_SIGMA_VALUE);
        e_collision[0].e_collisions::set_crossSec_ionization(&e_collisions::sigma_const_ionization);
    }
    else if (SPECIES == "ARGON")  //implement collisions with energy dependent cross-section
        e_collision[0].e_collisions::set_crossSec_ionization(&e_collisions::sigma_Ar_ionize);
    else   // hydrogen, again with energy dependent cross-section
        e_collision[0].e_collisions::set_crossSec_ionization(&e_collisions::sigma_H_ionize);
    
}

/*
* Initialize boundary data using user defined input.
*  All the data is defined in the user input file.
*/
void pic::initialize_boundary(void)
{
    
    //Set the implement_x_boundary function pointers:
    /*
     *  There are currently two options:
     *  1)e_boundary and ion_boundary
     *         This includes reflection, backscattering, absorption, secondary
     *          -electrons.
     *  2)steady_e_boundary and steady_i_boundary
     *         Implements a rapid steady-state by inserting an electron-ion pair
     *         into the simulation domain everytime an electron arrives at the
     *         surface. Ionization should be turned off, and it is recommended
     *         that the ion and electron time-step are the same. The particle
     *         number remains constant throughtout the simulation so it is
     *         better to make INIT_PARTICLE_NO large.
     */

     //Set parameters for electron induced secondary emmision, electron reflection
     //and backscattering

     //secondary emmision threshold in Vaughan secondary emmision model
    sec_e.boundary::set_E0(4.32);

    //max yield cooresponding to max incoming energy Emax0
    sec_e.boundary::set_max_yield(1.4);
    sec_e.boundary::set_Emax0(650.0);
    sec_e.boundary::set_incident_mass(ELECTRON_MASS);

    //Ionization threshold of ion species
    sec_e.boundary::set_threshold(THRESHOLD);

    //Work function of surface
    sec_e.boundary::set_work_function(WORK_FUNCTION);
    sec_e.boundary::set_vt(VTE);  //electron thermal velocity

    //set parameters for ion induced secondary electron emmission
    //max yield cooresponding to max incoming energy Emax0
    sec_ions[0].boundary::set_max_yield(1.4);
    sec_ions[0].boundary::set_E0(4.32); // secondary emmision threshold
    sec_ions[0].boundary::set_Emax0(650.0);

    sec_ions[0].boundary::set_incident_mass(ION_MASS); //ion mass
    sec_ions[0].boundary::set_threshold(THRESHOLD); // threshold energy

    //Work function of surface
    sec_ions[0].boundary::set_work_function(WORK_FUNCTION);
    sec_ions[0].boundary::set_vt(VTI); // ion thermal velocity

    //set function pointers to calculate yields
    /* current options include:

         e_impact_see_yield_v();  // get electron impact secondary electrion yield
         i_impact_see_yield_Ar(); // get ion impact secondary electrion yield
         calc_reflection_yield_v(); // get reflection yield
         calc_backscatter_yield_v(); // get backscatter yield
    */

    if (E_IMPACT_SEE == "YES" )
        sec_e.boundary::set_e_impact_yield(&boundary::e_impact_see_yield_v);

    if (SPECIES == "ARGON" && ION_IMPACT_SEE == "YES")
        sec_ions[0].boundary::set_i_impact_yield(&boundary::i_impact_see_yield_Ar);

    if (REFLECTION == "YES")
        sec_e.boundary::set_reflection_yield(&boundary::calc_reflection_yield_v);

    if (BACKSCATTER == "YES")
    sec_e.boundary::set_backscatter_yield(&boundary::calc_backscatter_yield_v);

   /* 
    * Select YES for steady state boundary where an electron-ion pair are inserted into 
    * the central plasma when an ion reaches a surface. Note when an electron reaches
    * the surface it is absorbed and deleted. This boundary results in a rapid
    * steady state being achieved. 
    * WARNING: BOUNDARY EFFECTS, IONIZATION AND PARTICLE INJECTION MUST BE TURNED OFF WITH THIS OPTION
    * as the number of particles is treated as a constant.
    */
    if (STEADYSTATE == "YES")
    {
        sec_e.boundary::set_ion_proportion(100.0);
        sec_ions[0].boundary::set_ion_proportion(100.0);
        sec_e.boundary::set_implement_e_boundary(&boundary::steady_e_boundary);
        sec_ions[0].boundary::set_implement_i_boundary(&boundary::steady_i_boundary);
    }
    /*Implement e_boundary and and i_boundary. This option includes secondary electrons.
     */
    else
    {
        sec_e.boundary::set_implement_e_boundary(&boundary::e_boundary);
        sec_ions[0].boundary::set_implement_i_boundary(&boundary::ion_boundary);
    }
}

/*
* Initialize graph data using user defined input. All the data is defined
* in the user input file.
* Graph data can be obtained when running in batch mode. The output is put
* in text files that can be imported into applications to generate graphs
* for post simulation output.
*/
void pic::initialize_graph_parameters(void)
{

    //graph parameters
    graphs_ions[0].graphs::set_type_in(ION_TYPE); // a label to identify ion output
    graphs_ions[0].graphs::set_vt(VTI); // ion thermal velocity
    graphs_ions[0].graphs::set_qc(CHARGE); // electron charge

    // set the number of graph steps over which to take averages over
    graphs_ions[0].graphs::set_no_graph_steps(NO_GRAPH_STEPS);

    graphs_electrons.graphs::set_type_in(ELECTRON_TYPE); // a label to identify electron output
    graphs_electrons.graphs::set_vt(VTE); // electron thermal velocity
    graphs_electrons.graphs::set_qc(ELECTRON_CHARGE); // electron charge

    // set the number of graph steps over which to take averages over
    graphs_electrons.graphs::set_no_graph_steps(NO_GRAPH_STEPS);

    //don't need a graph lable for fields
    graphs_fields.graphs::set_type_in(NO_TYPE);

    // set the number of graph steps over which to take averages over
    graphs_fields.graphs::set_no_graph_steps(NO_GRAPH_STEPS);

    //parameters for particle velocity distribution plots

    /* The number of cells around the dist_location from which velocity data is pulled
    */
    graphs_electrons.graphs::set_spread(SPREAD);
    graphs_ions[0].graphs::set_spread(SPREAD);

    /* Locations on the grid to take the particle velocity distribution */
    //for 2d define x and y location
    graphs_electrons.graphs::set_dist_location(DIST_LOCATIONX, DIST_LOCATIONY);
    graphs_electrons.graphs::set_dist_location2(DIST_LOCATION2X, DIST_LOCATION2Y);
    
    //for 2d define x and y location
    graphs_ions[0].graphs::set_dist_location(DIST_LOCATIONX, DIST_LOCATIONY);
    graphs_ions[0].graphs::set_dist_location2(DIST_LOCATION2X, DIST_LOCATION2Y);

}

/*
 * Initialize parameters needed by OpenCL and set OpenCL parameters in the particle
 * classes. The code for setting up OpenCL is courtesy of Mike Bailey from the
 * parallel programming course and he has given permission to use it.
 * 
 * For the case of using the OpenCL commands in class structure for functions
 * that a are called repeatidly in a loop the following additional resources
 * where used:
 * https://www.olcf.ornl.gov/tutorials/opencl-vector-addition/
 * https://github.com/conanwu777/particle_system/blob/master/main.cpp
 * 
*/
void pic::initialize_OpenCL(void)  
{
    const char* CL_FILE_NAME = { "pic.cl" };  // filename of the OpenCL kernel program

     // If unable to open the file holding the kernel code then exit. 

    FILE* fp;

    fp = fopen(CL_FILE_NAME, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Cannot open OpenCL source file '%s'\n", CL_FILE_NAME);
        printf("Unable to open OpenCL kernel file\n");
        exit(2);
    }
    
    // Get the platform id. 
    cl_platform_id platform;
    status = clGetPlatformIDs(1, &platform, NULL);
    if (status != CL_SUCCESS)
        fprintf(stderr, "clGetPlatformIDs failed (2)\n");

    // Get the device id and exit if unable to locate a GPU device
    status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
    if (status != CL_SUCCESS) {
        fprintf(stderr, "clGetDeviceIDs failed (2)\n");
        printf("Could not find GPU device\n");
        exit(2);
    }

    // 3. create an opencl context:

    // CONTEXT
    context = clCreateContext(NULL, 1, &device, NULL, NULL, &status);
    if (status != CL_SUCCESS)
        fprintf(stderr, "clCreateContext failed\n");

    // 4. create an opencl command queue:
    cmdQueue = clCreateCommandQueue(context, device, 0, &status);
    if (status != CL_SUCCESS)
        fprintf(stderr, "clCreateCommandQueue failed\n");
    

    
    // 7. read the kernel code from a file:
    // Courtesy of Mike Bailey from CS475 course
    fseek(fp, 0, SEEK_END);
    size_t fileSize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char* clProgramText = new char[fileSize + 1];	// leave room for '\0'
    size_t n = fread(clProgramText, 1, fileSize, fp);
    clProgramText[fileSize] = '\0';
    fclose(fp);
    if (n != fileSize)
        fprintf(stderr, "Expected to read %d bytes read from '%s' -- actually read %d.\n", fileSize, CL_FILE_NAME, n);

    // create the text for the kernel program:
    // Courtesy of Mike Bailey from CS475 course
    char* strings[1];
    strings[0] = clProgramText;
    program = clCreateProgramWithSource(context, 1, (const char**)strings, NULL, &status);
    if (status != CL_SUCCESS)
        fprintf(stderr, "clCreateProgramWithSource CC failed\n");
    delete[] clProgramText;

    // 8. compile and link the kernel code:
    // Courtesy of Mike Bailey from CS475 course
    const char* options = { "" };
    status = clBuildProgram(program, 1, &device, options, NULL, NULL);
    if (status != CL_SUCCESS)
    {
        size_t size;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &size);
        cl_char* log = new cl_char[size];
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, size, log, NULL);
        fprintf(stderr, "clBuildProgram failed:\n%s\n", log);
        delete[] log;
    }

    // 9. create the kernel objects:  

    // replaces the collect charge function in the particle class
    kernelCC = clCreateKernel(program, "collectCharge", &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateKernel failed\n");
        cout << "Status = " << status << endl;
    }
    
    // replaces the calculate particle electric field in the particle class
    kernelE = clCreateKernel(program, "get_particle_electric_field", &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateKernel failed\n");
        cout << "Status = " << status << endl;
    }

    // replaces the move particles function in the particle class
    kernelV = clCreateKernel(program, "move_particles", &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateKernel failed\n");
        cout << "Status = " << status << endl;
    }

    // replaces the sum average velocities function in the graphs class
    kernelVg = clCreateKernel(program, "sum_av_velocities", &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateKernel failed\n");
        cout << "Status = " << status << endl;
    }

    // replaces the sum average dists in the graphs class
    kernelVd = clCreateKernel(program, "sum_av_dist", &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateKernel failed\n");
        cout << "Status = " << status << endl;
    }

    // Since the grid size is constant throughout a simulation memory for grid
    // sized arrays that are to be passed to/from the GPU is reserved

    // size of the area required for passing grid arrays
    size_t qSize = Y_NO_OF_CELLS * X_NO_OF_CELLS * sizeof(cl_double);

    // will hold the particle charge collected at each grid point  (particles) 
    dQ = clCreateBuffer(context, CL_MEM_READ_WRITE, qSize, NULL, &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateBuffer failed (2)\n");
        cout << "Status = " << status << endl;
    }

    // size of the area required for passing electric field grid arrays
    size_t EgridSize = (X_NO_OF_CELLS + 1) * (Y_NO_OF_CELLS + 1) * sizeof(cl_double);

    // x component of the electric field at the grid (particles)
    dEx = clCreateBuffer(context, CL_MEM_READ_ONLY, EgridSize, NULL, &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateBuffer failed (4)\n");
        cout << "Status = " << status << endl;
    }

    // y component of the electric field at the grid (particles)
    dEy = clCreateBuffer(context, CL_MEM_READ_ONLY, EgridSize, NULL, &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateBuffer failed (5)\n");
        cout << "Status = " << status << endl;
    }

    // reserve memory for graph quantities to be passed to/from the gpu

    // one dimension of the grid size
    size_t gridSize = X_NO_OF_CELLS * sizeof(cl_double);

    // number of particles per cell where y = length/2 (graphs)
    dCountpX = clCreateBuffer(context, CL_MEM_READ_WRITE, gridSize, NULL, &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateBuffer failed (3)\n");
        cout << "Status = " << status << endl;
    }

    // number of particles per cell where x = length/2 (graphs)
    dCountpY = clCreateBuffer(context, CL_MEM_READ_WRITE, gridSize, NULL, &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateBuffer failed (4)\n");
        cout << "Status = " << status << endl;
    }

    //holds temp values for summing average velocities (graphs)
    dv_grid_tempX = clCreateBuffer(context, CL_MEM_READ_WRITE, gridSize, NULL, &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateBuffer failed (5)\n");
        cout << "Status = " << status << endl;
    }

    //holds temp values for summing average velocities (graphs)
    dv_grid_tempY = clCreateBuffer(context, CL_MEM_READ_WRITE, gridSize, NULL, &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateBuffer failed (6)\n");
        cout << "Status = " << status << endl;
    }

    size_t binSize = (NO_BINS + 1) * sizeof(cl_double);

    //holds temp values for summing average vx distribution (graphs)
    dtemp_distX1 = clCreateBuffer(context, CL_MEM_READ_WRITE, binSize, NULL, &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateBuffer failed (3)\n");
        cout << "Status = " << status << endl;
    }

    //holds temp values for summing average vx distribution (graphs)
    dtemp_distX2 = clCreateBuffer(context, CL_MEM_READ_WRITE, binSize, NULL, &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateBuffer failed (4)\n");
        cout << "Status = " << status << endl;
    }

    //holds temp values for summing average vy distribution (graphs)
    dtemp_distY1 = clCreateBuffer(context, CL_MEM_READ_WRITE, binSize, NULL, &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateBuffer failed (5)\n");
        cout << "Status = " << status << endl;
    }

    //holds temp values for summing average vy distribution (graphs)
    dtemp_distY2 = clCreateBuffer(context, CL_MEM_READ_WRITE, binSize, NULL, &status);
    if (status != CL_SUCCESS)
    {
        fprintf(stderr, "clCreateBuffer failed (6)\n");
        cout << "Status = " << status << endl;
    }

    // if steady state boundary is used then particle arrays don't change size
    // so a consistent amount of memory can be reserved at the beginning of the
    // simulation
    if (STEADYSTATE == "YES")
    {
        // size of the memory areas to reserve on the GPU
        size_t particleSizeVec = INIT_PARTICLE_NO * sizeof(particleVec);
        size_t fieldParticleSize = INIT_PARTICLE_NO * sizeof(cl_double);

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
    

    /* 
    * Pass the above obtained openCL  program variables to the particles class
    * where these will be utilized by functions to perform vectorized
    * computations.
    */

    electrons.particle::init_openCL(context, cmdQueue,kernelCC, kernelE, kernelV,dQ, dEx, dEy );
    graphs_electrons.graphs::init_openCL(context,cmdQueue,kernelVg, kernelVd, dCountpX,
               dCountpY, dv_grid_tempX, dv_grid_tempY, dtemp_distX1, dtemp_distX2,dtemp_distY1, dtemp_distY2);
               

   if (STEADYSTATE == "YES")
   {
       electrons.particle::init_openCL_steady(dP,dExp,dEyp);
       graphs_electrons.graphs::init_openCL_steady(dP);
   }
    

    for (int i = 0; i < NO_ION_SPECIES; i++) {
        ions[i].particle::init_openCL(context, cmdQueue,kernelCC, kernelE, kernelV,dQ, dEx, dEy );

        graphs_ions[i].graphs::init_openCL(context,cmdQueue,kernelVg, kernelVd, dCountpX,
               dCountpY, dv_grid_tempX, dv_grid_tempY, dtemp_distX1, dtemp_distX2,dtemp_distY1, dtemp_distY2);

        if (STEADYSTATE == "YES")
        {
            ions[i].particle::init_openCL_steady(dP,dExp,dEyp);
            graphs_ions[i].graphs::init_openCL_steady(dP);
        }

    }   

}

 /*
* Set the number of timesteps to run the simulation.
*/
void pic::setRuntime(int runtime)
{
    run_time = runtime;
}

/*
* Used to retrieve particle data for output plots by the graphics output
* functions.
*/

// get electrons
vector<particleVec> pic::getElectrons()
{
    return electrons.particle::get_particle();
}

// get ions
vector<particleVec> pic::getIons()
{
    return  ions[0].particle::get_particle();
}

/*
* Used to retrieve averaged graph class data for output plots rendered  
* in the graphics output functions.
*/

// electrostatic potential
double** pic::get_phiPlot(void)
{
    return graphs_fields.graphs::get_phiPlot();
}

// max value of the electrostatic potential
double pic::getPhiPlotMax(void)
{
    return graphs_fields.graphs::getPhiMax();
}

// x component of the electric field
double** pic::get_ExPlot(void)
{
    return graphs_fields.graphs::getExPlot();
}

// max value of the x component of the electric field
double pic::getExPlotMax(void)
{
    return graphs_fields.graphs::getExMax();
}

// y component of the electric field
double** pic::get_EyPlot(void)
{
    return graphs_fields.graphs::getEyPlot();
}

// max value of the y component of the electric field
double pic::getEyPlotMax(void)
{
    return graphs_fields.graphs::getEyMax();
}

// electron density
double** pic::get_eQPlot(void)
{
    return graphs_electrons.graphs::getQPlot();
}

// max value of the electron density
double pic::get_eQPlotMax(void)
{
    return graphs_electrons.graphs::getQMax();
}

// ion density
double** pic::get_ionQPlot(void)
{
    return graphs_ions[0].graphs::getQPlot();
}

// max value of the ion density
double pic::get_ionQPlotMax(void)
{
    return graphs_ions[0].graphs::getQMax();
}

/* True if want to write data to files. Default is true for
* batch mode and is always false for interactive mode.
*/
void pic::setWriteGraphDataToFiles(bool setFlag)
{
    writeGraphDataToFiles = setFlag;
}

/*
* Initialize the number of timesteps over which to take
* averaged data that is displayed in the graphics output
* window.
*/
void pic::initializeGraphsForInteractive(int interval)
{
    graphs_ions[0].graphs::set_no_graph_steps(interval);
    graphs_electrons.graphs::set_no_graph_steps(interval);
    graphs_fields.graphs::set_no_graph_steps(interval);
}

/* Pass false to turn off dumping simulation data into text files at the 
* end of the simulation. The default is true. These files can be used
* to start a new simulation where the previous one left off.
*/
void pic::set_dump(bool dump_in){
     dump=dump_in;
}

/* An alternative way to set the number of iterations that the simulation 
* steps through. 
*/
void pic::set_run_time(int run_time_in){
     run_time=run_time_in;
}

/* An alternative way to set the iteration number from which to start
* summing average data for the graphs class.
*/
void pic::set_start_graph(int start_graph_in){
     start_graph=start_graph_in;
}

/* An alternative way to set the iteration number to stop summing
* average data for the graphs class.
*/
void pic::set_end_graph(int end_graph_in){
     end_graph=end_graph_in;
}

/*
 * Implementation Notes: pic_run
 * -------------------------------
 * Implements the simulation cycle as described in Birdsall and Langdon,
 * C K Birdsall and A B Langdon,Plasma Physics via Computer Simulation, 
 * Chapter 2 and
 * C K Birdsall, IEEE Trans. Plasma Sci.,vol 19, (1991) 65:
 *
 * Integration of equations of motion, moving particles -> Collisions 
 *   ^                                                        |
 *   |                                                        v
 * Wieghting <- Integration of field equations on grid <- Weighting
 * 
 */
void pic::pic_run(void){
    int i,j,time_step,plot_elem=0;

    /* Catch any particles that might accidently get loaded outside of the
    * boundaries due to a poor choice of starting values by implementing 
    * the boundary functions.
    */  

    //cout << "before boundary before time start\n";
    for(i=0;i<NO_ION_SPECIES;i++){
    
                (sec_ions[i].*sec_ions[i].boundary::implement_i_boundary)(
                                        electrons.particle::get_particle(),
                                        ions[i].particle::get_particle());
                (sec_e.*sec_e.boundary::implement_e_boundary)(
                                        electrons.particle::get_particle(),
                                        ions[i].particle::get_particle());
    }

    //cout << "after boundary before time start\n";
	
	chrono::high_resolution_clock::time_point time_before, time_after;
    chrono::high_resolution_clock::time_point time_beforeEwrite, time_afterEwrite;
    chrono::high_resolution_clock::time_point time_beforeIonWrite, time_afterIonWrite;

    /* timestamp before function runs */
    time_before = chrono::high_resolution_clock::now();
	
	
    for(time_step=0; time_step < run_time; time_step++){

        if(haveCollisions)
        {
            for(i=0;i<NO_ION_SPECIES;i++)
            e_collision[i].e_collisions::collision(electrons.particle::get_particle(), ions[i].particle::get_particle());
        }
     
        //cout << "after e collisions\n";

        //Optional particle source during simulation, whereby particles are injected into the simulation domain
        electrons.particle::inject_more_particles();
        for(i=0;i<NO_ION_SPECIES;i++){
            ions[i].particle::inject_more_particles();
        }
    
        // Calculate the value of the grid electric field at the position of the electrons with charge-in-cloud weighting
        electrons.particle::get_particle_electric_field(grid_fields.fields::get_elec_fieldX(), grid_fields.fields::get_elec_fieldY());
        //cout << "after get electrons particle electric field\n";

        // Move the electrons in electric and (if set) a magnetic field
        electrons.particle::move_particles();
        // cout << "after electrons move particles\n";


    //***************************move ions******************************//

        // If it is time to move the ions
        if(time_step%DT_DTI == 0){

            for(i=0;i<NO_ION_SPECIES;i++){

                // if haveCollisions is set then implement collisions
                if (haveCollisions)
                {
                        ion_collision[i].ion_collisions::collision(ions[i].particle::get_particle());
                }
                
                // Calculate the grid electric field at the ions
                ions[i].particle::get_particle_electric_field(grid_fields.fields::get_elec_fieldX(), grid_fields.fields::get_elec_fieldY());

                // Move the ions in the electric and if set magnetic fields
                ions[i].particle::move_particles();
            
                // Implement the electron and ion boundaries together
                (sec_ions[i].*sec_ions[i].boundary::implement_i_boundary)(
                                        electrons.particle::get_particle(),
                                        ions[i].particle::get_particle());
                (sec_e.*sec_e.boundary::implement_e_boundary)(
                                        electrons.particle::get_particle(),
                                        ions[i].particle::get_particle());

            }
        }  
        else  // implement the electron boundary 
            for(i=0;i<NO_ION_SPECIES;i++){
                (sec_e.*sec_e.boundary::implement_e_boundary)(
                                        electrons.particle::get_particle(),
                                        ions[i].particle::get_particle());
            }
 //*******************************************************************//

        // If account for ionization the number of particles can grow large so
        // recalculate the weighting and reduce the particle number
		if (STEADYSTATE == "NO")
			electrons.particle::contain_particle_growth();

        for(i=0;i<NO_ION_SPECIES;i++){
            
            // If account for ionization the number of particles can grow large so
            // recalculate the weighting and reduce the particle number
			if (STEADYSTATE == "NO")
				ions[i].particle::contain_particle_growth();
    
            
            // Collect the charge at the grid point due to the ions
            ions[i].particle::collect_charge();
        }

        
        // Collect the charge at the grid point due to the electrons
        electrons.particle::collect_charge();
        //cout << "after electrons collect charge\n";

        // calculate the net charge density at each grid point using the above collected charge
        grid_fields.fields::grid_rho(ions,electrons.particle::get_q());
        //cout << "after grid rho and before grid potential\n";

        // calculate the potential at the grid points using the above net charge density
        grid_fields.fields::grid_potential(PHI_LEFT,PHI_RIGHT,PHI_BOTTOM,PHI_TOP, time_step, PHI_PERIOD);
        // cout << "after grid potential\n";

        // calculate the electric field from the above potential
        grid_fields.fields::grid_electric_field();

        //Collect the variation in some of the physical quantities with time by collecting
        //values at specific positions per time step. These values are output to a text
        //file at the end of the simulation. 
        if (time_step%TIME_DATA_FREQUENCY == 0){
            time_plot[plot_elem]=time_step;
            no_e[plot_elem]=electrons.particle::particle_size();
            phi_left[plot_elem] = grid_fields.get_phi(POINT1Y,POINT1X);
            phi_mid[plot_elem] = grid_fields.get_phi(POINT2Y,POINT2X);
            phi_right[plot_elem] = grid_fields.get_phi(POINT3Y,POINT3X);
            eden_right[plot_elem] = electrons.get_q(POINT3Y,POINT3X);

            for(i=0;i<NO_ION_SPECIES;i++){
                no_i[i][plot_elem]=ions[i].particle::particle_size();
                iden_right[i][plot_elem] = ions[i].get_q(POINT3Y,POINT3X);
            }
            plot_elem++;
        }

        /* Print particle number and time-step to screen during simulation */
        cout<<"time is "<<time_step<<" no of electrons= "<<electrons.particle::particle_size();
        for(i=0;i<NO_ION_SPECIES;i++)
        cout<<" no ions of species"<<i+1<<"  "<<ions[i].particle::particle_size();
        cout<<endl;

        
        if(time_step >= start_graph && time_step< end_graph){

            graphs_fields.graphs::sum_fields(grid_fields.fields::get_rho(),grid_fields.fields::get_phi(),
            grid_fields.fields::get_elec_fieldX(), grid_fields.fields::get_elec_fieldY());

            
            graphs_electrons.graphs::sum_density(electrons.particle::get_q());
            graphs_electrons.graphs::sum_av_velocities(electrons.particle::get_particle());
            graphs_electrons.graphs::sum_av_dist(electrons.particle::get_particle());
            
            for(i=0;i<NO_ION_SPECIES;i++){
         
            graphs_ions[i].graphs::sum_density(ions[i].particle::get_q());
            graphs_ions[i].graphs::sum_av_velocities(ions[i].particle::get_particle());
            graphs_ions[i].graphs::sum_av_dist(ions[i].particle::get_particle());
           
            }
        }

    if(time_step == end_graph-1){

        graphs_fields.graphs::average_fields();

        graphs_electrons.graphs::average_velocities();
        graphs_electrons.graphs::average_density();
        graphs_electrons.graphs::average_dist1();
        graphs_electrons.graphs::average_dist2();

        if (writeGraphDataToFiles)
        {
            graphs_fields.graphs::write_fields_to_file();
            graphs_fields.graphs::write_rho_to_file();
            graphs_fields.graphs::write_phi_to_file();
            graphs_fields.graphs::write_EgridX_to_file();
            graphs_fields.graphs::write_EgridY_to_file();

            graphs_electrons.graphs::write_density_to_file();
            graphs_electrons.graphs::write_v_to_file();
            graphs_electrons.graphs::write_dist1_to_file();
            graphs_electrons.graphs::write_dist2_to_file();
        }
        

        for(i=0;i<NO_ION_SPECIES;i++){
           graphs_ions[i].graphs::average_velocities();
           graphs_ions[i].graphs::average_density();
           graphs_ions[i].graphs::average_dist1();
           graphs_ions[i].graphs::average_dist2();

           if (writeGraphDataToFiles)
           {
               graphs_ions[i].graphs::write_density_to_file();
               graphs_ions[i].graphs::write_v_to_file();
               graphs_ions[i].graphs::write_dist1_to_file();
               graphs_ions[i].graphs::write_dist2_to_file();
           }
      
       }
    }


} /*timestep loop*/

    /* timestamp after function runs */
    time_after = chrono::high_resolution_clock::now();

    chrono::minutes interval = chrono::duration_cast<chrono::minutes>(time_after - time_before);
    

    cout << "duration of PIC run loop = " << interval.count() << " minutes" << endl;
    

  //put time data in file 
    if (INTERACTIVE == "NO")
    {
        for (j = 0; j < SIZE_TIME_DATA; j++) {

            timedata << time_plot[j] << " " << no_e[j];
            for (i = 0; i < NO_ION_SPECIES; i++)
                timedata << " " << no_i[i][j];
            timedata << " " << phi_left[j] << " " << phi_mid[j] << " " << phi_right[j] << " " << eden_right[j];
            for (i = 0; i < NO_ION_SPECIES; i++)
                timedata << " " << iden_right[i][j];
            timedata << endl;
        }
    }
   
   
   cout<<"no electron reductions= "<<electrons.particle::count_reduce<<endl;
   cout<<"no ion reductions= "<<ions[0].particle::count_reduce<<endl;
}

/*
* Dump the current state of each particle type to a separate output file.
* The contents of the struct particleVec for each particle is output.
*/
void pic::dump_data()
{

    cout << "Dumping Data... \n";
    
    electrons.particle::dump_data();
    for (int i = 0; i < NO_ION_SPECIES; i++)
        ions[i].particle::dump_data();

    cout << "Done\n";  
}

pic::~pic(){
    int i;

    // release memory objects for openCL
    clReleaseMemObject(dQ);
    clReleaseMemObject(dEx);
    clReleaseMemObject(dEy);
    clReleaseMemObject(dCountpX);
    clReleaseMemObject(dCountpY);
    clReleaseMemObject(dv_grid_tempX);
    clReleaseMemObject(dv_grid_tempY);
    clReleaseMemObject(dtemp_distX1);
    clReleaseMemObject(dtemp_distX2);
    clReleaseMemObject(dtemp_distY1);
    clReleaseMemObject(dtemp_distY2);

    // release openCL objects
    clReleaseProgram(program);
    clReleaseCommandQueue(cmdQueue);
    clReleaseContext(context);
    clReleaseKernel(kernelCC);
    clReleaseKernel(kernelE);  
    clReleaseKernel(kernelV);
    clReleaseKernel(kernelVg);
    clReleaseKernel(kernelVd);

    delete [] ions;
    delete [] e_collision;
    delete [] ion_collision;
    delete [] graphs_ions;
    delete [] time_plot;
    delete [] no_e;
    delete [] phi_mid;
    delete [] phi_right;
    delete [] phi_left;
    delete [] eden_right;
    for(i=NO_ION_SPECIES;i>0;i--){
        delete[] no_i[i-1];
        delete[] iden_right[i-1];
    }
    delete[] no_i;
    delete[] iden_right;
    delete[] sec_ions;
}

