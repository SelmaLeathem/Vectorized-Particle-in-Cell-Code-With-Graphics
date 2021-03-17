/* 
 * File: main.cpp
 * ----------------
 * 
 * 
 * Date: 2/25/2021
 * 
 * Runs a two-dimensional particle-in-cell (PIC) simulation interactively or in.
 * batch mode. More information about the simulation is provided in the instructions.  
 * In summary a plasma between four boundaries of arbitrary electrostatic potential
 * is modelled. Optionally, secondary emmision from the surfaces and collisions
 * can be accounted for by setting the appropriate flags in the parameters.h file.
 * 
 * User input to the simulation is done via the paramters.h file.
 * 
 * The overall runtime of the simulation run for batch mode is output.
 *
 */


#include <time.h>
#include "constants.h"
#include "parameters.h"
#include "pic.h"
#include <chrono>   /* use the high resolution clock */
#include <iostream>
#include "graphics.h"

using namespace std;


int main(int argc, char* argv[]) {

  /* Set the random seed. */
   srand(RAND_SEED);

	/* If running in interactive mode then the simulation is started and run from within
	 * the graphics.cpp file. The following initializes OpenGL and starts the glut
	 * main loop.
	  */
   if (INTERACTIVE == "YES")
   {
	   glutInit(&argc, argv);
	   InitGraphics();  // Initialize the graphics
	   InitLists(); // store static display objects in gpu memory
	   InitPIC(); // initialize the PIC simulation
	   InitParameters(); //initalize parameters used when displaying graphs
	   InitGlui(); // initialize the user interface
	   glutMainLoop(); // run the graphics main loop
   }
   else  // run in batch mode
   {
	   /* time the simulation in batch mode */
	   chrono::high_resolution_clock::time_point time_before, time_after;

	   /* Need to initialize restart to true via the constructor in order to
	   * start a new simulation from dump files.
	   *
	   * e.g. pic pic1(true);
	   */
	   pic pic1(RESTART);

	   /* set the total runtime of the simulation */
	   pic1.pic::setRuntime(TOTAL_RUN_TIME);

	   /* allow graph data to be dumped into text files */
	   pic1.pic::setWriteGraphDataToFiles(true);

	   /* timestamp before function runs */
	   time_before = chrono::high_resolution_clock::now();

		/* run the simulation */
	   pic1.pic::pic_run();

	   /* timestamp after function runs */
	   time_after = chrono::high_resolution_clock::now();
	   chrono::minutes interval = chrono::duration_cast<chrono::minutes>(time_after - time_before);

	   cout << "duration of run = " << interval.count() << " minutes" << endl;
	   
		/* Dump the current state of the simulation. These files can be used to
		 * start a future simulation from this point.
		 */
	   pic1.pic::dump_data();
   }

   return 0;
}
