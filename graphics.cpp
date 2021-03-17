/*
 * Date: 2/26/2021
 * Description: Graphics function and variable definitions for the purpose of
 *              implementing a glui user interface with simulation graph output
 *              using OpenGL for the graphics.
 * 
 *              References for GLUI:
 *              Mike Bailey's class code (CS475)
 *              http://www-h.eng.cam.ac.uk/help/tpl/graphics/using_glui.html
 *              user manual: https://github.com/libglui/glui
 *       
 *              Note: Indicated parts of the code are from files given in Mike
 *              Bailey's graphics and parallel programming classes and have been
 *              replicated here with his permission.
*/

#include "graphics.h"
#include <iostream>
#include <fstream>

// title of these windows:
const char *WINDOWTITLE = { "Particle-In-Cell Graphics Output" };
const char *GLUITITLE   = { "PIC 2D" };

// define glui true/false constants
const int GLUITRUE  = { true  };
const int GLUIFALSE = { false };


// window size in pixels
const int INIT_WINDOW_WIDTH = { 700 };
const int INIT_WINDOW_HEIGHT = { 700 };

// Rotation angle increment
const float ANGFACT = { 1. };
// Minimum size an object being rendered can by scaled by
const float MINSCALE = { 0.001f };

// Mouse buttons
const int LEFT   = { 4 };
const int MIDDLE = { 2 };
const int RIGHT  = { 1 };

// Image projection options: orthographic or perspective
enum Projections
{
	ORTHO,
	PERSP
};

// Indicates whether or not to dump data
enum DumpData
{
	DO_NOT_DUMP,
	DUMP
};

// User Interface button options
enum ButtonVals
{
	START,
	PAUSE,
	QUIT
};

// Indicates which graph a radio button links to for plotting for the user interface
enum WhichGraph
{
	PARTICLE_XY,
	ELECTRON_VX_VY,
	ION_VX_VY,
	PHI,
	EFIELD_X,
	EFIELD_Y,
	E_DENSITY,
	ION_DENSITY
};

// Which radio button the user selected
int whichGraphToPlot;

// Background color of the main window
const float BACKCOLOR[ ] = { 0., 0., 0., 0. };

// Properies of the main axes
const GLfloat AXES_COLOR[ ] = { 1.0, 1.0, 1.0 };
const GLfloat AXES_WIDTH   = { 2. };
const GLfloat AXES_LENGTH = 2.0;

// Position of the graph titles in the main window
const GLfloat plotTitlePosition[3] = {0.2f, 2.8f, 0.0f };

// Used for scaling the particle and field graphs within the main window
float boxLength = 4.0;
float scaleFields = 1.0; 

// Position of field z-axis label in field graph plots
const GLfloat fieldPlotLabel[3] = { 0.2f, scaleFields * 2.0, 0.f };

// Origin of the field graphs
const GLfloat fieldAxisPlotPosition[3] = { -scaleFields / 2.0, -scaleFields / 2.0 , -scaleFields / 2.0 };

// Initial angle of the field graphs
const GLfloat fieldAxisPlotAngle = -45.0;

// Position of the x-y labels in particle graph plots
const GLfloat particlePlotLabelX[3] = { boxLength/2.0 + 0.1, 0.f, 0.f };
const GLfloat particlePlotLabelY[3] = {-0.4f, boxLength / 2.0 + 0.1, 0.f };
const GLfloat particlePlotLabelXm[3] = {-boxLength / 2.0 - 1.5, 0.f, 0.f };
const GLfloat particlePlotLabelYm[3] = { -0.5f, -boxLength / 2.0 - 0.3, 0.f };



// non-constant global variables:
// These are taken from Mike Bailey's graphics class files
int	ActiveButton;			// current button that is down
GLuint	AxesList;			// list to hold the axes
int	AxesOn;					// ON or OFF
int	GluiWindow;				// the glut id for the glui window
int	MainWindow;				// window id for main graphics window
int	Paused;					// switch that turns pause on/off
GLfloat	RotMatrix[4][4];	// set by glui rotation widget
float	Scale, Scale2;		// scaling factors
int	WhichProjection;		// ORTHO or PERSP
int	Xmouse, Ymouse;			// mouse values
float	Xrot, Yrot;			// rotation angles in degrees
float	TransXYZ[3];		// set by glui translation widgets

// An instance of a glui user interface window
GLUI* GluiRight;			

// An instance of the pic (particle-in-cell) class
pic pic1(RESTART);

// A handle to refer to the frame that is drawn around particle plots
GLuint	boxScaleList;		

// Size of the grid, assuming a square grid
const int GRIDSIZE = X_NO_OF_CELLS;

// id of the selected radio button
int radio_button_id = 0;

// number of timesteps into the simulations
int num_time_steps = 0;

// the simulation is run in run_interval lots of time before returning
int run_interval = PLOT_UPDATE_INTERVAL;

// flag to indicate whether to dump data upon quit
int dumpData = 0;

/*
* Through glut callbacks gets the latest data for the graphs by calling the 
* simulation's pic class run function. This function is periodically called
* by glut Idle. 
*/
void Animate( )
{
	glutSetWindow( MainWindow );
	glFinish( );

	// Run the simulation for another x number of timesteps. Note that this
	// updates the quantities being plotted in the display.
	pic1.pic::pic_run();
	num_time_steps += run_interval;
	printf("number of timesteps = %d\n", num_time_steps);

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}

/* Glui callback for the radio buttons which are used to decide which 
 * graph to display next.
 */
void GraphRadioGroup(int id)
{
	switch (radio_button_id)
	{
	case PARTICLE_XY:
		whichGraphToPlot = PARTICLE_XY;
		resetRotation();
		break;
	case ELECTRON_VX_VY:
		whichGraphToPlot = ELECTRON_VX_VY;
		resetRotation();
		break;
	case ION_VX_VY:
		whichGraphToPlot = ION_VX_VY;
		resetRotation();
		break;
	case PHI:
		whichGraphToPlot = PHI;
		break;
	case EFIELD_X:
		whichGraphToPlot = EFIELD_X;
		break;
	case EFIELD_Y:
		whichGraphToPlot = EFIELD_Y;
		break;
	case E_DENSITY:
		whichGraphToPlot = E_DENSITY;
		break;
	case ION_DENSITY:
		whichGraphToPlot = ION_DENSITY;
		break;
	default:
		break;
	}
}


/* What happens when a user presses a button on the glui interface.
 * Code from Mike Bailey's graphics class files.
 */
void Buttons( int id )
{
	switch( id )
	{
		case START:
			GLUI_Master.set_glutIdleFunc( Animate );
			break;
		case PAUSE:
			Paused = ! Paused;
			if( Paused )
				GLUI_Master.set_glutIdleFunc( NULL );
			else
				GLUI_Master.set_glutIdleFunc( Animate );
			break;
		case QUIT:
			Quit( );
			break;

		default:
			fprintf( stderr, "Don't know what to do with Button ID %d\n", id );
	}

}

/* Draws field graph quantities to the main window assuming two dimensional data
 * on a square grid. A series of openGL quadrilaterals are drawn over the grid.
 * The color of the quadrilateral indicates the z value of the field quantity at
 * that location on the grid.
 */
void drawField(double **field, double fieldMax, int gridStart, char *strLabel, char *strTitle)
{
	char strGraph[128];
	float z0, z1, z2, z3; // z values at the four corners of the quadrilateral
	float zAv; // average of the z values used to pick a color for the quadrilateral
	float alpha = 0.f; // used for mixing boundary colors
	// z is scaled to fall between 0 and 1. The range divider value indicates
	// when a different color is used.
	float rangeDivider = 0.2; 
	float x, y;  // values of the x and y coordinates
	float dxg = 2.0 / GRIDSIZE;  // for scaling the grid to fall between -1 and 1

	// for making the quadrilateral slighly smaller than the cell size to distinguish
	// contours/colors more easily
	float dc = 1.0 / (GRIDSIZE-1);
	float dc_offset = dc*0.2; 

	// use glm::vec3 types to hold the following colors in r,g,b formats
	glm::vec3 yellow = glm::vec3(1.000, 1.0, 0.0);
	glm::vec3 orange = glm::vec3(1.0, 0.396, 0.0);
	glm::vec3 purple = glm::vec3(0.396, 0.0, 1.0);
	glm::vec3 red = glm::vec3(1.0, 0.0, 0.0);
	glm::vec3 green = glm::vec3(0.0, 1.0, 0.0);
	glm::vec3 blue = glm::vec3(0.0, 0.0, 1.0);

	// holds a mixture of the two colors on either side of a color range
	// where the value of z decides how much of each is mixed in
	glm::vec3 mixedColor;

	// maximum value of the field to scale everything to lie between -1 and 1
	fieldMax = abs(fieldMax);
	
	if (fieldMax < 1.0)
	fieldMax = 1.0;

	// draw the axies
	glCallList(AxesList);

	// draw the field
	glPushMatrix();

	// translate to align with the axes
	glTranslatef(fieldAxisPlotPosition[0], fieldAxisPlotPosition[1], fieldAxisPlotPosition[2]);

	// rotate for a better starting view
	glRotatef(fieldAxisPlotAngle, 0.0, 1.0, 0.0);

	// print the graph title
	glColor3f(0.8, 0.9, 0.9);
	DoRasterString(plotTitlePosition[0], plotTitlePosition[1], plotTitlePosition[2], strTitle);

	// print the graph label with the latest value of the maximum field
	glColor3f(0.99, 0.99, 0.99);
	sprintf(strGraph, strLabel, (float)fieldMax);
	DoRasterString(fieldPlotLabel[0], fieldPlotLabel[1], fieldPlotLabel[2], strGraph);

	// draw the quadrilaterals
	glBegin(GL_QUADS);
		for (int j = gridStart; j < GRIDSIZE-gridStart-1; j++)
		{
			for (int i = gridStart; i < GRIDSIZE-gridStart-1; i++)
			{
				// z values of the four corners of the quadrilateral
				z0 = field[j + 1][i] / fieldMax ;
				z1 = field[j + 1][i + 1] / fieldMax;
				z2 = field[j][i + 1] / fieldMax;
				z3 = field[j][i] / fieldMax;

				// the average value is used to decide the color of the square
				zAv = abs((z0 + z1 + z2 + z3) / 4.0);
				
				// zAv needs to be less than one
				alpha = zAv < 1.0 ? zAv : 1.0;

				// if the z value falls within a certain range then get a
				// color for the quadrilateral in that range using the
				// mix function which calculates: 
				//               alpha*color1 + (1-alpha)*color2
				if (zAv < rangeDivider)
				{
					alpha /= rangeDivider;
					mixedColor = glm::mix(red, orange, alpha);
					
				}	
				else if (zAv < 2.0 * rangeDivider)
				{
					alpha = (alpha - rangeDivider) / rangeDivider;
					mixedColor = glm::mix(orange, yellow, alpha);
					
				}
				else if (zAv < 3.0 * rangeDivider)
				{
					alpha = (alpha - 2.0*rangeDivider) / rangeDivider;
					mixedColor = glm::mix(yellow, green, alpha);
					
				}
				else if (zAv < 4.0 * rangeDivider)
				{
					alpha = (alpha - 3.0*rangeDivider) / rangeDivider;
					mixedColor = glm::mix(green, blue, alpha);
					
				}
				else
				{
					alpha = (alpha - 4.0*rangeDivider) / rangeDivider;
					mixedColor = glm::mix(blue, purple, alpha);
					
				}
					
				// select the color
				glColor3f(mixedColor.x, mixedColor.y, mixedColor.z);

				// Get the x and y values by scaling to ensure the graph
				// fits within -1 to 1. The offset results in a slightly
				// smaller quadrilateral making contours easier to see
				x = i * dxg * scaleFields + dc_offset; 
				y = (j + 1) * dxg * scaleFields + dc_offset; 
				glVertex3f(x, z0 * scaleFields, y);

				x = (i + 1) * dxg * scaleFields - dc_offset;
				y = (j + 1) * dxg * scaleFields + dc_offset;
				glVertex3f(x, z1 * scaleFields, y);

				x = (i + 1) * dxg * scaleFields - dc_offset;
				y = j * dxg * scaleFields - dc_offset;
				glVertex3f(x, z2 * scaleFields, y);

				x = i * dxg * scaleFields + dc_offset;
				y = j * dxg * scaleFields - dc_offset;
				glVertex3f(x, z3 * scaleFields, y);
				
			}
		}
		glEnd();
	glPopMatrix();
}

/* Resets the glut and glui rotation matrices back to their default positions */
void resetRotation(void)
{
	Xrot = Yrot = 0.;
	TransXYZ[0] = TransXYZ[1] = TransXYZ[2] = 0.;

	RotMatrix[0][1] = RotMatrix[0][2] = RotMatrix[0][3] = 0.;
	RotMatrix[1][0] = RotMatrix[1][2] = RotMatrix[1][3] = 0.;
	RotMatrix[2][0] = RotMatrix[2][1] = RotMatrix[2][3] = 0.;
	RotMatrix[3][0] = RotMatrix[3][1] = RotMatrix[3][3] = 0.;
	RotMatrix[0][0] = RotMatrix[1][1] = RotMatrix[2][2] = RotMatrix[3][3] = 1.;
}

/* Draws particle phase graphs (vx-vy) to the main window */
void drawParticlePhase(vector<particleVec> particle, glm::vec3 particleColor, char * strTitle)
{
	int no_particles;  // size of particle vector
	int index, num_to_plot; // number of particles to include in plot
	int maxVx, maxVy; // maximum vx and vy values
	float x, y; // x,y positions of the particles
	float vxp, vyp; // vx, vy values of each particle
	float length = (float)XLENGTH;

	no_particles = particle.size();

	num_to_plot = NUM_PARTICLES_TO_PLOT;

	// print the plot title
	glColor3f(0.797, 0.8945, 0.996);
	DoRasterString(plotTitlePosition[0], plotTitlePosition[1], plotTitlePosition[2], strTitle);

	// draw the frame that goes around the particles
	glCallList(boxScaleList);

	// randomly select num_to_plot amount of particles to plot
	index = (int)(((double)rand() / ((double)RAND_MAX + 1.0)) * no_particles);
	
	// get maximum absolute velocity values to ensure vx, vy values fall 
	// between -1 to 1
	maxVx = abs((float)(particle[index].vx));
	maxVy = abs((float)(particle[index].vy));

	num_to_plot = (num_to_plot > no_particles) ? no_particles : num_to_plot;

	// first verify that the particles is inside the boundaries before 
	// calculating the maximum value of the velocities
	for (int i = 0; i < num_to_plot; i++)
	{
		index = (int)(((double)rand() / ((double)RAND_MAX + 1.0)) * no_particles);
		index = (num_to_plot == no_particles) ? i : index;

		// ensure plotting a particle within the boundaries
		x = (float)(particle[index].x) * boxLength / length - boxLength / 2.0;
		y = (float)(particle[index].y) * boxLength / length - boxLength / 2.0;
		
		if (abs(x) < boxLength / 2.0 && abs(y) < boxLength / 2.0)
		{
			if (maxVx < abs((float)(particle[index].vx)))
				maxVx = abs((float)(particle[index].vx));

			if (maxVy < abs((float)(particle[index].vy)))
				maxVy = abs((float)(particle[index].vy));

		}	
	}

	if (maxVx < 1.0)
		maxVx = 1.0;

	if (maxVy < 1.0)
		maxVy = 1.0;

	// Chart plot labels show the current minumum and maximum velocity values on the border

	char strEvx[128];
	glColor3f(0.99, 0.99, 0.99);
	sprintf(strEvx, "%3.2e m/s", (float)maxVx);
	DoRasterString(particlePlotLabelX[0], particlePlotLabelX[1]+ 0.3, particlePlotLabelX[2], (char *)"Max");
	DoRasterString(particlePlotLabelX[0], particlePlotLabelX[1], particlePlotLabelX[2], strEvx);

	char strEvxm[128];
	glColor3f(0.99, 0.99, 0.99);
	sprintf(strEvxm, "%3.2e m/s", -(float)maxVx);
	DoRasterString(particlePlotLabelXm[0]+0.2, particlePlotLabelXm[1] + 0.3, particlePlotLabelXm[2], (char *)"Min");
	DoRasterString(particlePlotLabelXm[0], particlePlotLabelXm[1], particlePlotLabelXm[2], strEvxm);

	char strEvy[128];
	glColor3f(0.99, 0.99, 0.99);
	sprintf(strEvy, "Max %3.2e m/s", (float)maxVy);
	DoRasterString(particlePlotLabelY[0]-0.5, particlePlotLabelY[1], particlePlotLabelY[2], strEvy);

	char strEvym[128];
	glColor3f(0.99, 0.99, 0.99);
	sprintf(strEvym, "Min %3.2e m/s", -(float)maxVy);
	DoRasterString(particlePlotLabelYm[0]-0.5, particlePlotLabelYm[1], particlePlotLabelYm[2], strEvym);

	// plot the actual particle velocities. It doesn't matter if different random particles
	// are selected from those that decided the maximum velocity values as if a 
	// particle velocity falls outside of this range then it is not plotted and it is
	// anticpated most particles will fall within range.		
	glPushMatrix();

	glPointSize(2.0f);
	
	glColor3f(particleColor.x, particleColor.y, particleColor.z);
	glBegin(GL_POINTS);
	
	for (int i = 0; i < num_to_plot; i++)
	{
		index = (int)(((double)rand() / ((double)RAND_MAX + 1.0)) * no_particles);
		index = (num_to_plot == no_particles) ? i : index;

		// ensure plotting a particle within the boundaries
		x = (float)(particle[index].x) * boxLength / length - boxLength / 2.0;
		y = (float)(particle[index].y) * boxLength / length - boxLength / 2.0;
		if (abs(x) < boxLength / 2.0 && abs(y) < boxLength / 2.0)
		{
			// scale the vertex coordinates to be between -1 and 1
			vxp = (float)(particle[index].vx) * boxLength / (2.0 * maxVx);
			vyp = (float)(particle[index].vy) * boxLength / (2.0 * maxVy);
			if (abs(vxp) < boxLength / 2.0 && abs(vyp) < boxLength / 2.0)
				glVertex3f(vxp, vyp, 0.0);
		}	
	}

	glEnd();
	glPopMatrix();
}

/* Draws particle X-Y graphs to the main window */
void draw_particleXY(vector<particleVec> electrons, vector<particleVec> ions)
{
	float length = (float)XLENGTH;  // length of grid side
	int no_e_particles;  // vector size
	int no_i_particles;
	int index, num_to_plot;  // number of particles to include in plot
	float x, y;
	// vec3 types holding the r,g,b values for the indicated colors
	glm::vec3 hotPink = glm::vec3(1.0, 0.114, 0.557);
	glm::vec3 lime = glm::vec3(0.557, 1.0, 0.114);

	num_to_plot = NUM_PARTICLES_TO_PLOT;  // number of particles to include

	no_e_particles = electrons.size();
	no_i_particles = ions.size();

	// plot title
	glColor3f(0.797, 0.8945, 0.996);
	DoRasterString(plotTitlePosition[0], plotTitlePosition[1], plotTitlePosition[2], (char*)"Particle X-Y Plot");

	// draw the frame that goes around the particles
	glCallList(boxScaleList);
		
	// chart plot labels show the dimensions of the domain

	char strLx[128];
	glColor3f(0.99, 0.99, 0.99);
	sprintf(strLx, "X = %3.2e m", (float)length);
	DoRasterString(particlePlotLabelX[0], particlePlotLabelX[1], particlePlotLabelX[2], strLx);

	char strLxm[128];
	glColor3f(0.99, 0.99, 0.99);
	sprintf(strLxm, "X = %2.1f m", (float)0.0);
	DoRasterString(particlePlotLabelXm[0]+0.4, particlePlotLabelXm[1], particlePlotLabelXm[2], strLxm);

	char strLy[128];
	glColor3f(0.99, 0.99, 0.99);
	sprintf(strLy, "Y = %3.2e m", (float)length);
	DoRasterString(particlePlotLabelY[0], particlePlotLabelY[1], particlePlotLabelY[2], strLy);

	char strLym[128];
	glColor3f(0.99, 0.99, 0.99);
	sprintf(strLym, "Y = %2.1f m", (float)0.0);
	DoRasterString(particlePlotLabelYm[0], particlePlotLabelYm[1], particlePlotLabelYm[2], strLym);

	glPushMatrix();

	// plot the electrons

	//point properties
	glPointSize(2.0f);
	glColor3f(lime.x, lime.y, lime.z);

	glBegin(GL_POINTS);

	num_to_plot = (num_to_plot > no_e_particles) ? no_e_particles : num_to_plot;
	
	for (int i = 0; i < num_to_plot; i++)
	{
		// randomly select num_to_plot amount of particles to plot
		index = (int)(((double)rand() / ((double)RAND_MAX + 1.0)) * no_e_particles);
		index = (num_to_plot == no_e_particles) ? i : index;
		
		// ensure plotting a particle within the boundaries
		x = (float)(electrons[index].x) * boxLength / length - boxLength / 2.0;
		y = (float)(electrons[index].y) * boxLength / length - boxLength / 2.0;
		if (abs(x) < boxLength / 2.0 && abs(y) < boxLength / 2.0)
		{
			glVertex3f(x, y, 0.0);
			
		}		
	}	

	// plot the ions
	num_to_plot = NUM_PARTICLES_TO_PLOT;  // reset number of particles to include

	// point properties
	glPointSize(2.0f);
	glColor3f(hotPink.x, hotPink.y, hotPink.z);

	num_to_plot = (num_to_plot > no_i_particles) ? no_i_particles : num_to_plot;
	
	for (int i = 0; i < num_to_plot; i++)
	{
		// randomly select num_to_plot amount of particles to plot
		index = (int)(((double)rand() / ((double)RAND_MAX + 1.0)) * no_i_particles);
		index = (num_to_plot == no_i_particles) ? i : index;

		// ensure plotting a particle within the boundaries
		x = (float)(ions[index].x) * boxLength / length - boxLength / 2.0;
		y = (float)(ions[index].y) * boxLength / length - boxLength / 2.0;
		if (abs(x) < boxLength / 2.0 && abs(y) < boxLength / 2.0)
			glVertex3f(x, y, 0.0);
	}
	glEnd();
	glPopMatrix();
}

/* Where everything is drawn to the main graphics window. A lot of the
 * openGL window setup code is from Mike Bailey's graphics class files
 */
void Display( )
{
	int gridStart; // starting and ending grid position to plot
	// arrays and vectors holding quantites being plotted
	vector<particleVec> electrons; 
	vector<particleVec> ions; 
	double** phi;  // electrostatic potential
	double phiMax;  // maximum  electrostatic potential
	double** EfieldX; // electric field in x
	double** EfieldY; // electric field in y
	double ExMax; // maximum  electric field in x
	double EyMax; // maximum  electric field in y
	double** eQdensity; // electron density
	double** ionQdensity; // ion density
	double eQMax;  // maximum electron density
	double ionQMax; // maximum ion density

	// vec3 types holding the r,g,b values for the indicated colors
	glm::vec3 hotPink = glm::vec3(1.0, 0.114, 0.557);
	glm::vec3 lime = glm::vec3(0.557, 1.0, 0.114);

	// Setting up the main openGL window
	glutSetWindow( MainWindow );
	glDrawBuffer( GL_BACK );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glEnable( GL_DEPTH_TEST );
	glShadeModel( GL_FLAT );
	GLsizei vx = glutGet( GLUT_WINDOW_WIDTH );
	GLsizei vy = glutGet( GLUT_WINDOW_HEIGHT );
	GLsizei v = vx < vy ? vx : vy;			// minimum dimension
	GLint xl = ( vx - v ) / 2;
	GLint yb = ( vy - v ) / 2;
	glViewport( xl, yb,  v, v );

	// set the viewing volume:
	// remember that the Z clipping  values are actually
	// given as DISTANCES IN FRONT OF THE EYE
	// USE gluOrtho2D( ) IF YOU ARE DOING 2D !

	glMatrixMode( GL_PROJECTION );
	glLoadIdentity( );
	if( WhichProjection == ORTHO )
		glOrtho(-3., 3., -3., 3., 0.1, 1000.);
	else
		gluPerspective(90., 1., 0.1, 1000.);
		
	// set up model view
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity( );

	// set the eye position, look-at position, and up-vector:
	gluLookAt(0., 0., 3.5, 0., 0., 0., 0., 1., 0.);

	// rotate the scene:
	glRotatef( (GLfloat)Yrot, 0., 1., 0. );
	glRotatef( (GLfloat)Xrot, 1., 0., 0. );
	glMultMatrixf( (const GLfloat *) RotMatrix );

	// uniformly scale the scene:
	glScalef( (GLfloat)Scale, (GLfloat)Scale, (GLfloat)Scale );
	float scale2 = 1.f + Scale2;		// because glui translation starts at 0.
	if( scale2 < MINSCALE )
		scale2 = MINSCALE;
	glScalef( (GLfloat)scale2, (GLfloat)scale2, (GLfloat)scale2 );

	// since using glScalef( ), be sure normals get unitized:
	glEnable(GL_NORMALIZE);

	// select a graph to plot based upon the value of the radio button selected by
	// the user
	switch(whichGraphToPlot)
	{
		case PARTICLE_XY:  // Plot the position of the particles in the domain
			electrons = pic1.pic::getElectrons();
			ions = pic1.pic::getIons();
			draw_particleXY(electrons, ions);
			break;
		case ELECTRON_VX_VY: // plot the electron phase vx-vy
			electrons = pic1.pic::getElectrons();
			drawParticlePhase(electrons, lime, (char*)"Electron Vx-Vy Plot");
			break;
		case ION_VX_VY:  // plot the ion phase vx-vy
			ions = pic1.pic::getIons();
			drawParticlePhase(ions, hotPink, (char*)"Ion Vx-Vy Plot");
			break;
		case PHI:  // plot the electrostatic potential
			gridStart = 0;
			phi = pic1.pic::get_phiPlot();
			phiMax = pic1.pic::getPhiPlotMax();
			drawField(phi, phiMax, gridStart,
					(char *)"Max  %3.1f Volts", (char *)"Electrostatic Potential");
			break;
		case EFIELD_X:  // plot the x component of the electric field
			gridStart = 0;
			EfieldX = pic1.pic::get_ExPlot();
			ExMax = pic1.pic::getExPlotMax();
			drawField(EfieldX, ExMax, gridStart,
				(char *)"Max  %6.1f Volts/meter", (char *)"Electric Field in X");	
			break;
		case EFIELD_Y:	// plot the y component of the electric field
			gridStart = 0;
			EfieldY = pic1.pic::get_EyPlot();
			EyMax = pic1.pic::getEyPlotMax();
			drawField(EfieldY, EyMax, gridStart,
				(char *)"Max  %6.1f Volts/meter", (char *)"Electric Field in Y");
			break;
		case E_DENSITY:  // plot the electron density
			gridStart = 1;
			eQdensity = pic1.pic::get_eQPlot();
			eQMax = pic1.pic::get_eQPlotMax();
			drawField(eQdensity, eQMax, gridStart,
			(char *)"Max  %3.2e  1/meter^3", (char *)"Electron Density");		
		break;
		case ION_DENSITY:  // plot the ion density
			gridStart = 1;
			ionQdensity = pic1.pic::get_ionQPlot();
			ionQMax = pic1.pic::get_ionQPlotMax();
			drawField(ionQdensity, ionQMax, gridStart,
			(char *)"Max  %3.2e  1/meter^3", (char *)"Ion Density");
		break;
		default:
			// do nothing
		break;
	}

	glutSwapBuffers( );
	glFlush( );
}


// Initialize the particle-in-cell code
void InitPIC(void)
{ 
	int zero = 0;

	// graph data is not written to output files after taking an average
	// since it is instead rendered to the screen
	pic1.pic::setWriteGraphDataToFiles(false);

	// how long to run the simulation before getting more data
	pic1.pic::setRuntime(PLOT_UPDATE_INTERVAL);

	// how many timesteps to take a graph average over
	pic1.pic::set_start_graph(zero);
	pic1.pic::set_end_graph(PLOT_UPDATE_INTERVAL);
	pic1.pic::initializeGraphsForInteractive(PLOT_UPDATE_INTERVAL);

}

/* Sets up the Glui user interface and its controls 
*  Some of the the code is from graphics class files.
*  References for GLUI:
*              Mike Bailey's class code (CS475)
*              http://www-h.eng.cam.ac.uk/help/tpl/graphics/using_glui.html
*              user manual: https://github.com/libglui/glui
*/
void InitGlui( )
{
	int selectedRadioButton;  // radio button selected by user
	int gluiType=0; // not used value indicates which button type is being used
	glutInitWindowPosition( INIT_WINDOW_WIDTH + 50, 0 );
	GluiRight = GLUI_Master.create_glui( (char *) GLUITITLE );

	// add a panel for the radio buttons
	GLUI_Panel* panel_select = GluiRight->add_panel((char *)"Select Graph");
	GLUI_RadioGroup* whichGraph = GluiRight->add_radiogroup_to_panel(panel_select, &radio_button_id, gluiType, (GLUI_Update_CB)GraphRadioGroup);
	GluiRight->add_radiobutton_to_group(whichGraph, (char *)"Particle X-Y");
	GluiRight->add_radiobutton_to_group(whichGraph, (char*)"Electron Vx-Vy");
	GluiRight->add_radiobutton_to_group(whichGraph, (char*)"Ion Vx-Vy");
	GluiRight->add_radiobutton_to_group(whichGraph, (char*)"Electrostatic Potential");
	GluiRight->add_radiobutton_to_group(whichGraph, (char*)"Electric Field in X");
	GluiRight->add_radiobutton_to_group(whichGraph, (char*)"Electric Field in Y");
	GluiRight->add_radiobutton_to_group(whichGraph, (char*)"Electron Charge Density");
	GluiRight->add_radiobutton_to_group(whichGraph, (char*)"Ion Charge Density");

	// add a panel that enables graph rotation with a widget
	GLUI_Panel *panel_transform = GluiRight->add_panel((char*)"Rotate Graph" );
	GLUI_Rotation *rot = GluiRight->add_rotation_to_panel( panel_transform, (char*)" ", (float *) RotMatrix );
	rot->set_spin( 1.0 );

	// add a panel for user control buttons: start, pause, quit
	GLUI_Panel* panel_controls= GluiRight->add_panel((char*)"", false );
	GluiRight->add_button_to_panel( panel_controls, (char*)"Start", START, (GLUI_Update_CB) Buttons );
	GluiRight->add_button_to_panel( panel_controls, (char*)"Pause", PAUSE, (GLUI_Update_CB) Buttons );
	GluiRight->add_button_to_panel( panel_controls, (char*)"Quit", QUIT, (GLUI_Update_CB) Buttons );

	// add a check box that enables the user to dump data when quit
	GluiRight->add_statictext( (char *) "Dump Data When Quit");
	GluiRight->add_separator( );
	GluiRight->add_checkbox((char*)" ",  &dumpData);

	// set the main window
	GluiRight->set_main_gfx_window( MainWindow );	
	GLUI_Master.set_glutIdleFunc( NULL );
}

/* Initialize the glut and OpenGL libraries and setup callback functions
 * Code from Mike Bailey's graphics class files.
 */
void InitGraphics( )
{
	glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH );
	glutInitWindowPosition( 0, 0 );
	glutInitWindowSize(INIT_WINDOW_WIDTH, INIT_WINDOW_HEIGHT);

	// set and create the main window
	MainWindow = glutCreateWindow( WINDOWTITLE );
	glutSetWindowTitle( WINDOWTITLE );
	glClearColor( BACKCOLOR[0], BACKCOLOR[1], BACKCOLOR[2], BACKCOLOR[3] );

	// setup the callback routines:
	glutSetWindow( MainWindow );
	glutDisplayFunc( Display );
	glutReshapeFunc( Resize );
	glutMouseFunc( MouseButton );
	glutMotionFunc( MouseMotion );
	glutVisibilityFunc( Visibility );

#ifdef WIN32
	GLenum err = glewInit();
	if( err != GLEW_OK )
	{
		fprintf( stderr, "glewInit Error\n" );
	}
#endif
}

/* Setup the static graphics objects and store in GPU memory */
void InitLists( )
{
	// the frame that goes around the particle plots
	boxScaleList = glGenLists(1);
	glNewList(boxScaleList, GL_COMPILE);
		glColor3f(0.199,0.199,0.996);
		glLineWidth(AXES_WIDTH);
		gridBox(boxLength);
		glLineWidth(1.);
	glEndList();

	// the axes used in the field plots
	AxesList = glGenLists( 1 );
	glNewList( AxesList, GL_COMPILE );
		glPushMatrix();
		// align with the field graphs
		glTranslatef(fieldAxisPlotPosition[0], fieldAxisPlotPosition[1], fieldAxisPlotPosition[2]);
		glRotatef(fieldAxisPlotAngle, 0.0, 1.0, 0.0);
		glColor3fv( AXES_COLOR );
		glLineWidth( AXES_WIDTH );
			Axes( AXES_LENGTH );
		glLineWidth( 1. );
		glPopMatrix();
	glEndList( );
}

/* Used by glut callbacks for mouse button actions. This function is
 * taken straight from Mike Bailey's graphics class files
 */ 
void MouseButton( int button, int state, int x, int y )
{
	int b;			// LEFT, MIDDLE, or RIGHT
	
	switch( button )
	{
		case GLUT_LEFT_BUTTON:
			b = LEFT;		break;

		case GLUT_MIDDLE_BUTTON:
			b = MIDDLE;		break;

		case GLUT_RIGHT_BUTTON:
			b = RIGHT;		break;

		default:
			b = 0;
			fprintf( stderr, "Unknown mouse button: %d\n", button );
	}

	// button down sets the bit, up clears the bit:

	if( state == GLUT_DOWN )
	{
		Xmouse = x;
		Ymouse = y;
		ActiveButton |= b;		// set the proper bit
	}
	else
	{
		ActiveButton &= ~b;		// clear the proper bit
	}
}

/* Used by glut callbacks for mouse button actions. This function is
 * taken straight from Mike Bailey's graphics class files
 */ 
void MouseMotion( int x, int y )
{
	int dx = x - Xmouse;		// change in mouse coords
	int dy = y - Ymouse;

	if( ActiveButton & LEFT )
	{
			Xrot += ( ANGFACT*dy );
			Yrot += ( ANGFACT*dx );
	}

	Xmouse = x;			// new current position
	Ymouse = y;

	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}

/* Initialize parameters used for displaying graphs of simulation data
 * and global variables used by OpenGL.
 */
void InitParameters( )
{
	ActiveButton = 0;  // for mouse controls
	AxesOn = GLUIFALSE;
	Paused = GLUIFALSE;
	Scale  = 1.0;
	Scale2 = 0.0;		// because add 1. to it in Display( )
	
	WhichProjection = PERSP;
	dumpData = DO_NOT_DUMP;

	// rotation from Mike Bailey's graphics class files
	Xrot = Yrot = 0.;
	TransXYZ[0] = TransXYZ[1] = TransXYZ[2] = 0.;

	                  RotMatrix[0][1] = RotMatrix[0][2] = RotMatrix[0][3] = 0.;
	RotMatrix[1][0]                   = RotMatrix[1][2] = RotMatrix[1][3] = 0.;
	RotMatrix[2][0] = RotMatrix[2][1]                   = RotMatrix[2][3] = 0.;
	RotMatrix[3][0] = RotMatrix[3][1] = RotMatrix[3][3]                   = 0.;
	RotMatrix[0][0] = RotMatrix[1][1] = RotMatrix[2][2] = RotMatrix[3][3] = 1.;
}

/* Functions for handling window resize and visibibilty glut callbacks.
 * Code from Mike Bailey's graphics class files.
 */

void Resize( int width, int height )
{
	glutSetWindow( MainWindow );
	glutPostRedisplay( );
}

void Visibility ( int state )
{
	if( state == GLUT_VISIBLE )
	{
		glutSetWindow( MainWindow );
		glutPostRedisplay( );
	}
	
}

/* Generates a frame to enclose particle data graphs */
void gridBox(float length)
{
	float halfLength = length / 2.0;
	glBegin(GL_LINE_STRIP);
		glVertex3f(-halfLength, -halfLength, 0.);
		glVertex3f(-halfLength, halfLength, 0.);
		glVertex3f(halfLength, halfLength, 0.);
		glVertex3f(halfLength, -halfLength, 0.);
		glVertex3f(-halfLength, -halfLength, 0.);
	glEnd();
}

//
//	 Draw a set of 3D axes:
//  (length is the axis length in world coordinates)
//   Code from Mike Bailey's graphics class files.
void Axes( float length )
{
	// fraction of the length to use as height of the characters:
	const float LENFRAC = 0.10;

	// fraction of length to use as start location of the characters:
	const float BASEFRAC = 1.10;

	glBegin( GL_LINE_STRIP );
		glVertex3f( length, 0., 0. );
		glVertex3f( 0., 0., 0. );
		glVertex3f( 0., length, 0. );
	glEnd( );
	glBegin( GL_LINE_STRIP );
		glVertex3f( 0., 0., 0. );
		glVertex3f( 0., 0., length );
	glEnd( );

	float fact = (float)LENFRAC * length;
	float base = (float)BASEFRAC * length;
	
	glBegin( GL_LINE_STRIP );
		for( int i = 0; i < 4; i++ )
		{
			int j = xorder[i];
			if( j < 0 )
			{
				
				glEnd( );
				glBegin( GL_LINE_STRIP );
				j = -j;
			}
			j--;
			glVertex3f( base + fact*xx[j], fact*xy[j], 0.0 );
		}
	glEnd( );

	glBegin( GL_LINE_STRIP );
		for( int i = 0; i < 5; i++ )
		{
			int j = yorder[i];
			if( j < 0 )
			{
				
				glEnd( );
				glBegin( GL_LINE_STRIP );
				j = -j;
			}
			j--;
			glVertex3f( fact*yx[j], base + fact*yy[j], 0.0 );
		}
	glEnd( );

	glBegin( GL_LINE_STRIP );
		for( int i = 0; i < 6; i++ )
		{
			int j = zorder[i];
			if( j < 0 )
			{
				
				glEnd( );
				glBegin( GL_LINE_STRIP );
				j = -j;
			}
			j--;
			glVertex3f( 0.0, fact*zy[j], base + fact*zx[j] );
		}
	glEnd( );
	

}

/* Close down glut and glui in addition dump data if this option is selected
 * by the user. 
 */
void Quit( )
{	
	GluiRight->close();
	glutSetWindow( MainWindow );
	glFinish( );
	glutDestroyWindow( MainWindow );

	if (dumpData == DUMP)  // saves the current state of the simulation
		pic1.pic::dump_data();

	exit( 0 );
}









