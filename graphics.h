/*
 * Date: 2/26/2021
 * Description: Graphics function and variable definitions for the purpose of
 *              implementing a glui user interface with simulation graph output
 *              using OpenGL for the graphics.
 *       
 *              Note: Indicated parts of the code are from files given in Mike
 *              Bailey's graphics and parallel programming classes and have been
 *              replicated here with his permission.
 * 
 *              References for GLUI:
 *              Mike Bailey's class code (CS475)
 *              http://www-h.eng.cam.ac.uk/help/tpl/graphics/using_glui.html
 *              user manual: https://github.com/libglui/glui
*/

#pragma once

#ifndef GRAPHICS_H
#define GRAPHICS_H

#define _CRT_SECURE_NO_WARNINGS 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <cmath>
#include <omp.h>
#include <vector>
#include "parameters.h"
#include "pic.h"
#include "particleVec.h"

#define _USE_MATH_DEFINES
#include <math.h>

#ifdef WIN32
#include <windows.h>
#pragma warning(disable:4996)
#endif

#ifdef WIN32
#include "glew.h"
#endif

#include <GL/gl.h>
#include <GL/glu.h>
#include "glut.h"
#include "glui.h"

#define GLM_FORCE_RADIANS
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"

#include "glslprogram.h"

#include "cl.h"
#include "cl_gl.h"

#include "parameters.h"
#include "particleVec.h"
#include "graphicsFonts.h"

/*
 * For Graphics Mode set the time between plot updates and time to take averages over
 * eg the amount of time over which to take averages for plots that are averaged over time
 * which includes all the field plots
 * WARNING: NOT RECOMMENDED TO GO ABOVE 5-10 FOR SMOOTH GRAPHICS OUTPUT 
*/
#define PLOT_UPDATE_INTERVAL (int)5 

/*
* The number of particles to display in plots of particle quantities.
*/
#define NUM_PARTICLES_TO_PLOT (int)5000

/*
* Through glut callbacks gets the latest data for the graphs by calling the 
* simulation's pic class run function.
*/
void	Animate();  

/* Implements the axes */
void	Axes(float);

/* Coordinates flags from pressing glui buttons to the graphics code */
void	Buttons(int);

/* Where everything is drawn to the main graphics window. */
void	Display();

/* Sets up the Glui user interface and its controls 
*  Some of the the code is from graphics class files.
*/
void	InitGlui();

/* Initialize the glut and OpenGL libraries and setup callback functions
 * Code from Mike Bailey's graphics class files.
 */
void	InitGraphics();

/* Setup the static graphics objects and store in GPU memory */
void	InitLists();

/* Used by glut callbacks for mouse button actions. These functions are
 * taken straight from Mike Bailey's graphics class files
 */ 
void	MouseButton(int, int, int, int);
void	MouseMotion(int, int);

/* Close down glut and glui in addition dump data if this option is selected
 * by the user. 
 */
void	Quit();

/* Initialize parameters used for displaying graphs of simulation data
 * and global variables used by OpenGL.
 */
void InitParameters(void); 

/* Functions for handling window resize and visibibilty glut callbacks.
 * Code from Mike Bailey's graphics class files.
 */
void	Resize(int, int);
void	Visibility(int);

/* Initialize the particle in cell simulation */
void InitPIC(void);

/* Generates a frame to enclose particle data graphs */
void gridBox(float length);

/* Handles glui radio buttons */
void GraphRadioGroup(int);

/* Draws field graph quantities to the main window */
void drawField(double **field, double fieldMax,int gridStart, char *strLabel, char *strTitle);

/* Draws particle phase graphs to the main window */
void drawParticlePhase(double **particle, glm::vec3 particleColor, char * strTitle);

/* Draws particle X-Y graphs to the main window */
void draw_particleXY(vector<particleVec> electrons, vector<particleVec> ions);

/* Resets the rotation of drawn quantities inside the main window back to zero */
void resetRotation(void);

#endif
