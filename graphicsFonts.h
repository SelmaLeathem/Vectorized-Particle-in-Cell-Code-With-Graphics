/*
* Date: 2/26/2021
* Description: Functions and variable definitions used to draw fonts on the 
*             	main OpenGL graphics window. The code in this file is taken 
*               from Mike Bailey's graphics and parallel programming class 
*               files and has been replicated here with his permission. 
*/

#ifndef GRAPHICSFONTS_H
#define GRAPHICSFONTS_H

#define _CRT_SECURE_NO_WARNINGS 
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <cmath>

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

// used for the axes
// the stroke characters 'X' 'Y' 'Z' :

static float xx[] = {
		0., 1., 0., 1.
};

static float xy[] = {
		-.5, .5, .5, -.5
};

static int xorder[] = {
		1, 2, -3, 4
};


static float yx[] = {
		0., 0., -.5, .5
};

static float yy[] = {
		0., .6f, 1., 1.
};

static int yorder[] = {
		1, 2, 3, -2, 4
};


static float zx[] = {
		1., 0., 1., 0., .25, .75
};

static float zy[] = {
		.5, .5, -.5, -.5, 0., 0.
};

static int zorder[] = {
		1, 2, 3, 4, -5, 6
};

//Draws a raster string to the main graphics window
void	DoRasterString(float, float, float, char*);

//Draws a stroke sting to the main graphics window
void	DoStrokeString(float, float, float, float, char*);

#endif