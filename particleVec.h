/*
 * File: particleVec.h
 * -------------------
 * This file defines a class that holds each particle's position, velocities,
 * and wieght. 
 *
 * 4/15
 */

#ifndef particleVec_h
#define particleVec_h

#define _CRT_SECURE_NO_WARNINGS   //Visual Studio Uses to Ignore VS related warnings
 //#ifdef WIN32
#include <windows.h>
//#else
//#include <unistd.h>
//#endif
#include <omp.h>

#include "CL/cl.h"
#include "CL/cl_platform.h"

using namespace std;

struct particleVec
{
      cl_double x;        /*position along x axis*/
      cl_double y;        /*position along y axis*/
      cl_double vx;        /*velocity along x axis*/
      cl_double vy;       /*y axis velocity*/
      cl_double vz;       /*z axis velocity*/
      cl_double weight;   /*wieght of particle*/

};
#endif /*particleVec*/

 
