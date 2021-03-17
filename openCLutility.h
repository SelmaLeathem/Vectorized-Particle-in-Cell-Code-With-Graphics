/*
 * Utility files for openCL. Includes a Wait() function for waiting
 * on all tasks to finish. These functions are from Mike Bailey's
 * CS 475 class files and are replicated here with permission.
 */

#ifndef OPENCL_UTILITY_H
#define OPENCL_UTILITY_H

#define _CRT_SECURE_NO_WARNINGS 
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
//#ifdef WIN32
#include <windows.h>
//#else
//#include <unistd.h>
//#endif
#include <omp.h>

#include "CL/cl.h"
#include "CL/cl_platform.h"


#define	LOCAL_SIZE	32
using namespace std;


// wait until all queued tasks have taken place:
void Wait( cl_command_queue queue );

// utility function to examine bits 
int LookAtTheBits( float fp );

#endif