/*
 * Utility files for openCL. Includes a Wait function for waiting
 * on a process to finish. These functions are from Mike Bailey's
 * CS 475 class files and are replicated here with permission.
 */

#include "openCLutility.h"

// wait until all queued tasks have taken place:

void
Wait( cl_command_queue queue )
{
      cl_event wait;
      cl_int      status;

      status = clEnqueueMarker( queue, &wait );
      if( status != CL_SUCCESS )
	      fprintf( stderr, "Wait: clEnqueueMarker failed\n" );

      status = clWaitForEvents( 1, &wait );
      if( status != CL_SUCCESS )
	      fprintf( stderr, "Wait: clWaitForEvents failed\n" );
}

// utility function to examine bits 
int
LookAtTheBits( float fp )
{
	int *ip = (int *)&fp;
	return *ip;
}