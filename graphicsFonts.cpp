/*
* Date: 2/26/2021
* Description: Functions and implementations used to draw fonts on the 
*             	main OpenGL graphics window. The code in this file is taken 
*               from Mike Bailey's graphics and parallel programming class 
*               files and has been replicated here with his permission. 
*/

#include "graphicsFonts.h"

//
// use glut to display a string of characters using a stroke font:
//
void DoStrokeString( float x, float y, float z, float ht, char *s )
{
	char c;			// one character to print

	glPushMatrix( );
		glTranslatef( (GLfloat)x, (GLfloat)y, (GLfloat)z );
		float sf = ht / ( 119.05f + 33.33f );
		glScalef( (GLfloat)sf, (GLfloat)sf, (GLfloat)sf );
		for( ; ( c = *s ) != '\0'; s++ )
		{
			glutStrokeCharacter( GLUT_STROKE_ROMAN, c );
		}
	glPopMatrix( );
}

//
// use glut to display a string of characters using a raster font:
//
void DoRasterString( float x, float y, float z, char *s )
{
	char c;			// one character to print

	glRasterPos3f( (GLfloat)x, (GLfloat)y, (GLfloat)z );
	for( ; ( c = *s ) != '\0'; s++ )
	{
		glutBitmapCharacter( GLUT_BITMAP_TIMES_ROMAN_24, c );
	}
}