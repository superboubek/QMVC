/*
    QMVC - Reference Implementation of paper: 
    
    "Mean Value Coordinates for Quad Cages in 3D", 
    jean-Marc Thiery, Pooran Memari and Tamy Boubekeur
    SIGGRAPH Asia 2018
    
    This program allows to compute QMVC for a set of 3D points contained 
    in a cage made of quad and triangles, as well as other flavors of 
    space coordinates for cages (MVC, SMVC, GC, MEC). It comes also with 
    a 3D viewer which helps deforming a mesh with a cage. 
    
    Copyright (C) 2018  jean-Marc Thiery, Pooran Memari and Tamy Boubekeur

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef OPENGLINCLUDEQTCOMP_H
#define OPENGLINCLUDEQTCOMP_H



#include <GL/glew.h>
#include <GL/glu.h>

#include <iostream>
#include <cstdlib>
#include <cstdio>

#if 1
        #define GetOpenGLError() t__GetOpenGLError( (char*)__FILE__, (int)__LINE__ )
#else
        #define GetOpenGLError()
#endif

inline int t__GetOpenGLError ( char* szFile, int iLine )
{
        int    retCode = 0;
        GLenum glErr = glGetError();
        while ( glErr != GL_NO_ERROR ) {
            std::cout << "GLError in file << " << szFile << " @ line " << iLine << ":" << gluErrorString( glErr ) << std::endl;
                retCode = 1;
                glErr = glGetError();
        }
        return retCode;
}




#define printOpenGLError() printOglError(__FILE__, __LINE__)

/// Returns 1 if an OpenGL error occurred, 0 otherwise.
static int printOglError (const char * file, int line) {
    GLenum glErr;
    int    retCode = 0;
    glErr = glGetError ();
    while (glErr != GL_NO_ERROR) {
        printf ("glError in file %s @ line %d: %s\n", file, line, gluErrorString(glErr));
        retCode = 1;
        glErr = glGetError ();
    }
    return retCode;
}





#endif // OPENGLINCLUDEQTCOMP_H
