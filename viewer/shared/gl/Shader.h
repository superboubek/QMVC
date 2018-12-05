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


#pragma once
#include <string>

#define GLEW_STATIC 1
#include "openglinclude.h"



#include <iostream>
#include <cstdlib>
#include <cstdio>

//#define printOpenGLError() printOglError(__FILE__, __LINE__)

///// Returns 1 if an OpenGL error occurred, 0 otherwise.
//static int printOglError (const char * file, int line) {
//    GLenum glErr;
//    int    retCode = 0;
//    glErr = glGetError ();
//    while (glErr != GL_NO_ERROR) {
//        printf ("glError in file %s @ line %d: %s\n", file, line, gluErrorString(glErr));
//        retCode = 1;
//        glErr = glGetError ();
//    }
//    return retCode;
//}



  
class Shader {
public:
	Shader ();
	virtual ~Shader ();
        inline bool hasVertexShader () const { return (vertexShaderSize > 0); }
        inline bool hasFragmentShader () const { return (fragmentShaderSize > 0); }
        inline bool hasGeometryShader () const { return (geometryShaderSize > 0); }
	inline GLuint getShaderProgram () { return shaderProgram; }
        inline GLuint getVertexShader () { return vertexShader; }
        inline GLuint getFragmentShader () { return fragmentShader; }
        inline GLuint getGeometryShader () { return geometryShader; }
        void loadFromFile (const std::string & vertexShaderFilename,
                           const std::string & fragmentShaderFilename,
                           const std::string & geometryShaderFilename,
                           bool hasVSFile , bool hasFSFile , bool hasGSFile);
        void loadFromFile (const std::string & vertexShaderFilename,
                           const std::string & fragmentShaderFilename,
                           const std::string & geometryShaderFilename) {
            loadFromFile(vertexShaderFilename , fragmentShaderFilename , geometryShaderFilename , true , true , true);
        }

        void loadFromFile (const std::string & vertexShaderFilename,
                           const std::string & fragmentShaderFilename) {
            loadFromFile (vertexShaderFilename,fragmentShaderFilename, "", true , true , false);
            }
	inline void loadFromFile (const std::string & vertexShaderFilename) { 
            loadFromFile (vertexShaderFilename, "", "" , true , false , false);
	}
        static bool CheckExtension( const std::string& strExt )
         {
          GLint iExtensionCount = 0;
          glGetIntegerv( GL_NUM_EXTENSIONS, &iExtensionCount );
          for( GLint i=0; i < iExtensionCount; ++i ) {

           // signature: const GLubyte* (GLAPIENTRY * PFNGLGETSTRINGIPROC) (GLenum, GLuint);
           const GLubyte* pBytes = glGetStringi( GL_EXTENSIONS, i );
           std::string strGlExt( ( const char* ) pBytes );
           if( strGlExt == strExt ) {
            return true;
           }
          }
          return false;
         }
    void bind ( bool andSay = false ) const;
    void unbind () const;
        inline GLint getUniformLocation(const GLchar*name) const  {return getUniLoc (shaderProgram, name); }
        inline GLint getAttributLocation(const GLchar*name)  {return getAttributLoc (shaderProgram, name); }
	
protected:
	GLchar * readShaderSource (const std::string & shaderFilename, unsigned int & shaderSize);

        GLint getUniLoc (GLuint program, const GLchar *name) const  ;
        inline GLint getUniLoc( const GLchar*name)  { return getUniLoc (shaderProgram, name); }

	void compileAttach (GLuint & shader, GLenum type, const GLchar ** source);
	static void printShaderInfoLog (GLuint shader);
	static void printProgramInfoLog (GLuint program);
	/// Returns the size in bytes of the shader fileName. If an error occurred, it returns -1.
	static unsigned int getShaderSize (const std::string & shaderFilename);

        GLint getAttributLoc (GLuint program, const GLchar *name);
        inline GLint getAttributLoc( const GLchar*name) { return getAttributLoc(shaderProgram, name); }

	
private:
        GLuint shaderProgram, vertexShader, fragmentShader, geometryShader;
        unsigned int vertexShaderSize, fragmentShaderSize, geometryShaderSize;
};
  
  

// Some Emacs-Hints -- please don't remove:
//
//  Local Variables:
//  mode:C++
//  tab-width:4
//  End:
