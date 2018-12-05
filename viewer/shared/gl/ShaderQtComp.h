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

#include "openglincludeQtComp.h"

#include <QOpenGLFunctions_4_3_Core>
#include <QOpenGLFunctions>
#include <QGLContext>


#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cassert>


class DummyShader : protected QOpenGLFunctions{
public:
    DummyShader() : QOpenGLFunctions()
    {
        std::cout << "1" << std::endl;
        QGLContext * contextflkj = new QGLContext( QGLFormat::defaultFormat() );
        assert(contextflkj  && "dfgfd");
        std::cout << "2" << std::endl;
        contextflkj->makeCurrent();
        std::cout << "3" << std::endl;
        initializeOpenGLFunctions();
        std::cout << "4" << std::endl;
    }
};




  
class ShaderQtComp  : protected QOpenGLFunctions_4_3_Core{
public:
    ShaderQtComp ();
    virtual ~ShaderQtComp ();
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
        bool CheckExtension( const std::string& strExt )
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
    void bind ( bool andSay = false );
    void unbind ();
        inline GLint getUniformLocation(const GLchar*name)  {return getUniLoc (shaderProgram, name); }
        inline GLint getAttributLocation(const GLchar*name)  {return getAttributLoc (shaderProgram, name); }
	
protected:
	GLchar * readShaderSource (const std::string & shaderFilename, unsigned int & shaderSize);

        GLint getUniLoc (GLuint program, const GLchar *name)  ;
        inline GLint getUniLoc( const GLchar*name)  { return getUniLoc (shaderProgram, name); }

	void compileAttach (GLuint & shader, GLenum type, const GLchar ** source);
    void printShaderInfoLog (GLuint shader);
    void printProgramInfoLog (GLuint program);
	/// Returns the size in bytes of the shader fileName. If an error occurred, it returns -1.
    unsigned int getShaderSize (const std::string & shaderFilename);

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
