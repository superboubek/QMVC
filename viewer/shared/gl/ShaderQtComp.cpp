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
#include "ShaderQtComp.h"

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <fcntl.h>

#include <sys/types.h>
#include <unistd.h>

#ifdef _WIN32 /*[*/
#include <io.h>
#endif /*]*/

using namespace std;


ShaderQtComp::ShaderQtComp () : QOpenGLFunctions_4_3_Core() , shaderProgram (0), vertexShader (0), fragmentShader (0),
vertexShaderSize (0), fragmentShaderSize (0), geometryShaderSize (0)
{
    initializeOpenGLFunctions();
}


ShaderQtComp::~ShaderQtComp () {
    if (hasVertexShader ())
        glDeleteShader (vertexShader);
    if (hasFragmentShader ())
        glDeleteShader (fragmentShaderSize);
    if (hasGeometryShader ())
        glDeleteShader (geometryShaderSize);
    glDeleteProgram (shaderProgram);
}


void ShaderQtComp::loadFromFile (const string & vertexShaderFilename, const string & fragmentShaderFilename, const string & geometryShaderFilename,
                           bool hasVSFile , bool hasFSFile , bool hasGSFile) {
    //printOpenGLError ();

    const GLchar * vertexShaderSource;
    const GLchar * fragmentShaderSource;
    const GLchar * geometryShaderSource;

    CheckExtension( "GL_EXT_geometry_shader4" );
    printOpenGLError ();

    if (hasVSFile)
    {
        vertexShaderSource = readShaderSource (vertexShaderFilename, vertexShaderSize);
        printOpenGLError ();
    }
    if (hasFSFile)
    {
        fragmentShaderSource = readShaderSource (fragmentShaderFilename, fragmentShaderSize);
        printOpenGLError ();
    }
    if (hasGSFile)
    {
        geometryShaderSource = readShaderSource (geometryShaderFilename, geometryShaderSize);
        printOpenGLError ();
    }
    shaderProgram = glCreateProgram ();
    printOpenGLError ();

    if (hasVertexShader () == true) {
//        std::cout << "compileAttach GL_VERTEX_SHADER" << std::endl;
        compileAttach (vertexShader, GL_VERTEX_SHADER, &vertexShaderSource);
        printOpenGLError ();
        delete [] vertexShaderSource;
    }
    if (hasFragmentShader () == true) {
//        std::cout << "compileAttach GL_FRAGMENT_SHADER" << std::endl;
        compileAttach (fragmentShader, GL_FRAGMENT_SHADER, &fragmentShaderSource);
        printOpenGLError ();
        delete [] fragmentShaderSource;
    }
    if (hasGeometryShader () == true) {
//        std::cout << "compileAttach GL_GEOMETRY_SHADER_EXT" << std::endl;
        compileAttach (geometryShader, GL_GEOMETRY_SHADER_EXT, &geometryShaderSource);
        printOpenGLError ();
        delete [] geometryShaderSource;
    }
//    std::cout << "glLinkProgram (shaderProgram)" << std::endl;
    glLinkProgram (shaderProgram);
    printOpenGLError ();
    GLint linked;
//    std::cout << "glGetProgramiv (shaderProgram, GL_LINK_STATUS, &linked)" << std::endl;
    glGetProgramiv (shaderProgram, GL_LINK_STATUS, &linked);
    printOpenGLError ();
//    std::cout << "printProgramInfoLog (shaderProgram)" << std::endl;
    printProgramInfoLog (shaderProgram);
    printOpenGLError ();
}


void ShaderQtComp::bind ( bool andSay ) {
    // glLinkProgram (shaderProgram);
    glUseProgram (shaderProgram);
    if( andSay )
        std::cout << "bind shaderProgram " << shaderProgram << std::endl;
}


void ShaderQtComp::unbind () {
    //glLinkProgram (0);
    glUseProgram (0);
}


// ------------------
// Protected methods.
// ------------------

unsigned int ShaderQtComp::getShaderSize (const string & filename) {
    int fd;
    unsigned int count = 0;
//#ifdef _WIN32
//    fd = _open (filename.c_str (), _O_RDONLY);
//    if (fd != -1) {
//        count = static_cast<unsigned int>(_lseek (fd, 0, SEEK_END) + 1);
//        _close(fd);
//    } else
//        throw ShaderException (string ("getShaderSize: bad Shader File Name") + filename);
//#else
    fd = open (filename.c_str (), O_RDONLY);
    if (fd != -1) {
        count = static_cast<unsigned int>(lseek (fd, 0, SEEK_END) + 1);
        close(fd);
    }
//#endif
    return count;
}


GLchar * ShaderQtComp::readShaderSource(const string & shaderFilename, unsigned int & shaderSize) {
    shaderSize = getShaderSize (shaderFilename);
    FILE * fh = fopen (shaderFilename.c_str (), "r");

    if(fh == NULL)
    {
        std::cout << shaderFilename << " could not be opened" << std::endl;
        // and then crash for now...
    }

    GLchar * shaderSource = new GLchar[shaderSize];
    fseek (fh, 0, SEEK_SET);
    int count = fread (shaderSource, 1, shaderSize, fh);
    shaderSource[count] = '\0';
    fclose (fh);

    return shaderSource;
}


GLint ShaderQtComp::getUniLoc (GLuint program, const GLchar *name)  {
    GLint loc = glGetUniformLocation (program, name);

    return loc;
}

GLint ShaderQtComp::getAttributLoc (GLuint program, const GLchar *name)
{
    GLint loc = glGetAttribLocation(program, name);

    return loc;
}



void ShaderQtComp::printShaderInfoLog (GLuint shader) {
    int infologLength = 0;
    int charsWritten  = 0;
    GLchar *infoLog;

    //  printOpenGLError ();  // Check for OpenGL errors
    glGetShaderiv (shader, GL_INFO_LOG_LENGTH, &infologLength);
    //printOpenGLError ();  // Check for OpenGL errors

    if (infologLength > 0) {
        infoLog = new GLchar[infologLength];
        if (infoLog == NULL) {
            printf("ERROR: Could not allocate InfoLog buffer\n");
            exit(1);
        }
        glGetShaderInfoLog (shader, infologLength, &charsWritten, infoLog);
        //    cerr << "InfoLog:" << endl << infoLog << endl << endl;
        delete [] infoLog;
    }
    // printOpenGLError();  // Check for OpenGL errors
}

void ShaderQtComp::printProgramInfoLog (GLuint program) {
    int infologLength = 0;
    int charsWritten  = 0;
    GLchar *infoLog;

    //   printOpenGLError ();  // Check for OpenGL errors
    glGetProgramiv (program, GL_INFO_LOG_LENGTH, &infologLength);
    //   printOpenGLError ();  // Check for OpenGL errors

    if (infologLength > 0) {
        infoLog = new GLchar[infologLength];
        if (infoLog == NULL) {
            printf("ERROR: Could not allocate InfoLog buffer\n");
            exit(1);
        }
        glGetProgramInfoLog (program, infologLength, &charsWritten, infoLog);
        //cerr << "InfoLog:" << endl << infoLog << endl << endl;
        delete [] infoLog;
    }
    //  printOpenGLError();  // Check for OpenGL errors
}

void ShaderQtComp::compileAttach (GLuint & shader, GLenum type, const GLchar ** source) {
    GLint shaderCompiled;
    shader = glCreateShader (type);
    glShaderSource (shader, 1, source, NULL);
    // printOpenGLError ();  // Check for OpenGL errors
    glCompileShader (shader);

//     int  iInfoLogLength = 0;
//    int  iCharsWritten  = 0;
//    char* szInfoLog = 0;
//    glGetObjectParameterivARB( shader, GL_OBJECT_INFO_LOG_LENGTH_ARB, &iInfoLogLength );
//    if ( iInfoLogLength > 1 ) {
//        szInfoLog = new char[ iInfoLogLength ];
//        std::cout << "You've got a problem in your shader : " << std::endl << (*source) << std::endl << std::endl << std::endl;
//        glGetInfoLogARB( shader, iInfoLogLength, &iCharsWritten, szInfoLog );
//        std::cout << szInfoLog << std::endl;
// //        char *m_strError = szInfoLog;

//        delete [] szInfoLog;
//        return;
//    }
    //printOpenGLError ();  // Check for OpenGL errors
    glGetShaderiv (shader, GL_COMPILE_STATUS, &shaderCompiled);
    //printOpenGLError ();  // Check for OpenGL errors
    printShaderInfoLog (shader);

    glAttachShader (shaderProgram, shader);
}
