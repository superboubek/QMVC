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
#ifndef GLUTILITYMETHODS_H
#define GLUTILITYMETHODS_H


#include <GL/glew.h>
#include <QGLViewer/qglviewer.h>

#include "BasicColors.h"


namespace RGB
{
    void calc_RGB( float val , float val_min , float val_max , float & r , float & g , float & b );
    void get_RGB_from_HSV( float & r , float & g , float & b , float H , float S , float V );
    void get_random_RGB_from_HSV( float & r , float & g , float & b );
    void get_random_RGB_from_HSV( float & r , float & g , float & b , float H);
    void fromLightGrayToRed( float & r , float & g , float & b , float val , float greyValue = 0.75f, float maxSat=0.8f );
    void fromLightGrayToBlue( float & r , float & g , float & b , float val , float greyValue = 0.75f, float maxSat=0.8f );
}


namespace TriangleCoverage
{
extern float weights15[15][3];
extern float weights45[45][3];
void coutBasisFunctions( unsigned int subdiv );

unsigned int nElements( unsigned int N_splits_on_each_edge );
void getBarycentricCoordinates( unsigned int N_splits_on_each_edge , unsigned int e , float & x , float & y , float & z );
}

namespace BasicGL
{
    void drawSphere(float x,float y,float z,float radius,int slices,int stacks);

    int optimalSlices( float radius , float referenceRadius );
    int optimalStacks( float radius , float referenceRadius );

            template< class point_t >
            void glBox( point_t const & b0 ,
                        point_t const & b1 ,
                        point_t const & b2 ,
                        point_t const & b3 ,
                        point_t const & b4 ,
                        point_t const & b5 ,
                        point_t const & b6 ,
                        point_t const & b7 )
    {
        // bottom :
        glVertex( b0 );  glVertex( b1 );  glVertex( b2 );  glVertex( b3 );

        // top:
        glVertex( b4 );  glVertex( b5 );  glVertex( b6 );  glVertex( b7 );

        // front:
        glVertex( b0 );  glVertex( b1 );  glVertex( b5 );  glVertex( b4 );

        // back:
        glVertex( b2 );  glVertex( b3 );  glVertex( b7 );  glVertex( b6 );

        // left:
        glVertex( b0 );  glVertex( b4 );  glVertex( b7 );  glVertex( b3 );

        // right:
        glVertex( b1 );  glVertex( b2 );  glVertex( b6 );  glVertex( b5 );
    }
}






namespace GLTools{
    void initLights ();
    void setSunsetLight ();
    void setSunriseLight ();
    void setSingleSpotLight ();
    void setDefaultMaterial ();
    void setInverseDefaultMaterial ();
    void SetDiffuse(float r, float g , float b,float alpha);
    void SetAmbient(float r, float g , float b,float alpha);



    void color4(unsigned int c_id);
}


namespace GLCheck{
    int checkErrors( std::string const & szFile , int iLine );
}


#endif // GLUTILITYMETHODS_H
