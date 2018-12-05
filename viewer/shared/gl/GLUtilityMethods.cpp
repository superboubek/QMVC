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
#include "GLUtilityMethods.h"

#if 1
#define GetOpenGLError() __GetOpenGLError( ( char* )__FILE__, ( int )__LINE__ )
#else
#define GetOpenGLError()
#endif



inline int __GetOpenGLError ( char* szFile, int iLine )
{
    int iRetCode = 0;
    GLenum glErr = glGetError();
    while ( glErr != GL_NO_ERROR ) {
        std::cout << "GLError in file << " << szFile << " @ line " << iLine << ":" << gluErrorString( glErr ) << std::endl;
        iRetCode = 1;
        glErr = glGetError();
    }
    return iRetCode;
}

namespace GLCheck{
    int checkErrors( std::string const & szFile , int iLine )
    {
        int iRetCode;
        GLenum glErr = glGetError();
        while ( glErr != GL_NO_ERROR ) {
            std::cout << "GLError in file << " << szFile << " @ line " << iLine << " : " << gluErrorString( glErr ) << std::endl;
            iRetCode = 1;
            glErr = glGetError();
        }
        return iRetCode;
    }
}


namespace TriangleCoverage
{
float weights15[15][3] =
{
    { 1       , 0       , 0 },
    { 0       , 1       , 0 },
    { 0       , 0       , 1 },
    { 3.f/4.f , 1.f/4.f , 0 },
    { 1.f/2.f , 1.f/2.f , 0 },
    { 1.f/4.f , 3.f/4.f , 0 },
    { 0       , 3.f/4.f , 1.f/4.f },
    { 0       , 1.f/2.f , 1.f/2.f },
    { 0       , 1.f/4.f , 3.f/4.f },
    { 3.f/4.f , 0       , 1.f/4.f },
    { 1.f/2.f , 0       , 1.f/2.f },
    { 1.f/2.f , 1.f/4.f , 1.f/4.f },
    { 1.f/4.f , 1.f/2.f , 1.f/4.f },
    { 1.f/4.f , 1.f/4.f , 1.f/2.f }
};


float weights45[45][3] =
{
    { 1   ,   0   ,   0 },
    { 0   ,   1   ,   0 },
    { 0   ,   0   ,   1 },
    { 0.875   ,   0.125   ,   0 },
    { 0.75   ,   0.25   ,   0 },
    { 0.625   ,   0.375   ,   0 },
    { 0.5   ,   0.5   ,   0 },
    { 0.375   ,   0.625   ,   0 },
    { 0.25   ,   0.75   ,   0 },
    { 0.125   ,   0.875   ,   0 },
    { 0.875   ,   0   ,   0.125 },
    { 0.75   ,   0.125   ,   0.125 },
    { 0.625   ,   0.25   ,   0.125 },
    { 0.5   ,   0.375   ,   0.125 },
    { 0.375   ,   0.5   ,   0.125 },
    { 0.25   ,   0.625   ,   0.125 },
    { 0.125   ,   0.75   ,   0.125 },
    { 0   ,   0.875   ,   0.125 },
    { 0.75   ,   0   ,   0.25 },
    { 0.625   ,   0.125   ,   0.25 },
    { 0.5   ,   0.25   ,   0.25 },
    { 0.375   ,   0.375   ,   0.25 },
    { 0.25   ,   0.5   ,   0.25 },
    { 0.125   ,   0.625   ,   0.25 },
    { 0   ,   0.75   ,   0.25 },
    { 0.625   ,   0   ,   0.375 },
    { 0.5   ,   0.125   ,   0.375 },
    { 0.375   ,   0.25   ,   0.375 },
    { 0.25   ,   0.375   ,   0.375 },
    { 0.125   ,   0.5   ,   0.375 },
    { 0   ,   0.625   ,   0.375 },
    { 0.5   ,   0   ,   0.5 },
    { 0.375   ,   0.125   ,   0.5 },
    { 0.25   ,   0.25   ,   0.5 },
    { 0.125   ,   0.375   ,   0.5 },
    { 0   ,   0.5   ,   0.5 },
    { 0.375   ,   0   ,   0.625 },
    { 0.25   ,   0.125   ,   0.625 },
    { 0.125   ,   0.25   ,   0.625 },
    { 0   ,   0.375   ,   0.625 },
    { 0.25   ,   0   ,   0.75 },
    { 0.125   ,   0.125   ,   0.75 },
    { 0   ,   0.25   ,   0.75 },
    { 0.125   ,   0   ,   0.875 },
    { 0   ,   0.125   ,   0.875 }
};


void coutBasisFunctions( unsigned int N )
{
    unsigned int nPoints = 0;
    std::cout << "{" << std::endl;
    for( unsigned int k = 0 ; k <= N ; ++k )
    {
        for( unsigned int j = 0 ; j <= N-k ; ++j )
        {
            std::cout << "    { " << ( 1.f - (float)(j)/(float)(N) - (float)(k)/(float)(N) )
                      << "   ,   " << ( (float)(j)/(float)(N) )
                      << "   ,   " << ( (float)(k)/(float)(N) )
                      << " }," << std::endl;
            ++nPoints;
        }
    }
    std::cout << "};" << std::endl;
    std::cout << nPoints << std::endl;
}



unsigned int nElements( unsigned int N_splits_on_each_edge )
{
    return ((N_splits_on_each_edge+1) * (N_splits_on_each_edge+2)) / 2;
}


void getBarycentricCoordinates( unsigned int N_splits_on_each_edge , unsigned int e , float & x , float & y , float & z )
{
    unsigned int i = (unsigned int) (  floor(  sqrt( 2.f * (float)(e) )  )  );
    if( (i * (i+1)) / 2 > e )
        --i;

    unsigned int k = N_splits_on_each_edge - i;
    unsigned int j = e  -  ((i * (i+1)) / 2);

    x = ( 1.f - (float)(j)/(float)(N_splits_on_each_edge) - (float)(k)/(float)(N_splits_on_each_edge) );
    y = ( (float)(j)/(float)(N_splits_on_each_edge) );
    z = ( (float)(k)/(float)(N_splits_on_each_edge) );
}
}



namespace RGB
{

    void get_random_RGB( float & r , float & g , float & b )
    {
        r = (float)(rand()) / (float)(RAND_MAX);
        g = (float)(rand()) / (float)(RAND_MAX);
        b = (float)(rand()) / (float)(RAND_MAX);
    }

    void get_random_RGB_from_HSV( float & r , float & g , float & b )
    {
        QColor c = QColor::fromHsvF( (float)(rand()) / (float)(RAND_MAX) , 0.3f + 0.5f * (float)(rand()) / (float)(RAND_MAX) , 0.5f + 0.5f * (float)(rand()) / (float)(RAND_MAX) );

        r = c.redF();
        g = c.greenF();
        b = c.blueF();
    }

    void get_RGB_from_HSV( float & r , float & g , float & b , float H , float S , float V )
    {
        QColor c = QColor::fromHsvF( H , S, V );

        r = c.redF();
        g = c.greenF();
        b = c.blueF();
    }

    void get_random_RGB_from_HSV( float & r , float & g , float & b , float H )
    {
        get_RGB_from_HSV( r , g , b , H , 0.3f + 0.5f * (float)(rand()) / (float)(RAND_MAX) , 0.5f + 0.5f * (float)(rand()) / (float)(RAND_MAX) );
    }

    void fromLightGrayToRed( float & r , float & g , float & b , float val , float greyValue , float maxSat  )
    {
        // non saturated red: (hue=0 -> red, value=1 -> light)
        get_RGB_from_HSV( r , g , b , 0.f , maxSat , 1.0 );

        r = val * r + (1.f-val) * greyValue;
        g = val * g + (1.f-val) * greyValue;
        b = val * b + (1.f-val) * greyValue;
    }

    void fromLightGrayToBlue( float & r , float & g , float & b , float val , float greyValue , float maxSat  )
    {
        // non saturated red: (hue=0 -> red, value=1 -> light)
        get_RGB_from_HSV( r , g , b , 0.5f , maxSat , 1.0 );

        r = val * r + (1.f-val) * greyValue;
        g = val * g + (1.f-val) * greyValue;
        b = val * b + (1.f-val) * greyValue;
    }

    void calc_RGB( float val , float val_min , float val_max , float & r , float & g , float & b )
    {
        // define uniform color intervalls [v0,v1,v2,v3,v4]
        float v0, v1, v2, v3, v4 ;

        v0 = val_min ;
        v1 = val_min + 1.0/4.0 * (val_max - val_min);
        v2 = val_min + 2.0/4.0 * (val_max - val_min);
        v3 = val_min + 3.0/4.0 * (val_max - val_min);
        v4 = val_max ;



        if (val < v0)
        {
            r = 0.f;
            g = 0.f;
            b = 1.f;
            return;
        }
        else if (val > v4)
        {
            r = 1.f;
            g = 0.f;
            b = 0.f;
            return;
        }
        else if (val <= v2)
        {
            if (val <= v1) // [v0, v1]
            {
                r = 0.f;
                g = (val - v0) / (v1 - v0);
                b = 1.f;
                return;
            }
            else // ]v1, v2]
            {
                r = 0.f;
                g = 1.f;
                b = 1.f - (val - v1) / (v2 - v1);
                return;
            }
        }
        else
        {
            if (val <= v3) // ]v2, v3]
            {
                r = (val - v2) / (v3 - v2);
                g = 1.f;
                b = 0.f;
                return;
            }
            else // ]v3, v4]
            {
                r = 1.f;
                g = 1.f - (val - v3) / (v4 - v3);
                b = 0.f;
                return;
            }
        }
    }
}










namespace BasicGL
{

int optimalSlices( float radius , float referenceRadius )
{
    return std::min( (int)(15) , std::max( (int)( 3 ) , (int)(10 * radius / referenceRadius ) ) );
}
int optimalStacks( float radius , float referenceRadius )
{
    return std::min( (int)(15) , std::max( (int)( 3 ) , (int)(10 * radius / referenceRadius ) ) );
}
    void drawSphere(float x,float y,float z,float radius,int slices,int stacks)
    {
        if(stacks < 2){stacks = 2;}
        if(stacks > 30){stacks = 30;}
        if(slices < 3){slices = 3;}
        if(slices > 30){slices = 30;}
        //Pas essentiel ...

        int Nb = slices*stacks +2;
        std::vector< qglviewer::Vec > points(Nb);

        qglviewer::Vec centre(x,y,z);

        float sinP , cosP , sinT , cosT , Phi , Theta;
        points[0] = qglviewer::Vec(0,0,1);
        points[Nb-1] = qglviewer::Vec(0,0,-1);

        for(int i=1; i<=stacks; i++)
        {
            Phi = 90 - (float)(i*180)/(float)(stacks+1);
            sinP = sinf(Phi*3.14159265/180);
            cosP = cosf(Phi*3.14159265/180);

            for(int j=1; j<=slices; j++)
            {
                Theta = (float)(j*360)/(float)(slices);
                sinT = sinf(Theta*3.14159265/180);
                cosT = cosf(Theta*3.14159265/180);

                points[ j + (i-1)*slices ] = qglviewer::Vec(cosT*cosP,sinT*cosP,sinP);
            }
        }

        int k1,k2;
        glBegin(GL_TRIANGLES);
        for(int i=1; i<=slices; i++)
        {
            k1 = i;
            k2 = (i%slices+1);
            glNormal3fv(points[0]);
            glVertex3fv((centre + radius*points[0]));
            glNormal3fv(points[k1]);
            glVertex3fv((centre + radius*points[k1]));
            glNormal3fv(points[k2]);
            glVertex3fv((centre + radius*points[k2]));

            k1 = (stacks-1)*slices+i;
            k2 = (stacks-1)*slices+(i%slices+1);
            glNormal3fv(points[k1]);
            glVertex3fv((centre + radius*points[k1]));
            glNormal3fv(points[Nb-1]);
            glVertex3fv((centre + radius*points[Nb-1]));
            glNormal3fv(points[k2]);
            glVertex3fv((centre + radius*points[k2]));
        }
        glEnd();

        glBegin(GL_QUADS);
        for(int j=1; j<stacks; j++)
        {
            for(int i=1; i<=slices; i++)
            {
                k1 = i + (j-1)*slices;
                k2 = (i%slices+1) + (j-1)*slices;
                glNormal3fv(points[k2]);
                glVertex3fv((centre + radius*points[k2]));
                glNormal3fv(points[k1]);
                glVertex3fv((centre + radius*points[k1]));

                k1 = i + (j)*slices;
                k2 = (i%slices+1) + (j)*slices;
                glNormal3fv(points[k1]);
                glVertex3fv((centre + radius*points[k1]));
                glNormal3fv(points[k2]);
                glVertex3fv((centre + radius*points[k2]));
            }
        }
        glEnd();
    }
}









namespace GLTools{
    void initLights () {
        GLfloat light_position0[4] = {0, 50, 50, 0};
        GLfloat light_position1[4] = {52, 16, 50, 0};
        GLfloat light_position2[4] = {26, 48, 50, 0};
        GLfloat light_position3[4] = {-16, 52, 50, 0};
        GLfloat light_position4[4] = {42, 374, 161, 0};
        GLfloat light_position5[4] = {473, -351, -259, 0};
        GLfloat light_position6[4] = {-438, 167, -48, 0};

        GLfloat direction1[3] = {-52,-16,-50};
        GLfloat direction2[3] = {-26,-48,-50};
        GLfloat direction3[3] = {16,-52,-50};
        GLfloat direction4[3] = {-42, -374, -161};
        GLfloat direction5[3] = {-473, 351, 259};
        GLfloat direction6[3] = {438, -167, 48};


        GLfloat color1[4] = {1,0, 0, 1};
        GLfloat color2[4] = {0, 1, 0, 1};
        GLfloat color3[4] = {0, 0, 1, 1};
        GLfloat color4[4] = {1, 1, 1, 1};
        GLfloat color5[4] = {0.28, 0.39, 1.0, 1};
        GLfloat color6[4] = {1.0, 0.69, 0.23, 1};

        GLfloat specularColor4[4] = {0.8, 0.8, 0.8, 1};
        GLfloat specularColor5[4] = {0.8, 0.8, 0.8, 1};
        GLfloat specularColor6[4] = {0.8, 0.8, 0.8, 1};

        GLfloat ambient[4] = {0.3f, 0.3f, 0.3f, 0.5f};

        glLightfv (GL_LIGHT0, GL_POSITION, light_position0);

        glLightfv (GL_LIGHT1, GL_POSITION, light_position1);
        glLightfv (GL_LIGHT1, GL_SPOT_DIRECTION, direction1);
        glLightfv (GL_LIGHT1, GL_DIFFUSE, color1);
        glLightfv (GL_LIGHT1, GL_SPECULAR, color1);

        glLightfv (GL_LIGHT2, GL_POSITION, light_position2);
        glLightfv (GL_LIGHT2, GL_SPOT_DIRECTION, direction2);
        glLightfv (GL_LIGHT2, GL_DIFFUSE, color2);
        glLightfv (GL_LIGHT2, GL_SPECULAR, color2);

        glLightfv (GL_LIGHT3, GL_POSITION, light_position3);
        glLightfv (GL_LIGHT3, GL_SPOT_DIRECTION, direction3);
        glLightfv (GL_LIGHT3, GL_DIFFUSE, color3);
        glLightfv (GL_LIGHT3, GL_SPECULAR, color3);

        glLightfv (GL_LIGHT4, GL_POSITION, light_position4);
        glLightfv (GL_LIGHT4, GL_SPOT_DIRECTION, direction4);
        glLightfv (GL_LIGHT4, GL_DIFFUSE, color4);
        glLightfv (GL_LIGHT4, GL_SPECULAR, specularColor4);

        glLightfv (GL_LIGHT5, GL_POSITION, light_position5);
        glLightfv (GL_LIGHT5, GL_SPOT_DIRECTION, direction5);
        glLightfv (GL_LIGHT5, GL_DIFFUSE, color5);
        glLightfv (GL_LIGHT5, GL_SPECULAR, specularColor5);

        glLightfv (GL_LIGHT6, GL_POSITION, light_position6);
        glLightfv (GL_LIGHT6, GL_SPOT_DIRECTION, direction6);
        glLightfv (GL_LIGHT6, GL_DIFFUSE, color6);
        glLightfv (GL_LIGHT6, GL_SPECULAR, specularColor6);

        glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambient);

        glEnable (GL_LIGHTING);
    }

    void setSunsetLight () {
        glDisable (GL_LIGHT0);
        glEnable (GL_LIGHT1);
        glEnable (GL_LIGHT2);
        glEnable (GL_LIGHT3);
        glDisable (GL_LIGHT4);
        glDisable (GL_LIGHT5);
        glDisable (GL_LIGHT6);
    }

    void setSunriseLight () {
        glDisable (GL_LIGHT0);
        glDisable (GL_LIGHT1);
        glDisable (GL_LIGHT2);
        glDisable (GL_LIGHT3);
        glEnable (GL_LIGHT4);
        glEnable (GL_LIGHT5);
        glEnable (GL_LIGHT6);
    }

    void setSingleSpotLight () {
        glEnable (GL_LIGHT0);
        glDisable (GL_LIGHT1);
        glDisable (GL_LIGHT2);
        glDisable (GL_LIGHT3);
        glDisable (GL_LIGHT4);
        glDisable (GL_LIGHT5);
        glDisable (GL_LIGHT6);
    }

    void setDefaultMaterial () {
        GLfloat material_color[4] = {1,1,1,1.0f};
        GLfloat material_specular[4] = {0.5,0.5,0.5,1.0};
        GLfloat material_ambient[4] = {0.0,0.0,0.0,0.0};

        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_color);
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128);
    }

    void setInverseDefaultMaterial () {
        GLfloat material_color[4] = {0,0,0,1.0f};
        GLfloat material_specular[4] = {0.5,0.5,0.5,1.0};
        GLfloat material_ambient[4] = {0.5,0.5,0.5,0.0};

        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_color);
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 128);
    }

    void SetDiffuse(float r, float g , float b,float alpha)
        {
            float mat[4] = {r,g,b,alpha};
            glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, mat );
        }
    void SetAmbient(float r, float g , float b,float alpha)
        {
            float mat[4] = {r,g,b,alpha};
            glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, mat );
        }



    void color4(unsigned int c_id)
    {
        glColor4fv(RGB::color4[c_id  %  RGB::nColor]);
    }

}

