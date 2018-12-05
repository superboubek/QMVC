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
#ifndef TEXTUREHANDLER_H
#define TEXTUREHANDLER_H

#include <QImage>

class ScalarTextureHandler{
    GLuint ui_text_id;

public:
    template< class glViewer >
    void initTexture( glViewer * viewer , const QString & filename = "./icons/texture.png")
    {
        glEnable(GL_TEXTURE_2D);

        glGenTextures(1, &(ui_text_id) );

        QImage textimg;
        QImage buf;
        buf.load(filename);
        textimg = viewer->convertToGLFormat( buf );
    //    textimg.load(filename);

        glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_MAG_FILTER , GL_LINEAR );
        glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_MIN_FILTER , GL_LINEAR );

        ui_text_id = viewer->bindTexture( textimg, GL_TEXTURE_2D, GL_RGB );

        glDisable(GL_TEXTURE_2D);
    }

    void enableTexturesAndBind()
    {
        glEnable(GL_TEXTURE_2D);
        glBindTexture( GL_TEXTURE_2D , ui_text_id );
        glColor3f(1,1,1);
    }
    void unBind()
    {
        glBindTexture( GL_TEXTURE_2D , 0 );
    }
    void bind()
    {
        glBindTexture( GL_TEXTURE_2D , ui_text_id );
        glColor3f(1,1,1);
    }
    void unBindAndDisableTextures()
    {
        glBindTexture( GL_TEXTURE_2D , 0 );
        glDisable(GL_TEXTURE_2D);
    }

    void setVal( float val )
    {
        glTexCoord2d(val, 0.f);
    }

    void setVal( float val , float valMin , float valMax )
    {
        val = std::min< float >( std::max< float >( val , valMin ) , valMax );
        glTexCoord2d((val-valMin) / (valMax - valMin), 0.f);
    }
};



class ScalarQuadTextureHandler{
    GLuint ui_text_id;

public:
    template< class glViewer >
    void initTexture( glViewer * viewer , QString textureString )
    {
        glEnable(GL_TEXTURE_2D);
        glEnable(GL_BLEND);

        glGenTextures(1, &(ui_text_id) );

        QImage textimg , buf;
        buf.load(textureString);

//        textimg = glViewer::convertToGLFormat( buf );
        textimg = buf;

        glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_MAG_FILTER , GL_LINEAR );
        glTexParameteri( GL_TEXTURE_2D , GL_TEXTURE_MIN_FILTER , GL_LINEAR );

        ui_text_id = viewer->bindTexture( textimg, GL_TEXTURE_2D, GL_RGBA );

        glDisable(GL_TEXTURE_2D);
    }

    void enableTexturesAndBind()
    {
        glEnable(GL_TEXTURE_2D);
        glColor4f(1,1,1,1);
        glBindTexture( GL_TEXTURE_2D , ui_text_id );
        glColor4f(1,1,1,1);
    }
    void unBind()
    {
        glBindTexture( GL_TEXTURE_2D , 0 );
    }
    void bind()
    {
        glBindTexture( GL_TEXTURE_2D , ui_text_id );
        glColor4f(1,1,1,1);
    }
    void unBindAndDisableTextures()
    {
        glBindTexture( GL_TEXTURE_2D , 0 );
        glDisable(GL_TEXTURE_2D);
    }

    void setCoords( float u , float v )
    {
        glTexCoord2d(u,v);
    }
};


#endif // TEXTUREHANDLER_H
