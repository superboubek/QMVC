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
#ifndef CSPHERESVIEWER_H
#define CSPHERESVIEWER_H


//--------------------------------    DISCLAIMER    -------------------------------//
// We got rid of glew in this version, since it is not supported anymore by Qt5
//--------------------------------    DISCLAIMER    -------------------------------//

#include "openglincludeQtComp.h"

#include "QGLViewer/mouseGrabber.h"

#include <QGLViewer/qglviewer.h>

#include <QColorDialog>
#include <QFileDialog>


// io:
#include <string>
#include <iostream>
#include <fstream>


// containers:
#include <vector>
#include <set>
#include <unordered_set>


// assert:
#include <cassert>


class CameraParams{
    qglviewer::Vec pos;
    qglviewer::Vec view;
    qglviewer::Vec up;
    float fov;

public:
    void setParamsFromCamera( qglviewer::Camera * cam )
    {
        pos = cam->position();
        view = cam->viewDirection();
        up = cam->upVector();
        fov = cam->fieldOfView();
    }


    void setCameraFromParams( qglviewer::Camera * cam )
    {
        cam->setPosition(pos);
        cam->setViewDirection(view);
        cam->setUpVector(up);
        cam->setFieldOfView(fov);

        cam->computeModelViewMatrix();
        cam->computeProjectionMatrix();
    }
};



class BasicQGLViewer : public QGLViewer
{
    Q_OBJECT


private:
    void saveCameraInFile(const QString &filename){
        std::ofstream out (filename.toUtf8());
        if (!out)
            exit (EXIT_FAILURE);

        out << camera()->position() << " " <<
               camera()->viewDirection() << " " <<
               camera()->upVector() << " " <<
               camera()->fieldOfView();
        out << std::endl;
        out.close ();
    }

    void openCameraFromFile(const QString &filename){

        std::ifstream file;
        file.open(filename.toStdString().c_str());

        qglviewer::Vec pos;
        qglviewer::Vec view;
        qglviewer::Vec up;
        float fov;

        file >> (pos[0]) >> (pos[1]) >> (pos[2]) >>
                (view[0]) >> (view[1]) >> (view[2]) >>
                (up[0]) >> (up[1]) >> (up[2]) >>
                fov;

        camera()->setPosition(pos);
        camera()->setViewDirection(view);
        camera()->setUpVector(up);
        camera()->setFieldOfView(fov);

        camera()->computeModelViewMatrix();
        camera()->computeProjectionMatrix();

        updateGL();
    }


public :
    BasicQGLViewer(QWidget * parent = NULL) : QGLViewer(parent)
    {
    }

    template< class point_t >
    void adjustCamera( point_t const & p_bb , point_t const & p_BB )
    {
        point_t const & center = ( p_bb + p_BB )/2.f;
        setSceneCenter( qglviewer::Vec( center[0] , center[1] , center[2] ) );
        setSceneRadius( 0.5f * ( p_BB - p_bb ).norm() );
        showEntireScene();
    }

    void pickBackgroundColor()
    {
        QColor _bc = QColorDialog::getColor( this->backgroundColor(), this);
        if( _bc.isValid() )
        {
            this->setBackgroundColor( _bc );
            this->updateGL();
        }
    }

    virtual
    void drawGeometryIndices()
    {
        assert( 0  &&  "You need to overload BasicQGLViewer::drawGeometryIndices() , HAVE A LOOK AT THE FOLLOWING EXAMPLE" );

        // Here is an example of what you should do, in order to draw the index index_to_draw:

        int redShift = 16;
        int greenShift = 8;
        int blueShift = 0;
        GLuint redMask = 0xFF << redShift;
        GLuint greenMask = 0xFF << greenShift;
        GLuint blueMask = 0xFF;

        // the color that you need to use to associate a given vertex with an index:
        int index_to_draw = 15682;
        glColor4ub((index_to_draw & redMask) >> redShift, (index_to_draw & greenMask) >> greenShift, (index_to_draw & blueMask) >> blueShift, 255);
        // draw gl primitive
    }

    void setupGlFrustumMatrixToFitFocusZone( int xMin , int xMax , int yMin , int yMax )
    {
        // use glFrustrum to display full screen the zone corresponding to the focus zone:
        // todo

        double far = this->camera()->zFar();
        double near = this->camera()->zNear();
        double fov = this->camera()->fieldOfView();
        double aspect = (double)(this->width()) / (double)(this->height());

        double screenWidth_by_2 = (double)(this->width()) / 2.0; // it's also xScreenCenter
        double screenHeight_by_2 = (double)(this->height()) / 2.0; // it's also yScreenCenter


        double right = 0.0 , top = 0.0 ;

        // first compute the top, bottom, left and right values corresponding the current camera:
        top = tan(fov*0.5) * near;
        right = aspect * top;

        // now update these values based on the zoom on focus operation:
        double newtop = top * ((double)(yMax) - screenHeight_by_2) / screenHeight_by_2;
        double newbottom = top * ((double)(yMin) - screenHeight_by_2) / screenHeight_by_2;
        double newright = right * ((double)(xMax) - screenWidth_by_2) / screenWidth_by_2;
        double newleft = right * ((double)(xMin) - screenWidth_by_2) / screenWidth_by_2;

        glLoadIdentity();
        glFrustum(newleft , newright , newbottom , newtop , near , far);
    }

    void drawGeometryIndicesOnFocusZone( int xMin , int xMax , int yMin , int yMax )
    {
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        setupGlFrustumMatrixToFitFocusZone(xMin,xMax,yMin,yMax);
        drawGeometryIndices();
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
    }

    int pickedGeometryIndex(int x, int y)
    {
        glPushAttrib( GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT );

        glDisable(GL_BLEND);
        glDisable(GL_LIGHTING);
        glShadeModel(GL_FLAT);

        float back[4];
        glGetFloatv( GL_COLOR_CLEAR_VALUE , back );

        glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        drawGeometryIndices();

        GLubyte data[4];

        glReadPixels(x, y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, data);

        int index = ((data[0] << 16) + (data[1] << 8) + data[2]);

        // do we have geometry under the mouse or not: if it's index is FFFFFF, it means we pick up
        // the clear color. unless there are FFFFFF quads in the mesh, this should be safe
        bool isGeometry = !(index == 0xFFFFFF);

        // std::cout << std::hex << index << std::endl;

        if (!isGeometry || index < 0)
        {
            index = -1;
        }

        glClearColor(back[0], back[1], back[2], back[3]);

        glPopAttrib();

        return index;
    }

    void pickedGeometryIndexInDisk(int xCenter, int yCenter , float radius , std::set< int > & pickedIndices)
    {
        int xStart = (int)(floor( xCenter - radius )) ;
        int xEnd = (int)(ceil( xCenter + radius )) ;
        int yStart = (int)(floor( yCenter - radius )) ;
        int yEnd = (int)(ceil( yCenter + radius )) ;

        xStart = std::max( (int)(0) , xStart );
        xEnd = std::min( xEnd , this->width()-1 );
        yStart = std::max( (int)(0) , yStart );
        yEnd = std::min( yEnd , this->height()-1 );

        pickedIndices.clear();

        glPushAttrib( GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT );

        glDisable(GL_BLEND);
        glDisable(GL_LIGHTING);
        glShadeModel(GL_FLAT);

        float back[4];
        glGetFloatv( GL_COLOR_CLEAR_VALUE , back );

        glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        drawGeometryIndices();

        GLubyte fbdata[4 * (xEnd - xStart + 1) * (yEnd - yStart + 1)];
        glReadPixels(xStart, yStart, (xEnd - xStart + 1), (yEnd - yStart + 1), GL_RGBA, GL_UNSIGNED_BYTE, fbdata);

        for( int x = xStart ; x <= xEnd ; ++x )
        {
            for( int y = yStart ; y <= yEnd ; ++y )
            {
                if( (x-xCenter)*(x-xCenter) + (y-yCenter)*(y-yCenter) > radius*radius )
                    continue;

                GLubyte * data = &(fbdata[4 * ( (y-yStart) * (xEnd - xStart + 1) + x-xStart)]);

                int index = ((data[0] << 16) + (data[1] << 8) + data[2]);

                // do we have geometry under the mouse or not: if it's index is FFFFFF, it means we pick up
                // the clear color. unless there are FFFFFF quads in the mesh, this should be safe
                bool isGeometry = !(index == 0xFFFFFF);

                // std::cout << std::hex << index << std::endl;

                if (!isGeometry || index < 0)
                {
                    index = -1;
                }
                else
                {
                    pickedIndices.insert(index);
                }
            }
        }

        glClearColor(back[0], back[1], back[2], back[3]);

        glPopAttrib();
    }
    void pickedGeometryIndexInSquare(int xStart, int xEnd, int yStart, int yEnd , std::set< int > & pickedIndices)
    {
        pickedIndices.clear();

        xStart = std::max( (int)(0) , xStart );
        xEnd = std::min( xEnd , this->width()-1 );
        yStart = std::max( (int)(0) , yStart );
        yEnd = std::min( yEnd , this->height()-1 );

        glPushAttrib( GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT );

        glDisable(GL_BLEND);
        glDisable(GL_LIGHTING);
        glShadeModel(GL_FLAT);

        float back[4];
        glGetFloatv( GL_COLOR_CLEAR_VALUE , back );

        glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        drawGeometryIndices();

        GLubyte fbdata[4 * (xEnd - xStart + 1) * (yEnd - yStart + 1)];
        glReadPixels(xStart, yStart, (xEnd - xStart + 1), (yEnd - yStart + 1), GL_RGBA, GL_UNSIGNED_BYTE, fbdata);

        for( int x = xStart ; x <= xEnd ; ++x )
        {
            for( int y = yStart ; y <= yEnd ; ++y )
            {
                GLubyte * data = &(fbdata[4 * ( (y-yStart) * (xEnd - xStart + 1) + x-xStart)]);

                int index = ((data[0] << 16) + (data[1] << 8) + data[2]);

                // do we have geometry under the mouse or not: if it's index is FFFFFF, it means we pick up
                // the clear color. unless there are FFFFFF quads in the mesh, this should be safe
                bool isGeometry = !(index == 0xFFFFFF);

                // std::cout << std::hex << index << std::endl;

                if (!isGeometry || index < 0)
                {
                    index = -1;
                }
                else
                {
                    pickedIndices.insert(index);
                }
            }
        }

        glClearColor(back[0], back[1], back[2], back[3]);

        glPopAttrib();
    }







    void pickedGeometryIndexInDiskWithFocusOnDisk(int xCenter, int yCenter , float radius , std::unordered_set< int > & pickedIndices)
    {
        int xStart = (int)(floor( xCenter - radius )) ;
        int xEnd = (int)(ceil( xCenter + radius )) ;
        int yStart = (int)(floor( yCenter - radius )) ;
        int yEnd = (int)(ceil( yCenter + radius )) ;

        xStart = std::max( (int)(0) , xStart );
        xEnd = std::min( xEnd , this->width()-1 );
        yStart = std::max( (int)(0) , yStart );
        yEnd = std::min( yEnd , this->height()-1 );

        double xScreenCenter = (double)(this->width()-1) / 2.0;
        double yScreenCenter = (double)(this->height()-1) / 2.0;

        pickedIndices.clear();

        glPushAttrib( GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT );

        glDisable(GL_BLEND);
        glDisable(GL_LIGHTING);
        glShadeModel(GL_FLAT);

        float back[4];
        glGetFloatv( GL_COLOR_CLEAR_VALUE , back );

        glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        drawGeometryIndicesOnFocusZone( xStart , xEnd , yStart , yEnd );

        // the allocation takes 0 ms, so it's not worth trying to be clever about it
        GLubyte fbdata[4 * this->width() * this->height()];
        // in comparison, glReadPixels takes 10 ms when it's full screen on my laptop
        glReadPixels(0, 0, this->width(), this->height(), GL_RGBA, GL_UNSIGNED_BYTE, fbdata);

        for( int x = 0 ; x <= this->width()-1 ; ++x )
        {
            for( int y = 0 ; y <= this->height()-1 ; ++y )
            {
                // adapt the algebraic equation of the disk by taking into account the transformation which was setup
                if( (x-xScreenCenter)*(x-xScreenCenter)/(xScreenCenter*xScreenCenter) +
                        (y-yScreenCenter)*(y-yScreenCenter)/(yScreenCenter*yScreenCenter) > 1.0 )
                    continue;

                GLubyte * data = &(fbdata[4 * (y *this->width() + x)]);

                int index = ((data[0] << 16) + (data[1] << 8) + data[2]);

                // do we have geometry under the mouse or not: if it's index is FFFFFF, it means we pick up
                // the clear color. unless there are FFFFFF quads in the mesh, this should be safe
                bool isGeometry = !(index == 0xFFFFFF);

                // std::cout << std::hex << index << std::endl;

                if (!isGeometry || index < 0)
                {
                    index = -1;
                }
                else
                {
                    pickedIndices.insert(index);
                }
            }
        }

        glClearColor(back[0], back[1], back[2], back[3]);

        glPopAttrib();
    }



    // for debug:
    void saveBarycentricCellsScreenshot()
    {
        glPushAttrib( GL_LIGHTING_BIT | GL_COLOR_BUFFER_BIT );

        glDisable(GL_BLEND);
        glDisable(GL_LIGHTING);
        glShadeModel(GL_FLAT);

        float back[4];
        glGetFloatv( GL_COLOR_CLEAR_VALUE , back );

        glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        drawGeometryIndices();

        //        GLubyte data[4];

        //        for( int x = xStart ; x <= xEnd ; ++x )
        //        {
        //            for( int y = yStart ; y <= yEnd ; ++y )
        //            {
        //                glReadPixels(x, y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, data);

        //                int index = ((data[0] << 16) + (data[1] << 8) + data[2]);

        //                // do we have geometry under the mouse or not: if it's index is FFFFFF, it means we pick up
        //                // the clear color. unless there are FFFFFF quads in the mesh, this should be safe
        //                bool isGeometry = !(index == 0xFFFFFF);

        //                // std::cout << std::hex << index << std::endl;

        //                if (!isGeometry || index < 0)
        //                {
        //                    index = -1;
        //                }
        //                else
        //                {
        //                    pickedIndices.insert(index);
        //                }
        //            }
        //        }

        saveSnapshot("./barycentricCells.png");

        glClearColor(back[0], back[1], back[2], back[3]);

        glPopAttrib();
    }

    QImage get_QImage_with_alpha_from_screen()
    {
        GLubyte data[4 * width() * height()];
        glReadPixels(0, 0, width(), height(), GL_RGBA, GL_UNSIGNED_BYTE, data);
        QImage image(width(), height(), QImage::Format_ARGB32);

        QImage alpha = image.alphaChannel();

        for (int x = 0; x < width(); ++x)
            for (int y = 0; y < height(); ++y)
            {
                image.setPixel(QPoint(x, height() - y - 1), qRgb(data[4 * (x + y * width()) + 0] , data[4 * (x + y * width()) + 1] , data[4 * (x + y * width()) + 2]));
                alpha.setPixel(QPoint(x, height() - y - 1), data[4 * (x + y * width()) + 3]);
            }

        image.setAlphaChannel(alpha);
        return image;
    }

    void save_PNG_with_alpha_from_screen(const QString & filename)
    {
        GLubyte data[4 * width() * height()];
        glReadPixels(0, 0, width(), height(), GL_RGBA, GL_UNSIGNED_BYTE, data);
        QImage image(width(), height(), QImage::Format_ARGB32);

        QImage alpha = image.alphaChannel();

        for (int x = 0; x < width(); ++x)
            for (int y = 0; y < height(); ++y)
            {
                int r = (int)(data[4 * (x + y * width()) + 0]);
                int g = (int)(data[4 * (x + y * width()) + 1]);
                int b = (int)(data[4 * (x + y * width()) + 2]);
                image.setPixel(QPoint(x, height() - y - 1), qRgb( r , g , b ));
                if( r == 255  &&  g == 255  &&  b == 255 )
                    alpha.setPixel(QPoint(x, height() - y - 1), 0);
                else
                    alpha.setPixel(QPoint(x, height() - y - 1), (int)(data[4 * (x + y * width()) + 3]));
//                if( (x + y * width()) < 50 )
//                    std::cout << "alpha channel: " << (int)(data[4 * (x + y * width()) + 3]) << std::endl;
            }

        image.setAlphaChannel(alpha);
        image.save(filename,"PNG");
    }

    void init()
    {
        setMouseTracking(true);// Needed for MouseGrabber.
        printOpenGLError ();

        // standard background color:
        this->setBackgroundColor( Qt::white );

        // standard scene parameters:
        setSceneCenter( qglviewer::Vec( 0 , 0 , 0 ) );
        setSceneRadius( 10.f );
        showEntireScene();
    }


    void draw()
    {

    }


public slots:

    void openCamera(){
        QString fileName = QFileDialog::getOpenFileName(NULL,"","*.cam");
        if ( !fileName.isNull() ) {                 // got a file name
            openCameraFromFile(fileName);
        }
    }
    void saveCamera(){
        QString fileName = QFileDialog::getSaveFileName(NULL,"","*.cam");
        if ( !fileName.isNull() ) {                 // got a file name
            saveCameraInFile(fileName);
        }
    }

    void saveSnapShotPlusPlus(){
        QString fileName = QFileDialog::getSaveFileName(NULL,"*.png","");
        if ( !fileName.isNull() ) {                 // got a file name
            setSnapshotFormat("PNG");
            saveSnapshot( fileName );
            saveCameraInFile( fileName+QString(".cam") );
        }
    }

};












#endif // CSPHERESVIEWER_H
