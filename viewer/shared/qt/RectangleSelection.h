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
#ifndef RECTANGLESELECTION_H
#define RECTANGLESELECTION_H


#include <gl/openglincludeQtComp.h>
#include <QGLViewer/qglviewer.h>
#include <QMouseEvent>



class RectangleSelection : public QObject
{
    Q_OBJECT

    // Current rectangular selection
    QRect rectangle_;

    // Different selection modes
    enum SelectionMode { INACTIVE , ACTIVE , ADD , REMOVE };
    SelectionMode selectionMode_;

public:
    RectangleSelection() : selectionMode_( INACTIVE ) {}

    bool isInactive() const { return ( this->selectionMode_ == INACTIVE ); }

    void mousePressEvent(QMouseEvent* e , qglviewer::Camera* const)
    {
        if( this->selectionMode_ == INACTIVE )
            return;

        // Start selection. Mode is ADD by default and REMOVE with Ctrl key.
        rectangle_ = QRect(e->pos(), e->pos());

        if ((e->button() == Qt::LeftButton) && (e->modifiers() & Qt::ControlModifier))
        {
            selectionMode_ = REMOVE;
        }
        else if ((e->button() == Qt::LeftButton))
        {
            selectionMode_ = ADD;
        }
        else if ((e->button() == Qt::RightButton))
        {
            selectionMode_ = INACTIVE;
            emit apply();
        }
    }

    void mouseMoveEvent(QMouseEvent* e , qglviewer::Camera* const )
    {
        if(selectionMode_ == ACTIVE)
        {
            rectangle_.setTopLeft(e->pos());
        }
        else if(selectionMode_ != INACTIVE)
        {
            rectangle_.setBottomRight(e->pos());
        }
    }

    void mouseReleaseEvent(QMouseEvent* , qglviewer::Camera* const camera )
    {
        if (selectionMode_ != INACTIVE)
        {
            // Actual selection on the rectangular area.
            // Possibly swap left/right and top/bottom to make rectangle_ valid.
            rectangle_ = rectangle_.normalized();

            QRectF toEmit(
                    (rectangle_.center().x() - ((float)(rectangle_.width())/2.f))/ camera->screenWidth() ,
                    (rectangle_.center().y() - ((float)(rectangle_.height())/2.f))/ camera->screenHeight() ,
                    (float)(rectangle_.width())/ camera->screenWidth() ,
                    (float)(rectangle_.height())/ camera->screenHeight()
                    );

            // Use rectangle_.width() , rectangle_.height() , rectangle_.center() to emit the signal you want
            if( this->selectionMode_ == ADD )
                emit add( toEmit );
            else if( this->selectionMode_ == REMOVE )
                emit remove( toEmit );

            this->selectionMode_ = ACTIVE;
        }
    }

    void activate()
    {
        this->selectionMode_ = ACTIVE;
    }

    void draw() const
    {
        if( this->selectionMode_ == INACTIVE || this->selectionMode_ == ACTIVE )
            return;

        glDisable( GL_CLIP_PLANE0 );

        float viewport[4];
        glGetFloatv( GL_VIEWPORT , viewport );

        float w = viewport[2] , h = viewport[3];
        float left = (float)(rectangle_.left()) / w;
        float right = (float)(rectangle_.right()) / w;
        float top = 1.f - (float)(rectangle_.top()) / h;
        float bottom = 1.f - (float)(rectangle_.bottom()) / h;

        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        glEnable(GL_BLEND);

        glPolygonMode( GL_FRONT_AND_BACK , GL_FILL );

        glMatrixMode( GL_PROJECTION );
        glPushMatrix();
        glLoadIdentity();
        glOrtho( 0.f , 1.f , 0.f , 1.f , -1.f , 1.f );
        glMatrixMode( GL_MODELVIEW );
        glPushMatrix();
        glLoadIdentity();

        if( this->selectionMode_ == ADD )
            glColor4f(0.2, 0.2, 0.8f , 0.3f);
        if( this->selectionMode_ == REMOVE )
            glColor4f(0.8, 0.2, 0.2f , 0.3f);

        glBegin(GL_QUADS);
        glVertex2f( left , top );
        glVertex2f( left , bottom );
        glVertex2f( right , bottom );
        glVertex2f( right , top );
        glEnd();


        glLineWidth(2.0);
        if( this->selectionMode_ == ADD )
            glColor4f(0.1, 0.1, 1.f , 0.5f);
        if( this->selectionMode_ == REMOVE )
            glColor4f(0.9, 0.1, 0.1f , 0.5f);
        glBegin(GL_LINE_LOOP);
        glVertex2f( left , top );
        glVertex2f( left , bottom );
        glVertex2f( right , bottom );
        glVertex2f( right , top );
        glEnd();


        glPopMatrix();
        glMatrixMode( GL_PROJECTION );
        glPopMatrix();
        glMatrixMode( GL_MODELVIEW );


        glDisable(GL_BLEND);
        glEnable(GL_LIGHTING);
        glEnable(GL_DEPTH_TEST);
    }



    public slots:


    signals:
    void add( QRectF rectangle );
    void remove( QRectF rectangle );
    void apply();
};







#endif // RECTANGLESELECTION_H
