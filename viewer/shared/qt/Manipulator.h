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
#ifndef SimpleManipulator_H
#define SimpleManipulator_H


#include <vector>
#include <map>

#include <QGLViewer/qglviewer.h>
#include <QGLViewer/mouseGrabber.h>
#include "gl/GLUtilityMethods.h"
#include <QMouseEvent>


class SimpleManipulator : public QObject , public qglviewer::MouseGrabber
{
    Q_OBJECT


    int State;
    //  0:  desactive.
    //  1:  il vient d'etre cree car l'user a clique sur une zone selectionnee quand aucun mode de selection n'etait actif

    int mode_modification;
    //  1:  tx
    //  2:  ty
    //  3:  tz
    //  4:  rx
    //  5:  ry
    //  6:  rz
    //  7:  sx
    //  8:  sy
    //  9:  sz

    bool mouse_released;

    qglviewer::Vec Origine;              // (1)
    qglviewer::Vec RepX, RepY, RepZ;    // (2)

    float uTeta , vTeta ;
    qglviewer::Vec PrevPos ;

    float display_scale;
    float Xscale , Yscale , Zscale;     // (3)

    std::vector< std::pair< int , qglviewer::Vec > > coordinates;
    std::map< int,int > idpoints;


public:
    SimpleManipulator()
    {
        mouse_released = true;
        Origine = qglviewer::Vec(0,0,0);
        RepX = qglviewer::Vec(1,0,0);
        RepY = qglviewer::Vec(0,1,0);
        RepZ = qglviewer::Vec(0,0,1);
        display_scale = 1.f;
        mode_modification = 0.f;
        Xscale = Yscale = Zscale = 1.f;
    }
    ~SimpleManipulator(){}


    void setState( int e ){ this->State = e; }
    int getState( ){ return this->State; }

    void activate(){ this->setState(1); }

    void setOrigine( qglviewer::Vec const & p ){ Origine = p; }
    void setRepX( qglviewer::Vec const & p ){ RepX = p; }
    void setRepY( qglviewer::Vec const & p ){ RepY = p; }
    void setRepZ( qglviewer::Vec const & p ){ RepZ = p; }


    void setDisplayScale(float ds){ display_scale = ds; }
    int getModification(){ return this->mode_modification; }

    void resetScales()
    {
        Xscale = Yscale = Zscale = 1.f;
    }

    void addPoint(int i , qglviewer::Vec const & p)
    {
        if( idpoints[i] == 0 )
        {
            qglviewer::Vec vr = p - Origine;
            coordinates.push_back( std::make_pair(i , qglviewer::Vec( vr * RepX , vr * RepY , vr * RepZ ) ) );
            idpoints[i] = coordinates.size();
        }
        else
        {
            qglviewer::Vec vr = p - Origine;
            coordinates[idpoints[i]-1] = std::make_pair(i , qglviewer::Vec( vr * RepX , vr * RepY , vr * RepZ ) );
        }
    }


    unsigned int n_points()
    {
        return this->coordinates.size();
    }


    void getTransformedPoint( unsigned int i , int & idx , qglviewer::Vec & p )
    {
        idx = this->coordinates[ i ].first;
        qglviewer::Vec vr = this->coordinates[ i ].second;
        p = Origine + vr[0] * Xscale * RepX + vr[1] * Yscale * RepY + vr[2] * Zscale * RepZ;
    }



    void draw()
    {
        glPolygonMode( GL_FRONT_AND_BACK , GL_FILL );

        // Dans BLENDER, ils affichent un modle de manipulateur diffrent suivant : mode_modification , mouse_released.
        if(State == 1)
        {
            float CX[3] = {0.f , 1.f , 0.f};
            float CY[3] = {0.f , 0.f , 1.f};
            float CZ[3] = {1.f , 1.f , 0.f};
            float Selection[3] = {1.f , 0.f , 0.f};

            glDisable(GL_LIGHTING);
            glLineWidth( 2.f );
            qglviewer::Vec p;
            glBegin( GL_LINES );
            if(mode_modification == 1)
                glColor3fv( Selection );
            else
                glColor3fv( CX );
            p = Origine - 2 * display_scale * RepX;
            glVertex3f(p[0],p[1],p[2]);
            p = Origine + 2 * display_scale * RepX;
            glVertex3f(p[0],p[1],p[2]);

            if(mode_modification == 2)
                glColor3fv( Selection );
            else
                glColor3fv( CY );
            p = Origine - 2 * display_scale * RepY;
            glVertex3f(p[0],p[1],p[2]);
            p = Origine + 2 * display_scale * RepY;
            glVertex3f(p[0],p[1],p[2]);

            if(mode_modification == 3)
                glColor3fv( Selection );
            else
                glColor3fv( CZ );
            p = Origine - 2 * display_scale * RepZ;
            glVertex3f(p[0],p[1],p[2]);
            p = Origine + 2 * display_scale * RepZ;
            glVertex3f(p[0],p[1],p[2]);
            glEnd();

            float teta;
            glBegin( GL_LINE_LOOP );
            if(mode_modification == 4)
                glColor3fv( Selection );
            else
                glColor3fv( CX );
            for(int i=0; i<360; i+=5)
            {
                teta = (float)(i) * 3.1415927 / float(180);
                p = Origine + display_scale*cosf(teta)*RepY + display_scale*sinf(teta)*RepZ;
                glVertex3f(p[0],p[1],p[2]);
            }
            glEnd();

            glBegin( GL_LINE_LOOP );
            if(mode_modification == 5)
                glColor3fv( Selection );
            else
                glColor3fv( CY );
            for(int i=0; i<360; i+=5)
            {
                teta = (float)(i) * 3.1415927 / float(180);
                p = Origine + display_scale*cosf(teta)*RepX + display_scale*sinf(teta)*RepZ;
                glVertex3f(p[0],p[1],p[2]);
            }
            glEnd();

            glBegin( GL_LINE_LOOP );
            if(mode_modification == 6)
                glColor3fv( Selection );
            else
                glColor3fv( CZ );
            for(int i=0; i<360; i+=5)
            {
                teta = (float)(i) * 3.1415927 / float(180);
                p = Origine + display_scale*cosf(teta)*RepY + display_scale*sinf(teta)*RepX;
                glVertex3f(p[0],p[1],p[2]);
            }
            glEnd();


            if(mode_modification == 7 || mode_modification == -7)
                glColor3fv( Selection );
            else
                glColor3fv( CX );
            p = Origine + (1.5 * Xscale * display_scale) * RepX;
            BasicGL::drawSphere(p[0],p[1],p[2],display_scale/15,5,5);
            p = Origine - (1.5 * Xscale * display_scale) * RepX;
            BasicGL::drawSphere(p[0],p[1],p[2],display_scale/15,5,5);


            if(mode_modification == 8 || mode_modification == -8)
                glColor3fv( Selection );
            else
                glColor3fv( CY );
            p = Origine + (1.5 * Yscale * display_scale) * RepY;
            BasicGL::drawSphere(p[0],p[1],p[2],display_scale/15,5,5);
            p = Origine - (1.5 * Yscale * display_scale) * RepY;
            BasicGL::drawSphere(p[0],p[1],p[2],display_scale/15,5,5);


            if(mode_modification == 9 || mode_modification == -9)
                glColor3fv( Selection );
            else
                glColor3fv( CZ );
            p = Origine + (1.5 * Zscale * display_scale) * RepZ;
            BasicGL::drawSphere(p[0],p[1],p[2],display_scale/15,5,5);
            p = Origine - (1.5 * Zscale * display_scale) * RepZ;
            BasicGL::drawSphere(p[0],p[1],p[2],display_scale/15,5,5);
        }
    }






    void checkIfGrabsMouse(int x, int y,const qglviewer::Camera* const cam)
    {
        qglviewer::Vec src;
        qglviewer::Vec img ;

        float lambda , epsilon_rotation_detect = 0.2f , epsilon_tranlation_detect = 0.005f;
        qglviewer::Vec eye,dir;
        qglviewer::Vec Eye , X , Z;
        qglviewer::Vec Dir , d , e;

        //Ce test marche, et permet maintenant de manipuler une sphre dans un gl, et un modle dans l'autre, ou la sphre de
        //l'autre gl en mme temps ... parfait :-) ...

        // Liste des tats successifs de la sphere:
        // 0 : absente -> alors setGrabsMouse(false);
        // 1 : positionning -> true; et on la dirige sur la surface du maillage.
        // 2 : (positionned and) directing -> true; son centre est fix, on la dirige.
        // 3 : fixed -> manhattanLength() <= 50 px; on devra rappuyer sur le bouton crer une sphre ou appuyer sur clic droit

        // Clic gauche : passer  l'tat suivant
        // Clic droit : reprendre l'tat prcdent
        switch(State)
        {
        case 0:
            setGrabsMouse(false);
            // Dans cet tat, la sphre est dsactive
            break;

        case 1:
            // Dans cet tat, on vrifie si l'utilisateur a cliqu dans une zone d'influence du SimpleManipulator
            if( ! mouse_released )
            {
                setGrabsMouse(true);
                return;
            }

            cam->convertClickToLine(QPoint(x,y),eye,dir);
            Eye = qglviewer::Vec(eye[0],eye[1],eye[2]);
            Dir = qglviewer::Vec(dir[0],dir[1],dir[2]);



            //--------  Dilatations:   --------//

            // Check on sx :
            X = Origine + (1.5 * Xscale*display_scale) * RepX;
            lambda = ( X-Eye )*( X-Eye ) - ( (X-Eye)*(Dir) ) * ( (X-Eye)*(Dir) )/(Dir * Dir);
            if( lambda < display_scale*display_scale / 100 )
            {
                mode_modification = 7;
                setGrabsMouse(true);
                return;
            }
            X = Origine - (1.5 * Xscale*display_scale) *RepX;
            lambda = ( X-Eye )*( X-Eye ) - ( (X-Eye)*(Dir) ) * ( (X-Eye)*(Dir) )/(Dir * Dir);
            if( lambda < display_scale*display_scale / 100 )
            {
                mode_modification = -7;
                setGrabsMouse(true);
                return;
            }

            // Check on sy :
            X = Origine + (1.5 * Yscale*display_scale) *RepY;
            lambda = ( X-Eye )*( X-Eye ) - ( (X-Eye)*(Dir) ) * ( (X-Eye)*(Dir) )/(Dir * Dir);
            if( lambda < display_scale*display_scale / 100 )
            {
                mode_modification = 8;
                setGrabsMouse(true);
                return;
            }
            X = Origine - (1.5 * Yscale*display_scale) *RepY;
            lambda = ( X-Eye )*( X-Eye ) - ( (X-Eye)*(Dir) ) * ( (X-Eye)*(Dir) )/(Dir * Dir);
            if( lambda < display_scale*display_scale / 100 )
            {
                mode_modification = -8;
                setGrabsMouse(true);
                return;
            }

            // Check on sz :
            X = Origine + (1.5 * Zscale*display_scale) *RepZ;
            lambda = ( X-Eye )*( X-Eye ) - ( (X-Eye)*(Dir) ) * ( (X-Eye)*(Dir) )/(Dir * Dir);
            if( lambda < display_scale*display_scale / 100 )
            {
                mode_modification = 9;
                setGrabsMouse(true);
                return;
            }
            X = Origine - (1.5 * Zscale*display_scale) *RepZ;
            lambda = ( X-Eye )*( X-Eye ) - ( (X-Eye)*(Dir) ) * ( (X-Eye)*(Dir) )/(Dir * Dir);
            if( lambda < display_scale*display_scale / 100 )
            {
                mode_modification = -9;
                setGrabsMouse(true);
                return;
            }







            //------  Rotations:   -------//

            // Check on rx :
            lambda = ( ( Origine - Eye )*( RepX ) )/( ( Dir )*( RepX ) );
            X = Eye + lambda*Dir;
            if( fabs( ( (X-Origine)*(X-Origine) )/(display_scale*display_scale) - 1 ) < epsilon_rotation_detect )
            {
                mode_modification = 4;
                setGrabsMouse(true);
                return;
            }

            // Check on ry :
            lambda = ( ( Origine - Eye )*( RepY ) )/( ( Dir )*( RepY ) );
            X = Eye + lambda*Dir;
            if( fabs( ( (X-Origine)*(X-Origine) )/(display_scale*display_scale) - 1 ) < epsilon_rotation_detect )
            {
                mode_modification = 5;
                setGrabsMouse(true);
                return;
            }

            // Check on rz :
            lambda = ( ( Origine - Eye )*( RepZ ) )/( ( Dir )*( RepZ ) );
            X = Eye + lambda*Dir;
            if( fabs( ( (X-Origine)*(X-Origine) )/(display_scale*display_scale) - 1 ) < epsilon_rotation_detect )
            {
                mode_modification = 6;
                setGrabsMouse(true);
                return;
            }






            //------  Translations:   ------//

            // Check on tx :
            d = cross( RepX , Dir );
            e = cross( Dir , d );
            lambda = ( ( Eye - Origine )*( e ) ) / ( RepX * e );
            X = Origine + lambda*RepX;
            if( lambda < 2.2*display_scale && lambda > -2.2*display_scale )
            {
                Z = Eye + ( Dir * ( X-Eye ) )*Dir/sqrt( (Dir*Dir) );
                if( ( ( Z-X )*( Z-X ) ) < display_scale*display_scale*epsilon_tranlation_detect )
                {
                    mode_modification = 1;
                    setGrabsMouse(true);
                    return;
                }
            }

            // Check on ty :
            d = cross( RepY , Dir );
            e = cross( Dir , d );
            lambda = ( ( Eye - Origine )*( e ) ) / ( RepY * e );
            X = Origine + lambda*RepY;
            if( lambda < 2.2*display_scale && lambda > -2.2*display_scale )
            {
                Z = Eye + ( Dir * ( X-Eye ) )*Dir/sqrt( (Dir*Dir) );
                if( ( ( Z-X )*( Z-X ) ) < display_scale*display_scale*epsilon_tranlation_detect )
                {
                    mode_modification = 2;
                    setGrabsMouse(true);
                    return;
                }
            }

            // Check on tz :
            d = cross( RepZ , Dir );
            e = cross( Dir , d );
            lambda = ( ( Eye - Origine )*( e ) ) / ( RepZ * e );
            X = Origine + lambda*RepZ;
            if( lambda < 2.2*display_scale && lambda > -2.2*display_scale )
            {
                Z = Eye + ( Dir * ( X-Eye ) )*Dir/sqrt( (Dir*Dir) );
                if( ( ( Z-X )*( Z-X ) ) < display_scale*display_scale*epsilon_tranlation_detect )
                {
                    mode_modification = 3;
                    setGrabsMouse(true);
                    return;
                }
            }

            mode_modification = 0;
            setGrabsMouse(false);
            break;
        }
    }

    void clear()
    {
        this->coordinates.clear();
        this->setState( 0 );
        this->idpoints.clear();
    }

    void manipulatedCallback()
    {
        emit moved();
    }


    void fakeMouseDoubleClickEvent( QMouseEvent* const )
    {
        if( mode_modification == 7 || mode_modification == -7 )
        {
            Xscale = 1.f;
            manipulatedCallback();
            return;
        }
        if( mode_modification == 8 || mode_modification == -8 )
        {
            Yscale = 1.f;
            manipulatedCallback();
            return;
        }
        if( mode_modification == 9 || mode_modification == -9 )
        {
            Zscale = 1.f;
            manipulatedCallback();
            return;
        }
    }

    void mousePressEvent( QMouseEvent* const event  , qglviewer::Camera* const cam )
    {
        mouse_released = false;

        if( event->buttons() & Qt::RightButton )
        {
            mouse_released = true;
            this->clear();
            this->setState( 0 );
            emit manipulatorReleased();
        }

        if( mode_modification > 6 || mode_modification < -6 )
        {
            return;
        }

        float lambda;
        qglviewer::Vec eye,dir;
        qglviewer::Vec Eye;
        qglviewer::Vec Ur,Vr,Dir , d , e;
        cam->convertClickToLine(event->pos(),eye,dir);
        Eye = qglviewer::Vec(eye[0],eye[1],eye[2]);
        Dir = qglviewer::Vec(dir[0],dir[1],dir[2]);

        if( mode_modification == 1 )
        {
            d = cross( RepX , Dir );
            e = cross( Dir , d );
            lambda = ( ( Eye - Origine )*( e ) ) / ( RepX * e );
            PrevPos = Origine + lambda*RepX;
            return;
        }
        if( mode_modification == 2 )
        {
            d = cross( RepY , Dir );
            e = cross( Dir , d );
            lambda = ( ( Eye - Origine )*( e ) ) / ( RepY * e );
            PrevPos = Origine + lambda*RepY;
            return;
        }
        if( mode_modification == 3 )
        {
            d = cross( RepZ , Dir );
            e = cross( Dir , d );
            lambda = ( ( Eye - Origine )*( e ) ) / ( RepZ * e );
            PrevPos = Origine + lambda*RepZ;
            return;
        }

        if( mode_modification == 4 )
        {
            // Alors on est en train de tourner autour de (Origine,RepX)
            lambda = ( ( Origine - Eye )* RepX )/( Dir * RepX );
            Ur = Eye + lambda*Dir - Origine;
            Ur.normalize();
            Vr = cross( RepX , Ur );
            Vr.normalize();
            // On a maintenant un repre Ur,Vr du plan orthogonal  RepX, et caractrisant X = le point cliqu sur le plan.

            uTeta = (RepY * Ur);
            vTeta = (RepY * Vr);
            return;
        }

        if( mode_modification == 5 )
        {
            // Alors on est en train de tourner autour de (Origine,RepY)
            lambda = ( ( Origine - Eye )* RepY )/( Dir * RepY );
            Ur = Eye + lambda*Dir - Origine;
            Ur.normalize();
            Vr = cross( RepY , Ur );
            Vr.normalize();
            // On a maintenant un repre Ur,Vr du plan orthogonal  RepY, et caractrisant X = le point cliqu sur le plan.

            uTeta = (RepZ * Ur);
            vTeta = (RepZ * Vr);
            return;
        }

        if( mode_modification == 6 )
        {
            // Alors on est en train de tourner autour de (Origine,RepZ)
            lambda = ( ( Origine - Eye )* RepZ )/( Dir * RepZ );
            Ur = Eye + lambda*Dir - Origine;
            Ur.normalize();
            Vr = cross( RepZ , Ur );
            Vr.normalize();
            // On a maintenant un repre Ur,Vr du plan orthogonal  RepZ, et caractrisant X = le point cliqu sur le plan.

            uTeta = (RepX * Ur);
            vTeta = (RepX * Vr);
            return;
        }
    }

    void mouseReleaseEvent( QMouseEvent* const , qglviewer::Camera* const  )
    {
        mouse_released = true;
        emit mouseReleased();
    }

    void mouseMoveEvent(QMouseEvent* const event, qglviewer::Camera* const cam)
    {
        if( ! mouse_released )
        {
            qglviewer::Vec eye,dir;
            qglviewer::Vec Eye , NewPos;
            qglviewer::Vec Dir , d , e ;
            qglviewer::Vec Ur , Vr;
            float lambda;

            cam->convertClickToLine(event->pos(),eye,dir);
            Eye = qglviewer::Vec(eye[0],eye[1],eye[2]);
            Dir = qglviewer::Vec(dir[0],dir[1],dir[2]);

            switch(mode_modification)
            {
            case 1:
                // Alors on doit trouver le point sur la droite (Origine,RepX) qui est le plus proche du rayon
                d = cross( RepX , Dir );
                e = cross( Dir , d );
                lambda = ( ( Eye - Origine ) * e ) / ( RepX * e );
                NewPos = Origine + lambda*RepX;
                Origine += NewPos - PrevPos;
                PrevPos = NewPos;
                manipulatedCallback();
                break;
            case 2:
                // Alors on doit trouver le point sur la droite (Origine,RepY) qui est le plus proche du rayon
                d = cross( RepY , Dir );
                e = cross( Dir , d );
                lambda = ( ( Eye - Origine ) * e ) / ( RepY * e );
                NewPos = Origine + lambda*RepY;
                Origine += NewPos - PrevPos;
                PrevPos = NewPos;
                manipulatedCallback();
                break;
            case 3:
                // Alors on doit trouver le point sur la droite (Origine,RepZ) qui est le plus proche du rayon
                d = cross( RepZ , Dir );
                e = cross( Dir , d );
                lambda = ( ( Eye - Origine ) * e ) / ( RepZ * e );
                NewPos = Origine + lambda*RepZ;
                Origine += NewPos - PrevPos;
                PrevPos = NewPos;
                manipulatedCallback();
                break;
            case 4:
                lambda = ( ( Origine - Eye ) * RepX )/( Dir * RepX );
                Ur = Eye + lambda*Dir - Origine;
                Ur.normalize();
                Vr = cross( RepX , Ur );
                Vr.normalize();

                RepY = uTeta * Ur + vTeta * Vr;
                RepY.normalize();
                RepZ = cross(RepX,RepY);
                manipulatedCallback();
                break;
            case 5:
                lambda = ( ( Origine - Eye ) * RepY )/( Dir * RepY );
                Ur = Eye + lambda*Dir - Origine;
                Ur.normalize();
                Vr = cross( RepY , Ur );
                Vr.normalize();

                RepZ = uTeta * Ur + vTeta * Vr;
                RepZ.normalize();
                RepX = cross(RepY,RepZ);
                manipulatedCallback();
                break;
            case 6:
                lambda = ( ( Origine - Eye ) * RepZ )/( Dir * RepZ );
                Ur = Eye + lambda*Dir - Origine;
                Ur.normalize();
                Vr = cross( RepZ , Ur );
                Vr.normalize();

                RepX = uTeta * Ur + vTeta * Vr;
                RepX.normalize();
                RepY = cross(RepZ,RepX);
                manipulatedCallback();
                break;

            case 7:
                d = cross( RepX , Dir );
                e = cross( Dir , d );
                lambda = ( ( Eye - Origine ) * e ) / ( RepX * e );
                Xscale = lambda / (1.5*display_scale);
                manipulatedCallback();
                break;
            case -7:
                d = cross( RepX , Dir );
                e = cross( Dir , d );
                lambda = ( ( Eye - Origine ) * e ) / ( RepX * e );
                Xscale = -lambda / (1.5*display_scale);
                manipulatedCallback();
                break;

            case 8:
                d = cross( RepY , Dir );
                e = cross( Dir , d );
                lambda = ( ( Eye - Origine ) * e ) / ( RepY * e );
                Yscale = lambda / (1.5*display_scale);
                manipulatedCallback();
                break;
            case -8:
                d = cross( RepY , Dir );
                e = cross( Dir , d );
                lambda = ( ( Eye - Origine ) * e ) / ( RepY * e );
                Yscale = -lambda / (1.5*display_scale);
                manipulatedCallback();
                break;

            case 9:
                d = cross( RepZ , Dir );
                e = cross( Dir , d );
                lambda = ( ( Eye - Origine ) * e ) / ( RepZ * e );
                Zscale = lambda / (1.5*display_scale);
                manipulatedCallback();
                break;
            case -9:
                d = cross( RepZ , Dir );
                e = cross( Dir , d );
                lambda = ( ( Eye - Origine ) * e ) / ( RepZ * e );
                Zscale = -lambda / (1.5*display_scale);
                manipulatedCallback();
                break;
            }
        }
    }







public slots:



signals:
    void mouseReleased();
    void moved();
    void manipulatorReleased();

};









#endif // SimpleManipulator_H
