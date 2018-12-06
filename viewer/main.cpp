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

#include "CageManip.h"
#include <QApplication>
#include <QMainWindow>
#include <QDockWidget>

double angle( point3d const & a1 , point3d const & a2 ) {
    return atan2( point3d::cross(a1,a2).norm() , point3d::dot(a1,a2) );
}

double area( point3d const & a1 , point3d const & a2 ) {
    return point3d::cross(a1,a2).norm();
}

double det(double a0 , double a1 , double b0 , double b1) {
    return a0 * b1 - b0 * a1;
}

void testBilinearQuantities() {
    point3d q0(0,0,0), q1(1,0,0) , q2(1.5,1,0) , q3(0,1,0.5);
    point3d p0(0,0,0) , p1(1,0,0), p2(1,1,0) , p3(0,1,0);

    unsigned int nSamples = 50;
    for( unsigned int s = 0 ; s < nSamples ; ++s ) {
        double u = (double)(rand()) / (double)(RAND_MAX);
        double v = (double)(rand()) / (double)(RAND_MAX);

        point3d quv = QuadUtilities::bilinearInterpolation(q0,q1,q2,q3,u,v);
        point3d puv = QuadUtilities::bilinearInterpolation(p0,p1,p2,p3,u,v);

        double a0 = angle(q1-quv,q0-quv) , a1 = angle(q2-quv,q1-quv) , a2 = angle(q3-quv,q2-quv) , a3 = angle(q0-quv,q3-quv) ;
        double b0 = angle(p1-puv,p0-puv) , b1 = angle(p2-puv,p1-puv) , b2 = angle(p3-puv,p2-puv) , b3 = angle(p0-puv,p3-puv) ;

        std::cout << det(a0,a1,b0,b1) << "    " << det(a0,a2,b0,b2) << "     " << det(a3,a1,b3,b1) << std::endl;

        a0 = area(q1-quv,q0-quv); a1 = area(q2-quv,q1-quv) ;a2 = area(q3-quv,q2-quv) ; a3 = area(q0-quv,q3-quv) ;
        b0 = area(p1-puv,p0-puv) ; b1 = area(p2-puv,p1-puv) ; b2 = area(p3-puv,p2-puv); b3 = area(p0-puv,p3-puv) ;

        std::cout << det(a0,a1,b0,b1) << "    " << det(a0,a2,b0,b2) << "     " << det(a3,a1,b3,b1) << std::endl;
    }
}


int main( int argc , char** argv )
{
    QApplication app( argc , argv );

    // testBilinearQuantities();
    QMainWindow * mainWindow = new QMainWindow;
    CMViewer * viewer = new CMViewer;
    viewer->setMinimumWidth (1920);
    viewer->setMinimumHeight (1080);
    mainWindow->setCentralWidget (viewer);
    auto * controlDockWidget = new QDockWidget (/*tr ("Informations"),*/ mainWindow);
    controlDockWidget->setWidget (viewer->getControlWidget ());
    mainWindow->addDockWidget (Qt::TopDockWidgetArea, controlDockWidget);
    controlDockWidget->setAllowedAreas (Qt::TopDockWidgetArea);
    controlDockWidget->setAllowedAreas (Qt::BottomDockWidgetArea);
    controlDockWidget->setFeatures (QDockWidget::NoDockWidgetFeatures);

    mainWindow->show();

    return app.exec();
}
