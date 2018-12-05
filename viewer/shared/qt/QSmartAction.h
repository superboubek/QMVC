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
#ifndef QSMARTACTION_H
#define QSMARTACTION_H


#include <QAction>
#include <QIcon>



class DetailedAction : public QAction
{
public :
    DetailedAction( const QIcon &icon, const QString &text, const QString &statusTip, QObject* parent , QObject * plugin , const char * connectedToTriggered) :  QAction( icon, text, parent )
    {
        this->setStatusTip(statusTip);
        connect(this, SIGNAL(triggered()), plugin, connectedToTriggered);
    }
};

class DetailedCheckableAction : public QAction
{
public :
    DetailedCheckableAction( const QIcon &icon, const QString &text, const QString &statusTip, QObject* parent , QObject * plugin , const char * connectedToTriggered) :  QAction( icon, text, parent )
    {
        this->setStatusTip(statusTip);
        this->setCheckable(true);
        connect(this, SIGNAL(toggled(bool)), plugin, connectedToTriggered);
    }
};

//DetailedAction( QIcon("./icons/open.png") , "Open OFF Mesh" , "Open OFF Mesh" , parent , this , SLOT(open_OFF_mesh()) );






#endif // QSMARTACTION_H
