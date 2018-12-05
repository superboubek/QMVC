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
#ifndef QSMARTCOLORPICKER_H
#define QSMARTCOLORPICKER_H

//#include <QtGui>
#include <QToolButton>
#include <QColorDialog>


class QSmartColorPicker : public QToolButton
{
    Q_OBJECT

private:
    QColor bgdColor;
    QWidget * m_parent;
    void changeBackgroundColor( QColor c )
    {
        bgdColor = c;
        setStyleSheet(QString("background-color: rgb(%1, %2, %3); color: rgb(255,255,255)").arg(c.red()).arg(c.green()).arg(c.blue()));
    }

public:
    QSmartColorPicker (QWidget * p = NULL) : QToolButton (p) , m_parent(p){
        setAutoFillBackground(true);
        connect( this, SIGNAL(clicked()) , this , SLOT(openColorDialog()) );
        setText("  ");
        changeBackgroundColor( QColor( 255,255,255 ) );
    }
    QSmartColorPicker ( QColor c ,QWidget * p = NULL) : QToolButton (p) , m_parent(p){
        setAutoFillBackground(true);
        connect( this, SIGNAL(clicked()) , this , SLOT(openColorDialog()) );
        setText("  ");
        changeBackgroundColor(c);
    }

signals:
    void colorChanged( QColor c );

public slots:
    void openColorDialog(  )
    {
        QColor color = QColorDialog::getColor( bgdColor , m_parent);
        if( color.isValid() )
        {
            changeBackgroundColor(color);
            emit colorChanged(color);
        }
    }
};







#endif // QSMARTCOLORPICKER_H
