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

#include "BasicColors.h"

namespace RGB
{
int nColor = 30;

float color[ 30 ][3] = {
    { 0.52459 , 0.370169 , 0.818006 } ,
    { 0.22797 , 0.736751 , 0.578866 } ,
    { 0.615946 , 0.840482 , 0.941527 } ,
    { 0.820386 , 0.886748 , 0.611765 } ,
    { 0.523095 , 0.277974 , 0.277974 } ,
    { 0.160586 , 0.399985 , 0.656397 } ,
    { 0.563622 , 0.205234 , 0.501854 } ,
    { 0.948363 , 0.724559 , 0.358358 } ,
    { 0.377066 , 0.299641 , 0.860731 } ,
    { 0.29955 , 0.770504 , 0.824872 } ,
    { 0.659052 , 0.409354 , 0.861982 } ,
    { 0.559197 , 0.574395 , 0.133989 } ,
    { 0.513832 , 0.308644 , 0.163775 } ,
    { 0.825589 , 0.496513 , 0.564584 } ,
    { 0.800473 , 0.423407 , 0.657435 } ,
    { 0.347524 , 0.813748 , 0.765499 } ,
    { 0.501244 , 0.87425 , 0.44152 } ,
    { 0.493843 , 0.905364 , 0.692546 } ,
    { 0.304311 , 0.506661 , 0.197818 } ,
    { 0.520516 , 0.124697 , 0.534661 } ,
    { 0.551827 , 0.655329 , 0.927062 } ,
    { 0.64477 , 0.907103 , 0.321981 } ,
    { 0.924224 , 0.835676 , 0.410391 } ,
    { 0.271122 , 0.30605 , 0.777478 } ,
    { 0.720943 , 0.366674 , 0.274281 } ,
    { 0.596414 , 0.334386 , 0.334386 } ,
    { 0.293492 , 0.801221 , 0.328527 } ,
    { 0.409415 , 0.607416 , 0.464027 } ,
    { 0.861845 , 0.425895 , 0.606302 } ,
    { 0.49691 , 0.196338 , 0.592523 }
};

float color4[ 30 ][4] = {
    { 0.52459 , 0.370169 , 0.818006 , 1.f } ,
    { 0.22797 , 0.736751 , 0.578866 , 1.f } ,
    { 0.615946 , 0.840482 , 0.941527 , 1.f } ,
    { 0.820386 , 0.886748 , 0.611765 , 1.f } ,
    { 0.523095 , 0.277974 , 0.277974 , 1.f } ,
    { 0.160586 , 0.399985 , 0.656397 , 1.f } ,
    { 0.563622 , 0.205234 , 0.501854 , 1.f } ,
    { 0.948363 , 0.724559 , 0.358358 , 1.f } ,
    { 0.377066 , 0.299641 , 0.860731 , 1.f } ,
    { 0.29955 , 0.770504 , 0.824872 , 1.f } ,
    { 0.659052 , 0.409354 , 0.861982 , 1.f } ,
    { 0.559197 , 0.574395 , 0.133989 , 1.f } ,
    { 0.513832 , 0.308644 , 0.163775 , 1.f } ,
    { 0.825589 , 0.496513 , 0.564584 , 1.f } ,
    { 0.800473 , 0.423407 , 0.657435 , 1.f } ,
    { 0.347524 , 0.813748 , 0.765499 , 1.f } ,
    { 0.501244 , 0.87425 , 0.44152 , 1.f } ,
    { 0.493843 , 0.905364 , 0.692546 , 1.f } ,
    { 0.304311 , 0.506661 , 0.197818 , 1.f } ,
    { 0.520516 , 0.124697 , 0.534661 , 1.f } ,
    { 0.551827 , 0.655329 , 0.927062 , 1.f } ,
    { 0.64477 , 0.907103 , 0.321981 , 1.f } ,
    { 0.924224 , 0.835676 , 0.410391 , 1.f } ,
    { 0.271122 , 0.30605 , 0.777478 , 1.f } ,
    { 0.720943 , 0.366674 , 0.274281 , 1.f } ,
    { 0.596414 , 0.334386 , 0.334386 , 1.f } ,
    { 0.293492 , 0.801221 , 0.328527 , 1.f } ,
    { 0.409415 , 0.607416 , 0.464027 , 1.f } ,
    { 0.861845 , 0.425895 , 0.606302 , 1.f } ,
    { 0.49691 , 0.196338 , 0.592523 , 1.f }
};
}
