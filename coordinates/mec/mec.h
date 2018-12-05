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
#ifndef MEC_H
#define MEC_H

#include <vector>
#include <cmath>
#include <cassert>




namespace MEC{
template< class mat_t , class point_t , class float_t , class float_t2 >
bool computeCoordinates(
        point_t eta ,
        std::vector< point_t > const & cage_vertices_in ,
        std::vector< float_t > const & i_masses,
        std::vector< float_t2 > & weights)
{
    assert( cage_vertices_in.size() == i_masses.size()   &&   "You have to provide masses w.r.t. the cage vertices" );
    unsigned int n_vertices = cage_vertices_in.size();
    point_t lambda(0,0,0); // initialization
    double epsilonTermination = 0.000001; // why not
    double errorTrigerringLineSearch = 0.0001; // why not
    unsigned int nbStepsForLineSearch = 10; // why not
    unsigned int max_iterations = 20; // why not...
    unsigned int it = 0;
    for(  ; it < max_iterations ; ++it ) {
        point_t gZ(0,0,0);
        mat_t HZ(0,0,0,0,0,0,0,0,0);
        double Z = 0.0;
        for( unsigned int i = 0 ; i < n_vertices ; ++i ) {
            point_t vi_bar = cage_vertices_in[i] - eta;
            double Zi = i_masses[i] * exp( - point_t::dot(lambda , vi_bar) );
            Z += Zi;
            gZ += - Zi * vi_bar;
            HZ += Zi * mat_t::tensor( vi_bar , vi_bar );
        }
        point_t gF = gZ / Z;
        mat_t HF = HZ / Z - mat_t::tensor( gF , gF );


        // Newton search direction:
        bool isPseudoInv;
        point_t dLambda = - HF.solveLinearSystem( gF , isPseudoInv );

        double errorCurrent = 0.0;
        for( unsigned int i = 0 ; i < n_vertices ; ++i ) {
            point_t vi_bar = cage_vertices_in[i] - eta;
            double Zi = i_masses[i] * exp( - point_t::dot(lambda + dLambda , vi_bar) );
            errorCurrent += Zi;
        }
        errorCurrent = fabs(log(errorCurrent));

        if(errorCurrent > errorTrigerringLineSearch) {
            double minError = FLT_MAX;
            point_t nextLambda;
            for( unsigned int step = 1 ; step < nbStepsForLineSearch ; ++step ) {
                double s = (double)(step) / (double)(nbStepsForLineSearch);
                point_t lambdaCurrent = lambda + s * dLambda;
                errorCurrent = 0.0;
                for( unsigned int i = 0 ; i < n_vertices ; ++i ) {
                    point_t vi_bar = cage_vertices_in[i] - eta;
                    double Zi = i_masses[i] * exp( - point_t::dot(lambdaCurrent , vi_bar) );
                    errorCurrent += Zi;
                }
                errorCurrent = fabs(log(errorCurrent));
                if( minError > errorCurrent ) {
                    errorCurrent = minError;
                    nextLambda = lambdaCurrent;
                }
            }
            lambda = nextLambda;
        }
        else{
            lambda += dLambda;
        }

        if( gF.norm() < epsilonTermination )
            break;
    }

    if(weights.size() != n_vertices) weights.resize(n_vertices);
    double sumW = 0.0;
    for( unsigned int i = 0 ; i < n_vertices ; ++i ) {
        point_t vi_bar = cage_vertices_in[i] - eta;
        double Zi = i_masses[i] * exp( - point_t::dot(lambda , vi_bar) );
        weights[i] = Zi;
        sumW += Zi;
    }

    bool isNan = false;
    for( unsigned int i = 0 ; i < n_vertices ; ++i ) {
        weights[i] /= sumW;
                if( std::isnan(weights[i]) ) {
                    weights[i] = 0.0;
                    isNan = true;
                }
    }
    return isNan;
}
}



#endif // MEC_H
