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
#ifndef SMVC_H
#define SMVC_H





#include <vector>
#include <cmath>
#include <cassert>








namespace SMVC3D {

template< class point_t >
inline
double getAngleBetweenUnitVectors( point_t const & u1 , point_t const & u2 ) {
    return 2.0 * asin( (u1 - u2).norm() / 2.0 );
}

template< class point_t >
inline
double getTangentOfHalfAngleBetweenUnitVectors( point_t const & u1 , point_t const & u2 ) {
    double factor = (u1 - u2).norm() / 2.0;
    return factor / sqrt( std::max<double>(1 - factor * factor , 0.0) );
}


template< class int_t , class float_t , class point_t >
bool computeCoordinates(
        point_t const & eta ,
        std::vector< std::vector< int_t > > const & cage_faces ,
        std::vector< point_t > const & cage_vertices ,
        std::vector< float_t > & weights ,
        std::vector< float_t > & w_weights)// unnormalized weights
{
    typedef typename point_t::type_t    T;

    T epsilon = 0.000000001;

    unsigned int n_vertices = cage_vertices.size();
    unsigned int n_faces = cage_faces.size();

    w_weights.clear();
    weights.clear();
    weights.resize( n_vertices , 0.0 );
    T sumWeights = 0.0;

    std::vector< T > d( n_vertices , 0.0 );
    std::vector< point_t > u( n_vertices );

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        d[ v ] = ( eta - cage_vertices[ v ] ).norm();
        if( d[ v ] < epsilon )
        {
            weights[v] = 1.0;
            return true;
        }
        u[ v ] = ( cage_vertices[v] - eta ) / d[v];
    }

    w_weights.resize(  n_vertices , 0.0 );


    for( unsigned int f = 0 ; f < n_faces ; ++f )
    {
        // the Norm is CCW :
        point_t faceMeanVector;
        for( unsigned int i = 0 ; i < cage_faces[f].size() ; ++i ) {
            unsigned int v0 = cage_faces[f][i];
            unsigned int v1 = cage_faces[f][(i+1)%cage_faces[f].size()];

            point_t u0 = u[v0];
            point_t u1 = u[v1];

           // double angle = atan2( std::max<double>(-1, std::min<double>(1,point_t::cross(u0,u1).norm())) , std::max<double>(-1, std::min<double>(1,point_t::dot(u0,u1))) );
           // double angle = atan2( point_t::cross(u0,u1).norm() , point_t::dot(u0,u1) );
            double angle = getAngleBetweenUnitVectors( u0 , u1 );
            point_t n = point_t::cross(u0,u1).direction(true);

            faceMeanVector += (angle/2.0) * n;
        }

        double denominator = 0.0;
        std::vector<double> lambdas(cage_faces[f].size() , 0.0);

        for( unsigned int i = 0 ; i < cage_faces[f].size() ; ++i ) {
            unsigned int vi = cage_faces[f][i];
            unsigned int viplus1 = cage_faces[f][(i+1)%cage_faces[f].size()];
            unsigned int viminus1 = cage_faces[f][(i+cage_faces[f].size()-1)%cage_faces[f].size()];

            point_t ui = u[vi];
            point_t uiplus1 = u[viplus1];
            point_t uiminus1 = u[viminus1];

            // thetai = angle( vf , ui )
            double faceMeanVectorSqrNorm = faceMeanVector.sqrnorm();
            double vfnormDividedBysinthetai = faceMeanVectorSqrNorm;
            if(faceMeanVectorSqrNorm > epsilon) {
                vfnormDividedBysinthetai /= point_t::cross(faceMeanVector,ui).norm();
            }


            double tanAlphaiBy2 = getTangentOfHalfAngleBetweenUnitVectors( point_t::cross(faceMeanVector,ui).direction() , point_t::cross(faceMeanVector,uiplus1).direction() );

            double tanAlphaiMinus1By2 = getTangentOfHalfAngleBetweenUnitVectors( point_t::cross(faceMeanVector,uiminus1).direction() , point_t::cross(faceMeanVector,ui).direction() );

            double tangents = tanAlphaiBy2 + tanAlphaiMinus1By2;

            lambdas[i] = vfnormDividedBysinthetai * tangents / d[ vi ];

            denominator += tangents * point_t::dot(faceMeanVector,ui) / point_t::cross(faceMeanVector,ui).norm();
        }

        if(fabs(denominator) < epsilon) {
            // then we are on the face, and we output the unnormalized weights of the face, as they dominate in the final sum:
            w_weights.clear();
            weights.clear();
            sumWeights = 0.0;
            for( unsigned int i = 0 ; i < cage_faces[f].size() ; ++i ) {
                unsigned int vi = cage_faces[f][i];
                double lambdai = lambdas[i];
                weights[vi] = lambdai;
                sumWeights += lambdai;
            }
            for( unsigned int i = 0 ; i < cage_faces[f].size() ; ++i ) {
                unsigned int vi = cage_faces[f][i];
                weights[vi] /= sumWeights;
            }
            return true;
        }

        for( unsigned int i = 0 ; i < cage_faces[f].size() ; ++i ) {
            unsigned int vi = cage_faces[f][i];
            double lambdai = lambdas[i] / denominator;
            w_weights[vi] += lambdai;
            sumWeights += lambdai;
        }

    }

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
        weights[v]  = w_weights[v] / sumWeights;

    return false;
}



template< class int_t , class float_t , class point_t >
bool computeCoordinates(
        point_t const & eta ,
        std::vector< std::vector< int_t > > const & cage_faces ,
        std::vector< point_t > const & cage_vertices ,
        std::vector< float_t > & weights)
{
    std::vector< float_t > w_weights;
    return computeCoordinates(eta,cage_faces,cage_vertices,weights,w_weights);
}


} // namespace SMVC3D










#endif // MVCOORDINATES3D_H
