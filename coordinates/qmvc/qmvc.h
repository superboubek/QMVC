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
#ifndef QMVC_H
#define QMVC_H

#include <vector>
#include <cmath>
#include <cassert>

#include "point4.h"
#include "point3.h"

#include "quadutilities.h"



namespace QMVC3D{
template< class tri_t , class quad_t , class float_t >
bool computeCoordinates(
        point3d eta ,
        std::vector< tri_t > const & cage_triangles ,
        std::vector< quad_t > const & cage_quads ,
        std::vector< point3d > const & cage_vertices_in ,
        std::vector< float_t > & weights,
        std::vector< float_t > & w_weights)
{
    double epsilon = 0.000001;

    std::vector< point3d > const & cage_vertices = cage_vertices_in;

    unsigned int n_vertices = cage_vertices.size();
    unsigned int n_triangles = cage_triangles.size();
    unsigned int n_quads = cage_quads.size();

    w_weights.clear();
    weights.clear();
    weights.resize( n_vertices , 0.0 );
    double sumWeights = 0.0;

    std::vector< double > d( n_vertices , 0.0 );
    std::vector< point3d > u( n_vertices );

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

    // CAGE TRIANGLES:
    {
        unsigned int vid[3];
        double l[3]; double theta[3] ; double w[3];

        for( unsigned int t = 0 ; t < n_triangles ; ++t )
        {
            // the Norm is CCW :
            for( unsigned int i = 0 ; i <= 2 ; ++i )
                vid[i] =  cage_triangles[t][i];

            for( unsigned int i = 0 ; i <= 2 ; ++i )
                l[ i ] = ( u[ vid[ ( i + 1 ) % 3 ] ] - u[ vid[ ( i + 2 ) % 3 ] ] ).norm();

            for( unsigned int i = 0 ; i <= 2 ; ++i )
                theta[i] = 2.0 * asin( l[i] / 2.0 );

            double determinant = point3d::dot( cage_vertices[vid[0]] - eta , point3d::cross( cage_vertices[vid[1]] - cage_vertices[vid[0]] , cage_vertices[vid[2]] - cage_vertices[vid[0]] ) );
            double sqrdist = determinant*determinant / (4 * point3d::cross( cage_vertices[vid[1]] - cage_vertices[vid[0]] , cage_vertices[vid[2]] - cage_vertices[vid[0]] ).sqrnorm() );
            double dist = sqrt( sqrdist );

            if( dist < epsilon ) {
                // then the point eta lies on the support plane of the triangle
                double h = ( theta[0] + theta[1] + theta[2] ) / 2.0;
                if( M_PI - h < epsilon ) {
                    // eta lies inside the triangle t , use 2d barycentric coordinates :
                    for( unsigned int i = 0 ; i <= 2 ; ++i ) {
                        w[ i ] = sin( theta[ i ] ) * l[ (i+2) % 3 ] * l[ (i+1) % 3 ];
                    }
                    sumWeights = w[0] + w[1] + w[2];

                    w_weights.clear();
                    weights.clear();
                    weights.resize( n_vertices , 0.0 );
                    weights[ vid[0] ] = w[0] / sumWeights;
                    weights[ vid[1] ] = w[1] / sumWeights;
                    weights[ vid[2] ] = w[2] / sumWeights;
                    return true;
                }
            }

            point3d pt[3] , N[3];

            for( unsigned int i = 0 ; i < 3 ; ++i )
                pt[i] = cage_vertices[ vid[i] ];

            for( unsigned int i = 0 ; i < 3 ; ++i )
                N[i] = point3d::cross( pt[(i+1)%3] - eta , pt[(i+2)%3] - eta );

            for( unsigned int i = 0 ; i <= 2 ; ++i ) {
                w[i] = 0.0;
                for( unsigned int j = 0 ; j <= 2 ; ++j )
                    w[i] += theta[j] * point3d::dot( N[i] , N[j] ) / ( 2.0 * N[j].norm() );

                w[i] /= determinant;
            }

            sumWeights += ( w[0] + w[1] + w[2] );
            w_weights[ vid[0] ] += w[0];
            w_weights[ vid[1] ] += w[1];
            w_weights[ vid[2] ] += w[2];
        }
    }

    // CAGE QUADS :
    {
        unsigned int vid[4];
        point4d w;

        for( unsigned int q = 0 ; q < n_quads ; ++q ) {
            // the Norm is CCW :
            point3d pt[4];
            for( unsigned int i = 0 ; i < 4 ; ++i ) {
                vid[i] =  cage_quads[q][i];
                pt[i] = cage_vertices[ vid[i] ];
            }

            if( fabs(point3d::getSpanDeterminant(pt[1] - pt[0] , pt[2] - pt[0] , pt[3] - pt[0])) < epsilon
                    &&   fabs(point3d::getSpanDeterminant(pt[0] - eta , pt[1] - eta , pt[2] - eta)) < epsilon) {
                // then the quad is planar, eta is on the support plane.
                // i) if eta is INSIDE the quad, then output barycentric functions of the quad's vertices
                // ii) if eta is OUTSIDE the quad, just set weights w.r.t. the quad's vertices to 0 and go on with other faces
                // find out if eta is inside the quad:
                point3d const & Ntri = point3d::cross(pt[1] - pt[0] , pt[2] - pt[0]).direction();
                double wn = fabs(QuadUtilities::compute2DWindingNumberInQuad(eta,pt[0],pt[1],pt[2],pt[3],Ntri));
                bool eta_is_inside_the_quad = ( M_PI < wn ) && ( wn < 3.0*M_PI );
                if(eta_is_inside_the_quad) {
                    // compute barycentric coordinates:
                    double u,v;
                    QuadUtilities::computeFloaterBarycentricCoordinatesInPlanarQuad(eta,pt[0],pt[1],pt[2],pt[3],Ntri,u,v);

                    w_weights.clear();
                    weights.clear();
                    weights.resize( n_vertices , 0.0 );
                    weights[ vid[0] ] = (1.0-u)*(1.0-v);
                    weights[ vid[1] ] = u*(1.0-v);
                    weights[ vid[2] ] = u*v;
                    weights[ vid[3] ] = (1.0-u)*v;
                    return true;
                }
                else {
                    for( unsigned int i = 0 ; i < 4 ; ++i )
                        w[i] = 0.0;
                    continue;
                }
            }
            else {
                // 1) identify null space:
                point4d nullSpaceBasis = point4d(
                        point3d::getSpanDeterminant(pt[1] - eta , pt[3] - eta , pt[2] - eta) ,
                        point3d::getSpanDeterminant(pt[0] - eta , pt[2] - eta , pt[3] - eta) ,
                        point3d::getSpanDeterminant(pt[1] - eta , pt[0] - eta , pt[3] - eta) ,
                        point3d::getSpanDeterminant(pt[0] - eta , pt[1] - eta , pt[2] - eta)
                            );

                // nullSpaceBasis.normalize(); // THIS IS HIGHLY UNSTABLE !!!

                // check that the point is not on the bilinear quad:
                double alphaCand = nullSpaceBasis[0] + nullSpaceBasis[1] + nullSpaceBasis[2] + nullSpaceBasis[3];
                double uCand = (nullSpaceBasis[1] + nullSpaceBasis[2]) / alphaCand;
                double vCand = (nullSpaceBasis[2] + nullSpaceBasis[3]) / alphaCand;
                if( uCand <= 1.0  && vCand <= 1.0  &&  uCand >= 0.0  &&  vCand >= 0.0 ) {
                    point4d buvCand ( (1-uCand)*(1-vCand) , uCand*(1-vCand) , uCand*vCand , (1-uCand)*vCand );
                    if( (alphaCand * buvCand - nullSpaceBasis).sqrnorm() < epsilon*epsilon ) {
                        // Then the point is on the bilinear quad.
                        w_weights.clear();
                        weights.clear();
                        weights.resize( n_vertices , 0.0 );
                        weights[ vid[0] ] = buvCand[0];
                        weights[ vid[1] ] = buvCand[1];
                        weights[ vid[2] ] = buvCand[2];
                        weights[ vid[3] ] = buvCand[3];
                        return true;
                    }
                }

                // the point is NOT on the bilinear quad.
                // 2) solve full rank system:
                double theta[4] ;
                for( unsigned int i = 0 ; i < 4 ; ++i ) {
                    theta[i] = 2.0 * asin( ( u[ vid[ ( i ) % 4 ] ] - u[ vid[ ( i + 1 ) % 4 ] ] ).norm() / 2.0 );
                }

                point3d N[4];
                for( unsigned int i = 0 ; i < 4 ; ++i ) {
                    N[i] = point3d::cross( pt[(i+1)%4] - eta , pt[(i)%4] - eta );
                }

                point3d m_quad = -0.5 * ( theta[0]*N[0]/N[0].norm() + theta[1]*N[1]/N[1].norm() + theta[2]*N[2]/N[2].norm() + theta[3]*N[3]/N[3].norm() );

                w[0] = nullSpaceBasis[1] * point3d::getSpanDeterminant(m_quad,pt[2]-eta,pt[3]-eta)
                        - nullSpaceBasis[2] * point3d::getSpanDeterminant(m_quad,pt[1]-eta,pt[3]-eta)
                        + nullSpaceBasis[3] * point3d::getSpanDeterminant(m_quad,pt[1]-eta,pt[2]-eta);
                w[1] = - nullSpaceBasis[0] * point3d::getSpanDeterminant(m_quad,pt[2]-eta,pt[3]-eta)
                        - nullSpaceBasis[2] * point3d::getSpanDeterminant(pt[0]-eta,m_quad,pt[3]-eta)
                        + nullSpaceBasis[3] * point3d::getSpanDeterminant(pt[0]-eta,m_quad,pt[2]-eta);
                w[2] = - nullSpaceBasis[0] * point3d::getSpanDeterminant(pt[1]-eta,m_quad,pt[3]-eta)
                        + nullSpaceBasis[1] * point3d::getSpanDeterminant(pt[0]-eta,m_quad,pt[3]-eta)
                        + nullSpaceBasis[3] * point3d::getSpanDeterminant(pt[0]-eta,pt[1]-eta,m_quad);
                w[3] = - nullSpaceBasis[0] * point3d::getSpanDeterminant(pt[1]-eta,pt[2]-eta,m_quad)
                        + nullSpaceBasis[1] * point3d::getSpanDeterminant(pt[0]-eta,pt[2]-eta,m_quad)
                        - nullSpaceBasis[2] * point3d::getSpanDeterminant(pt[0]-eta,pt[1]-eta,m_quad);

                w /= nullSpaceBasis.sqrnorm();

                // WE NEED TO FIND THE VALUE TO ADD TO W (the component along nullSpaceBasis) :
                double lambda = 0.0;

                double uCenter , vCenter;
                QuadUtilities::smoothProjectOnQuad( eta , pt[0],pt[1],pt[2],pt[3], uCenter , vCenter );

                assert( uCenter >= 0.0 );
                assert( uCenter <= 1.0 );
                assert( vCenter >= 0.0 );
                assert( vCenter <= 1.0 );

                std::vector< double > uValues;
                std::vector< double > vValues;
                {
                    // SAMPLING FOR SUBMISSION : n = 4
                    unsigned int n = 4;

                    uValues.clear();
                    vValues.clear();

                    {
                        if( uCenter > 0.0 ){
                            for( unsigned int i = 0 ; i < n ; ++i ) {
                                double x = (double)(i) / (double)(n);
                                uValues.push_back(x * uCenter);
                            }
                        }
                        uValues.push_back(uCenter);
                        if( uCenter < 1.0 ){
                            for( unsigned int i = 0 ; i < n ; ++i ) {
                                double x = 1.0 - (double)(i + 1) / (double)(n);
                                uValues.push_back(uCenter * x + 1.0 * (1.0 - x));
                            }
                        }

                        if( vCenter > 0.0 ){
                            for( unsigned int i = 0 ; i < n ; ++i ) {
                                double x = (double)(i) / (double)(n);
                                vValues.push_back(x * vCenter);
                            }
                        }
                        vValues.push_back(vCenter);
                        if( vCenter < 1.0 ){
                            for( unsigned int i = 0 ; i < n ; ++i ) {
                                double x = (double)(i + 1) / (double)(n);
                                vValues.push_back(vCenter * (1.0 - x) + 1.0 * x);
                            }
                        }
                    }
                }

                // compute APPROXIMATE weights :
                point4d integratedBilinearWeights(0,0,0,0);
                for( unsigned int uIt = 0 ; uIt < uValues.size() ; ++uIt ) {
                    for( unsigned int vIt = 0 ; vIt < vValues.size() ; ++vIt ) {
                        double u = uValues[uIt];
                        double v = vValues[vIt];

                        point4d bilinearWeights ( (1-u)*(1-v) , u*(1-v) , u*v , (1-u)*v );
                        point3d const & Puv = QuadUtilities::bilinearInterpolation(pt[0],pt[1],pt[2],pt[3],u,v);

                        {
                            // ESTIMATE dB USING ELEMENTARY SOLID ANGLE AT THE POINT
                            double du = (uValues[std::min<unsigned int>(uIt+1 , uValues.size()-1)] - uValues[std::max<int>((int)(uIt)-1 , 0)])/2.0;
                            double dv = (vValues[std::min<unsigned int>(vIt+1 , vValues.size()-1)] - vValues[std::max<int>((int)(vIt)-1 , 0)])/2.0; // OK, triple checked
                            double dist = (Puv-eta).norm();
                            point3d const & Nuv = QuadUtilities::bilinearInterpolationNormal(pt[0],pt[1],pt[2],pt[3],u,v);
                            double dB = point3d::dot(Puv-eta , Nuv) * du * dv / (dist*dist*dist);

                            integratedBilinearWeights += (dB / dist) * bilinearWeights;
                        }
                    }
                }

                // NORM-CORRECT lambda:
                {
                    lambda = (point4d::dot( integratedBilinearWeights , nullSpaceBasis ) * w.sqrnorm()) /
                            ( nullSpaceBasis.sqrnorm() * point4d::dot( integratedBilinearWeights , w ) ); // with norm correction
                }

                w += lambda * nullSpaceBasis;
            }

            sumWeights += ( w[0] + w[1] + w[2] + w[3] );
            w_weights[ vid[0] ] += w[0];
            w_weights[ vid[1] ] += w[1];
            w_weights[ vid[2] ] += w[2];
            w_weights[ vid[3] ] += w[3];
        }
    }

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
        weights[v]  = w_weights[v] / sumWeights;

    return false;
}







template< class float_t , class point_t >
bool computeUnnormalizedMVCForOneTriangle(
        point_t const & eta ,
        point_t * tri_vertices ,
        float_t * w_weights )
{
    typedef typename point_t::type_t    T;

    T epsilon = 0.000000001;

    T d[3];
    point_t u[3];

    for( unsigned int v = 0 ; v < 3 ; ++v )
    {
        d[ v ] = ( eta - tri_vertices[ v ] ).norm();
        u[ v ] = ( tri_vertices[v] - eta ) / d[v];
    }

    T l[3]; T theta[3] ;

    {
        // the Norm is CCW :
        for( unsigned int i = 0 ; i <= 2 ; ++i ) {
            l[ i ] = ( u[ ( i + 1 ) % 3 ] - u[ ( i + 2 ) % 3 ] ).norm();
        }

        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            theta[i] = 2.0 * asin( l[i] / 2.0 );
        }

        T determinant = point_t::dot( tri_vertices[0] - eta ,
                point_t::cross( tri_vertices[1] - tri_vertices[0] , tri_vertices[2] - tri_vertices[0] ) );
        T sqrdist = determinant*determinant / (4 * point_t::cross( tri_vertices[1] - tri_vertices[0] ,
                tri_vertices[2] - tri_vertices[0] ).sqrnorm() );
        T dist = sqrt( (T)sqrdist );

        if( dist < epsilon )
        {
            // then the point eta lies on the support plane of the triangle
            T h = ( theta[0] + theta[1] + theta[2] ) / 2.0;
            if( M_PI - h < epsilon )
            {
                // then eta lies inside the triangle t
                for( unsigned int i = 0 ; i <= 2 ; ++i )
                {
                    w_weights[ i ] = sin( theta[ i ] ) * l[ (i+2) % 3 ] * l[ (i+1) % 3 ];
                }
                return true;
            }
            else {
                w_weights[0] = w_weights[1] = w_weights[2] = 0.0;
                return false;
            }
        }

        point_t N[3];

        for( unsigned int i = 0 ; i < 3 ; ++i )
            N[i] = point_t::cross( tri_vertices[(i+1)%3] - eta , tri_vertices[(i+2)%3] - eta );

        for( unsigned int i = 0 ; i <= 2 ; ++i )
        {
            w_weights[i] = 0.0;
            for( unsigned int j = 0 ; j <= 2 ; ++j )
                w_weights[i] += theta[j] * point_t::dot( N[i] , N[j] ) / ( 2.0 * N[j].norm() );

            w_weights[i] /= determinant;
        }
    }
    return false;
}








template< class int_t , class float_t >
bool computeCoordinatesWithTesselation(
        point3d eta ,
        std::vector< std::vector< int_t > > const & cage_triangles ,
        std::vector< std::vector< int_t > > const & cage_quads ,
        std::vector< point3d > const & cage_vertices_in ,
        std::vector< float_t > & weights,
        std::vector< float_t > & w_weights,
        unsigned int nTesselations )
{
    typedef  double    T;
    T epsilon = 0.00000001;


    std::vector< point3d > const & cage_vertices = cage_vertices_in;

    unsigned int n_vertices = cage_vertices.size();
    unsigned int n_triangles = cage_triangles.size();
    unsigned int n_quads = cage_quads.size();

    w_weights.clear();
    weights.clear();
    weights.resize( n_vertices , 0.0 );
    T sumWeights = 0.0;

    std::vector< T > d( n_vertices , 0.0 );
    std::vector< point3d > u( n_vertices );

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

    // CAGE TRIANGLES:
    {
        unsigned int vid[3];
        T l[3]; T theta[3] ; T w[3];

        for( unsigned int t = 0 ; t < n_triangles ; ++t )
        {
            // the Norm is CCW :
            for( unsigned int i = 0 ; i <= 2 ; ++i )
                vid[i] =  cage_triangles[t][i];

            for( unsigned int i = 0 ; i <= 2 ; ++i )
                l[ i ] = ( u[ vid[ ( i + 1 ) % 3 ] ] - u[ vid[ ( i + 2 ) % 3 ] ] ).norm();

            for( unsigned int i = 0 ; i <= 2 ; ++i )
                theta[i] = 2.0 * asin( l[i] / 2.0 );

            // test in original MVC paper: (they test if one angle psi is close to 0: it is "distance sensitive" in the sense that it does not
            // relate directly to the distance to the support plane of the triangle, and the more far away you go from the triangle, the worse it is)
            // In our experiments, it is actually not the good way to do it, as it increases significantly the errors we get in the computation of weights and derivatives,
            // especially when evaluating Hfx, Hfy, Hfz which can be of norm of the order of 10^3 instead of 0 (when specifying identity on the cage, see paper)

            // simple test we suggest:
            // the determinant of the basis is 2*area(T)*d( eta , support(T) ), we can directly test for the distance to support plane of the triangle to be minimum
            T determinant = point3d::dot( cage_vertices[vid[0]] - eta , point3d::cross( cage_vertices[vid[1]] - cage_vertices[vid[0]] , cage_vertices[vid[2]] - cage_vertices[vid[0]] ) );
            T sqrdist = determinant*determinant / (4 * point3d::cross( cage_vertices[vid[1]] - cage_vertices[vid[0]] , cage_vertices[vid[2]] - cage_vertices[vid[0]] ).sqrnorm() );
            T dist = sqrt( (T)sqrdist );


            if( dist < epsilon ) {
                // then the point eta lies on the support plane of the triangle
                T h = ( theta[0] + theta[1] + theta[2] ) / 2.0;
                if( M_PI - h < epsilon )
                {
                    // eta lies inside the triangle t , use 2d barycentric coordinates :
                    for( unsigned int i = 0 ; i <= 2 ; ++i )
                    {
                        w[ i ] = sin( theta[ i ] ) * l[ (i+2) % 3 ] * l[ (i+1) % 3 ];
                    }
                    sumWeights = w[0] + w[1] + w[2];

                    w_weights.clear();
                    weights.clear();
                    weights.resize( n_vertices , 0.0 );
                    weights[ vid[0] ] = w[0] / sumWeights;
                    weights[ vid[1] ] = w[1] / sumWeights;
                    weights[ vid[2] ] = w[2] / sumWeights;
                    return true;
                }
            }

            point3d pt[3] , N[3];

            for( unsigned int i = 0 ; i < 3 ; ++i )
                pt[i] = cage_vertices[ cage_triangles[t][i] ];
            for( unsigned int i = 0 ; i < 3 ; ++i )
                N[i] = point3d::cross( pt[(i+1)%3] - eta , pt[(i+2)%3] - eta );

            for( unsigned int i = 0 ; i <= 2 ; ++i )
            {
                w[i] = 0.0;
                for( unsigned int j = 0 ; j <= 2 ; ++j )
                    w[i] += theta[j] * point3d::dot( N[i] , N[j] ) / ( 2.0 * N[j].norm() );

                w[i] /= determinant;
            }

            sumWeights += ( w[0] + w[1] + w[2] );
            w_weights[ vid[0] ] += w[0];
            w_weights[ vid[1] ] += w[1];
            w_weights[ vid[2] ] += w[2];
        }
    }

    // CAGE QUADS :
    {
        unsigned int vid[4];
        point3d qc[4];

        for( unsigned int q = 0 ; q < n_quads ; ++q ) {
            point4d w(0,0,0,0);

            for( unsigned int i = 0 ; i <= 3 ; ++i )
                vid[i] =  cage_quads[q][i];

            for( unsigned int i = 0 ; i <= 3 ; ++i )
                qc[i] =  cage_vertices[vid[i]];

            for( unsigned int lTessU = 0 ; lTessU < nTesselations ; ++lTessU ) {
                for( unsigned int lTessV = 0 ; lTessV < nTesselations ; ++lTessV ) {
                    double uMin = (double)( lTessU ) / (double)( nTesselations ),
                            uMax = (double)( lTessU+1 ) / (double)( nTesselations ),
                            vMin = (double)( lTessV ) / (double)( nTesselations ),
                            vMax = (double)( lTessV+1 ) / (double)( nTesselations );
                    point4d buv ( (1-uMin)*(1-vMin) , uMin*(1-vMin) , uMin*vMin , (1-uMin)*vMin );
                    point4d bUv ( (1-uMax)*(1-vMin) , uMax*(1-vMin) , uMax*vMin , (1-uMax)*vMin );
                    point4d bUV ( (1-uMax)*(1-vMax) , uMax*(1-vMax) , uMax*vMax , (1-uMax)*vMax );
                    point4d buV ( (1-uMin)*(1-vMax) , uMin*(1-vMax) , uMin*vMax , (1-uMin)*vMax );

                    point3d quv = buv[0] * qc[0] + buv[1] * qc[1] + buv[2] * qc[2] + buv[3] * qc[3];

                    point3d quV = buV[0] * qc[0] + buV[1] * qc[1] + buV[2] * qc[2] + buV[3] * qc[3];

                    point3d qUV = bUV[0] * qc[0] + bUV[1] * qc[1] + bUV[2] * qc[2] + bUV[3] * qc[3];

                    point3d qUv = bUv[0] * qc[0] + bUv[1] * qc[1] + bUv[2] * qc[2] + bUv[3] * qc[3];

                    point3d triTess[3];
                    double triTessW[3];
                    {
                        triTess[0] = quv; triTess[1] = qUv; triTess[2] = qUV;
                        if( computeUnnormalizedMVCForOneTriangle(eta, triTess , triTessW) ) {
                            weights.clear();
                            weights.resize( n_vertices , 0.0 );
                            weights[vid[0]] = buv[0] * triTessW[0] +  bUv[0] * triTessW[1] +  bUV[0] * triTessW[2];
                            weights[vid[1]] = buv[1] * triTessW[0] +  bUv[1] * triTessW[1] +  bUV[1] * triTessW[2];
                            weights[vid[2]] = buv[2] * triTessW[0] +  bUv[2] * triTessW[1] +  bUV[2] * triTessW[2];
                            weights[vid[3]] = buv[3] * triTessW[0] +  bUv[3] * triTessW[1] +  bUV[3] * triTessW[2];

                            w_weights.clear();
                            return true;
                        }
                        else {
                            w[0] += buv[0] * triTessW[0] +  bUv[0] * triTessW[1] +  bUV[0] * triTessW[2];
                            w[1] += buv[1] * triTessW[0] +  bUv[1] * triTessW[1] +  bUV[1] * triTessW[2];
                            w[2] += buv[2] * triTessW[0] +  bUv[2] * triTessW[1] +  bUV[2] * triTessW[2];
                            w[3] += buv[3] * triTessW[0] +  bUv[3] * triTessW[1] +  bUV[3] * triTessW[2];
                        }
                    }
                    {
                        triTess[0] = quv; triTess[1] = qUV; triTess[2] = quV;
                        if( computeUnnormalizedMVCForOneTriangle(eta, triTess , triTessW) ) {
                            weights.clear();
                            weights.resize( n_vertices , 0.0 );
                            weights[vid[0]] = buv[0] * triTessW[0] +  bUV[0] * triTessW[1] +  buV[0] * triTessW[2];
                            weights[vid[1]] = buv[1] * triTessW[0] +  bUV[1] * triTessW[1] +  buV[1] * triTessW[2];
                            weights[vid[2]] = buv[2] * triTessW[0] +  bUV[2] * triTessW[1] +  buV[2] * triTessW[2];
                            weights[vid[3]] = buv[3] * triTessW[0] +  bUV[3] * triTessW[1] +  buV[3] * triTessW[2];

                            w_weights.clear();
                            return true;
                        }
                        else {
                            w[0] += buv[0] * triTessW[0] +  bUV[0] * triTessW[1] +  buV[0] * triTessW[2];
                            w[1] += buv[1] * triTessW[0] +  bUV[1] * triTessW[1] +  buV[1] * triTessW[2];
                            w[2] += buv[2] * triTessW[0] +  bUV[2] * triTessW[1] +  buV[2] * triTessW[2];
                            w[3] += buv[3] * triTessW[0] +  bUV[3] * triTessW[1] +  buV[3] * triTessW[2];
                        }
                    }
                }
            }

            sumWeights += ( w[0] + w[1] + w[2] + w[3] );
            w_weights[ vid[0] ] += w[0];
            w_weights[ vid[1] ] += w[1];
            w_weights[ vid[2] ] += w[2];
            w_weights[ vid[3] ] += w[3];
        }
    }

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
        weights[v]  = w_weights[v] / sumWeights;

    return false;
}

















template< class int_t , class float_t , class point3d >
inline
bool computeCoordinates(
        point3d const & eta ,
        std::vector< std::vector< int_t > > const & cage_triangles ,
        std::vector< std::vector< int_t > > const & cage_quads ,
        std::vector< point3d > const & cage_vertices ,
        std::vector< float_t > & weights)
{
    std::vector< float_t > w_weights;
     return computeCoordinates(eta,cage_triangles,cage_quads,cage_vertices,weights,w_weights);
    // return computeCoordinatesWithTesselation(eta,cage_triangles,cage_quads,cage_vertices,weights,w_weights , 4);
}








}


#endif // QMVC_H
