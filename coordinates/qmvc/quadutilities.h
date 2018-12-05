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

#ifndef QUADUTILITIES_H
#define QUADUTILITIES_H



namespace QuadUtilities {
template< class point_t >
point_t triMeanVector(point_t const & p , point_t const & p0 , point_t const & p1 , point_t const & p2)
{
    point_t pttri[3] = {p0,p1,p2};
    point_t utri[3];
    for( unsigned int i = 0 ; i < 3 ; ++i ) {
        utri[i] = (pttri[i] - p).direction(true);
        assert( ! utri[i].isnan() );
    }

    double thetatri[3];
    for( unsigned int i = 0 ; i < 3 ; ++i ){
        assert( ! std::isnan(( utri[ ( i ) % 3 ] - utri[ ( i + 1 ) % 3 ] ).norm()) );
        thetatri[i] = 2.0 * asin( std::min<double>(1.0 , std::max<double>( -1.0 , ( utri[ ( i ) % 3 ] - utri[ ( i + 1 ) % 3 ] ).norm() / 2.0 ) ) );
        assert( ! std::isnan(thetatri[i]) );
    }

    point_t Ntri[3];
    for( unsigned int i = 0 ; i < 3 ; ++i ) {
        Ntri[i] = point_t::cross( pttri[(i+1)%3] - p , pttri[(i)%3] - p );
        assert( ! Ntri[i].isnan() );
    }

    point_t m_tri = -0.5 * ( thetatri[0]*Ntri[0].direction(true) + thetatri[1]*Ntri[1].direction(true) + thetatri[2]*Ntri[2].direction(true) );
    return m_tri;
}



template< class point_t >
point_t quadMeanVector(point_t const & p , point_t const & p0 , point_t const & p1 , point_t const & p2 , point_t const & p3)
{
    point_t mv(0,0,0);
    point_t t0,t1,t2 , pc = (p0+p1+p2+p3)/4.0;
    t0 = p0; t1 = p1; t2 = pc;
    {
        mv += triMeanVector( p , t0 , t1 , t2 );
    }
    t0 = p1; t1 = p2; t2 = pc;
    {
        mv += triMeanVector( p , t0 , t1 , t2 );
    }
    t0 = p2; t1 = p3; t2 = pc;
    {
        mv += triMeanVector( p , t0 , t1 , t2 );
    }
    t0 = p3; t1 = p0; t2 = pc;
    {
        mv += triMeanVector( p , t0 , t1 , t2 );
    }
    return mv;
}





template< class point_t >
double computeClosestPointOnLineFromLine( point_t const & A , point_t const & B , point_t const & eta , point_t const & d ) {
    point_t const & e = B - A;
    double eDotD = point_t::dot(e,d);
    double dDotD = point_t::dot(d,d);
    double eDotE = point_t::dot(e,e);

    double alpha = eDotE - eDotD*eDotD / dDotD;
    double beta = point_t::dot( A - eta , e )  -  point_t::dot( A - eta , d ) * eDotD / dDotD;

    return - beta / alpha;
}

template< class point_t >
void computeClosestPointOnQuadFromLine( point_t const & q0 , point_t const & q1 , point_t const & q2 , point_t const & q3,
                                        point_t const & origin , point_t const & direction ,
                                        double & uProj , double & vProj , bool ForceIntoUnitSquare = true ) {
    double epsilon = 0.0;

    double aV = point_t::dot( point_t::cross( q3-q0,q2-q1 ) , direction );
    double aU = point_t::dot( point_t::cross( q1-q0,q2-q3 ) , direction );

    if( fabs(aV)  >  fabs(aU) ){
        // solve for v :
        double a = point_t::dot( point_t::cross( q3-q0,q2-q1 ) , direction );
        double b = point_t::dot( point_t::cross( q0-origin,q2-q1 ) + point_t::cross( q3-q0,q1-origin ) , direction );
        double c = point_t::dot( point_t::cross( q0-origin,q1-origin ) , direction );
        if(fabs(a) > epsilon) {
            double Delta = b*b - 4*a*c;
            if(Delta >= 0.0) {
                double vPlus = (-b + sqrt(Delta)) / (2*a);
                double vMinus = (-b - sqrt(Delta)) / (2*a);
                if( fabs(vMinus - 0.5)  <  fabs(vPlus - 0.5) )
                    vProj = vMinus;
                else
                    vProj = vPlus;
            }
        }
        else {
            vProj = -c / b;
        }

        if(ForceIntoUnitSquare) vProj = std::min<double>( 1.0 , std::max<double>( 0.0 , vProj ) );

        point_t AV = (1.0 - vProj) * q0 + vProj * q3;
        point_t BV = (1.0 - vProj) * q1 + vProj * q2;

        uProj = computeClosestPointOnLineFromLine( AV , BV , origin , direction );
        if(ForceIntoUnitSquare) uProj = std::min<double>( 1.0 , std::max<double>( 0.0 , uProj ) );
    }
    else {
        // solve for u :
        double a = point_t::dot( point_t::cross( q1-q0,q2-q3 ) , direction );
        double b = point_t::dot( point_t::cross( q0-origin,q2-q3 ) + point_t::cross( q1-q0,q3-origin ) , direction );
        double c = point_t::dot( point_t::cross( q0-origin,q3-origin ) , direction );
        if(fabs(a) > epsilon) {
            double Delta = b*b - 4*a*c;
            if(Delta >= 0.0) {
                double uPlus = (-b + sqrt(Delta)) / (2*a);
                double uMinus = (-b - sqrt(Delta)) / (2*a);
                if( fabs(uMinus - 0.5)  <  fabs(uPlus - 0.5) )
                    uProj = uMinus;
                else
                    uProj = uPlus;
            }
        }
        else {
            uProj = -c / b;
        }

        if(ForceIntoUnitSquare) uProj = std::min<double>( 1.0 , std::max<double>( 0.0 , uProj ) );

        point_t AU = (1.0 - uProj) * q0 + uProj * q1;
        point_t BU = (1.0 - uProj) * q3 + uProj * q2;

        vProj = computeClosestPointOnLineFromLine( AU , BU , origin , direction );
        if(ForceIntoUnitSquare) vProj = std::min<double>( 1.0 , std::max<double>( 0.0 , vProj ) );
    }
}



template< class point_t >
bool isInConvexHull(point_t const & p , point_t const & p0 , point_t const & p1 , point_t const & p2 , point_t const & p3) {
    double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
    point_t t0,t1,t2;
    t0 = p0; t1 = p1; t2 = p3;
    {
        point_t nTri = point_t::cross( t1 - t0 , t2 - t0 );
        s0 = point_t::dot( t0 - p , nTri );
    }
    t0 = p1; t1 = p2; t2 = p3;
    {
        point_t nTri = point_t::cross( t1 - t0 , t2 - t0 );
        s1 = point_t::dot( t0 - p , nTri );
    }
    t0 = p0; t1 = p2; t2 = p1;
    {
        point_t nTri = point_t::cross( t1 - t0 , t2 - t0 );
        s2 = point_t::dot( t0 - p , nTri );
    }
    t0 = p0; t1 = p3; t2 = p2;
    {
        point_t nTri = point_t::cross( t1 - t0 , t2 - t0 );
        s3 = point_t::dot( t0 - p , nTri );
    }
    return (s0 >= 0.0 && s1 >= 0.0 && s2 >= 0.0 && s3 >= 0.0)  ||  (s0 <= 0.0 && s1 <= 0.0 && s2 <= 0.0 && s3 <= 0.0);
}


template< class point_t >
point_t bilinearInterpolation( point_t const & q0 , point_t const & q1 , point_t const & q2, point_t const & q3 , double u , double v ){
    return (1.0-u)*(1.0-v)*q0 + u*(1.0-v)*q1 + u*v*q2 + (1.0-u)*v*q3;
}



template< class point_t >
point_t bilinearInterpolationNormal( point_t const & q0 , point_t const & q1 , point_t const & q2, point_t const & q3 , double u , double v ){
    point_t const & N0 = point_t::cross(q1-q0 , q3-q0);
    point_t const & N1 = point_t::cross(q2-q1 , q0-q1);
    point_t const & N2 = point_t::cross(q3-q2 , q1-q2);
    point_t const & N3 = point_t::cross(q0-q3 , q2-q3);
    return (1.0-u)*(1.0-v)*N0 + u*(1.0-v)*N1 + u*v*N2 + (1.0-u)*v*N3;
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

        for( unsigned int i = 0 ; i <= 2 ; ++i ) {
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


template< class point_t >
point_t smoothProjectInsideTet( point_t const & eta , point_t const & q0 , point_t const & q1 , point_t const & q2 , point_t const & q3 ) {
    point_t proj(0,0,0);
    point_t tri[3];
    double sumWeights = 0.0;
    double phi[3];

    tri[0] = q0; tri[1] = q1; tri[2] = q2;
    {
        computeUnnormalizedMVCForOneTriangle(eta , tri , phi);
        sumWeights += fabs(phi[0]) + fabs(phi[1]) + fabs(phi[2]);
        proj += fabs(phi[0])*tri[0] + fabs(phi[1])*tri[1] + fabs(phi[2])*tri[2];
    }

    tri[0] = q0; tri[1] = q1; tri[2] = q3;
    {
        computeUnnormalizedMVCForOneTriangle(eta , tri , phi);
        sumWeights += fabs(phi[0]) + fabs(phi[1]) + fabs(phi[2]);
        proj += fabs(phi[0])*tri[0] + fabs(phi[1])*tri[1] + fabs(phi[2])*tri[2];
    }

    tri[0] = q0; tri[1] = q3; tri[2] = q2;
    {
        computeUnnormalizedMVCForOneTriangle(eta , tri , phi);
        sumWeights += fabs(phi[0]) + fabs(phi[1]) + fabs(phi[2]);
        proj += fabs(phi[0])*tri[0] + fabs(phi[1])*tri[1] + fabs(phi[2])*tri[2];
    }

    tri[0] = q3; tri[1] = q1; tri[2] = q2;
    {
        computeUnnormalizedMVCForOneTriangle(eta , tri , phi);
        sumWeights += fabs(phi[0]) + fabs(phi[1]) + fabs(phi[2]);
        proj += fabs(phi[0])*tri[0] + fabs(phi[1])*tri[1] + fabs(phi[2])*tri[2];
    }

    return proj / sumWeights;
}



template< class point_t >
point_t smoothProjectInsideTet_second( point_t const & eta , point_t const & q0 , point_t const & q1 , point_t const & q2 , point_t const & q3 ) {
    point_t proj(0,0,0);
    point_t tri[3];
    double sumWeights = 0.0;
    double phi[3];

    tri[0] = q0; tri[1] = q1; tri[2] = q2;
    {
        computeUnnormalizedMVCForOneTriangle(eta , tri , phi);
        double dist_to_t = point_t::dot( tri[0] - eta , point_t::cross( tri[1] - tri[0] , tri[2] - tri[0] ).direction() );
        if( fabs(dist_to_t) != 0.0 ) {
            phi[0] /= dist_to_t;
            phi[1] /= dist_to_t;
            phi[2] /= dist_to_t;
        }
        sumWeights += phi[0] + phi[1] + phi[2];
        proj += phi[0]*tri[0] + phi[1]*tri[1] + phi[2]*tri[2];
    }

    tri[0] = q0; tri[1] = q1; tri[2] = q3;
    {
        computeUnnormalizedMVCForOneTriangle(eta , tri , phi);
        double dist_to_t = point_t::dot( tri[0] - eta , point_t::cross( tri[1] - tri[0] , tri[2] - tri[0] ).direction() );
        if( fabs(dist_to_t) != 0.0 ) {
            phi[0] /= dist_to_t;
            phi[1] /= dist_to_t;
            phi[2] /= dist_to_t;
        }
        sumWeights += phi[0] + phi[1] + phi[2];
        proj += phi[0]*tri[0] + phi[1]*tri[1] + phi[2]*tri[2];
    }

    tri[0] = q0; tri[1] = q3; tri[2] = q2;
    {
        computeUnnormalizedMVCForOneTriangle(eta , tri , phi);
        double dist_to_t = point_t::dot( tri[0] - eta , point_t::cross( tri[1] - tri[0] , tri[2] - tri[0] ).direction() );
        if( fabs(dist_to_t) != 0.0 ) {
            phi[0] /= dist_to_t;
            phi[1] /= dist_to_t;
            phi[2] /= dist_to_t;
        }
        sumWeights += phi[0] + phi[1] + phi[2];
        proj += phi[0]*tri[0] + phi[1]*tri[1] + phi[2]*tri[2];
    }

    tri[0] = q3; tri[1] = q1; tri[2] = q2;
    {
        computeUnnormalizedMVCForOneTriangle(eta , tri , phi);
        double dist_to_t = point_t::dot( tri[0] - eta , point_t::cross( tri[1] - tri[0] , tri[2] - tri[0] ).direction() );
        if( fabs(dist_to_t) != 0.0 ) {
            phi[0] /= dist_to_t;
            phi[1] /= dist_to_t;
            phi[2] /= dist_to_t;
        }
        sumWeights += phi[0] + phi[1] + phi[2];
        proj += phi[0]*tri[0] + phi[1]*tri[1] + phi[2]*tri[2];
    }

    return proj / sumWeights;
}








template< class point_t >
void smoothProjectOnQuad( point_t const & eta , point_t const & q0 , point_t const & q1 , point_t const & q2 , point_t const & q3 , double & uProject , double & vProject ) {
    point_t pToProject = eta , m;
    bool pointIsStrictlyInConvexHull = true;

    if( ! isInConvexHull(pToProject,q0,q1,q2,q3) ) {
       {
           pToProject = smoothProjectInsideTet_second(eta , q0 , q1 , q2 , q3); // either smoothProjectInsideTet or smoothProjectInsideTet_second... try them
       }
    }

    {
        m = quadMeanVector(pToProject,q0,q1,q2,q3).direction();

        if(pointIsStrictlyInConvexHull)
            assert( ! m.isnan() &&   "m.isnan() : bilMeanVector seems to be problematic for a point strictly inside the convex hull"  );
        else
            assert( ! m.isnan() &&   "m.isnan() : bilMeanVector seems to be problematic for a point ON the convex hull"  );

        computeClosestPointOnQuadFromLine(q0 , q1 , q2 , q3 , pToProject , m , uProject , vProject );
    }
}







template< class point_t >
double compute2DWindingNumberInQuad( point_t const & eta , point_t const & q0 , point_t const & q1 , point_t const & q2 , point_t const & q3 , point_t const & Ntri ) {
    point_t const & u0 = (q0-eta).direction();
    point_t const & u1 = (q1-eta).direction();
    point_t const & u2 = (q2-eta).direction();
    point_t const & u3 = (q3-eta).direction();

    double t0 = asin(std::min<double>(1.0, std::max<double>(-1.0 , point_t::dot(Ntri,point_t::cross(u0,u1)))));
    double t1 = asin(std::min<double>(1.0, std::max<double>(-1.0 , point_t::dot(Ntri,point_t::cross(u1,u2)))));
    double t2 = asin(std::min<double>(1.0, std::max<double>(-1.0 , point_t::dot(Ntri,point_t::cross(u2,u3)))));
    double t3 = asin(std::min<double>(1.0, std::max<double>(-1.0 , point_t::dot(Ntri,point_t::cross(u3,u0)))));
    return t0+t1+t2+t3;
}





template< class point_t >
void computeFloaterBarycentricCoordinatesInPlanarQuad( point_t const & eta , point_t const & q0 , point_t const & q1 , point_t const & q2 , point_t const & q3 , point_t const & nq ,
                                                double & u , double & v ) {
    double A0 = point_t::dot(nq , point_t::cross(q0-eta,q1-eta));
    double A1 = point_t::dot(nq , point_t::cross(q2-eta,q3-eta));
    double A2 = point_t::dot(nq , point_t::cross(q3-eta,q0-eta));
    double A3 = point_t::dot(nq , point_t::cross(q0-eta,q1-eta));
    double B0 = point_t::dot(nq , point_t::cross(q3-eta,q1-eta));
    double B1 = point_t::dot(nq , point_t::cross(q1-eta,q3-eta));
    // double B2 = point_t::dot(nq , point_t::cross(q2-eta,q0-eta));
    double B3 = point_t::dot(nq , point_t::cross(q3-eta,q1-eta));

    double D = std::max<double>(0.0 , B0*B0 + B1*B1 + 2*A0*A2 + 2*A1*A3);
    double D_sqrt = sqrt(D);

    double E0  = 2*A0 - B0 - B1 + D_sqrt;
    double E3  = 2*A3 - B3 - B0 + D_sqrt;

    u = 2*A3 / E3;
    v = 2*A0 / E0;
}


}


#endif // QUADUTILITIES_H

