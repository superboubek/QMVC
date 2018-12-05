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

#ifndef GC_H
#define GC_H


//-----------------------------------   ASSOCIATED REFERENCES   ---------------------------------------------//
//                                                                                                           //
// "Green Coordinates". Lipman et al. 2008, ACM Siggraph:                                                    //
// this paper presents the use of Green coordinates for deforming 3D meshes                                  //
//                                                                                                           //
// "Variational harmonic maps for space deformation". Ben-Chen et al. 2009, ACM Siggraph                     //
// this paper contains an alternative computation of the coordinates, along with gradients and Hessians,     //
// based on a work of Urago (see annex of the paper)                                                         //
//                                                                                                           //
// You can find these expressions in URAGO 2000 :                                                            //
// "Analytical integrals of fundamental solution of three-dimensional Laplace equation and their gradients"  //
//-----------------------------------------------------------------------------------------------------------//




#include <vector>
#include <cmath>







namespace GreenCoordinates3D
{

template< class point_t >
double areaRatioSqrt( point_t const & _u , point_t const & _v , point_t const & u , point_t const & v ) {
    return sqrt( point_t::cross(u,v).norm() / point_t::cross(_u,_v).norm() );
}

template< class point_t >
double confFactor( point_t const & _u , point_t const & _v , point_t const & u , point_t const & v ) {
    return sqrt(
                ( (_u.sqrnorm()) * (v.sqrnorm())  -  2.0*( point_t::dot(_u , _v) )*( point_t::dot(u , v) ) + (u.sqrnorm()) * (_v.sqrnorm()) ) /
                ( 2.0 * point_t::cross( _u , _v ).sqrnorm() )
                ); // v1 : As in Lipman's paper
}


template< class point_t >
class GreenScalingFactor
{
    point_t _u;
    point_t _v;
    // u is the first halfedge of the triangle (v1 - v0) , and v is the second one (v2 - v1).

    typedef typename point_t::type_t T;
    T _s;

public:
    void setOriginalEdges( point_t const & u , point_t const & v ) { _u = u; _v = v; }
    GreenScalingFactor() {} // useless
    GreenScalingFactor( point_t const & u , point_t const & v ) : _u( u ) , _v( v ) , _s( 1.0 ) {}
    GreenScalingFactor( point_t const & v0 , point_t const & v1 , point_t const & v2 ) : _u( v1-v0 ) , _v( v2-v0 ) , _s( 1.0 ) {}

    void computeScalingFactor( point_t const & u , point_t const & v )
    {
        _s = confFactor(_u,_v,u,v);
    }
    void computeScalingFactor( point_t const & v0 , point_t const & v1 , point_t const & v2 )
    {
        computeScalingFactor(v1-v0 , v2-v0);
    }


    void computeScalingFactorBasedOnAreas(  point_t const & u , point_t const & v  )
    {
        _s = areaRatioSqrt(_u,_v,u,v);
    }
    void computeScalingFactorBasedOnAreas( point_t const & v0 , point_t const & v1 , point_t const & v2 )
    {
        computeScalingFactorBasedOnAreas(v1-v0 , v2-v0);
    }


    T scalingFactor() const { return _s; }
};


template< class point_t >
class GreenScalingFactorQuad {

    point_t v0 ,v1 ,v2 ,v3;

    typedef typename point_t::type_t T;
    T _s0, _s1, _s2, _s3;

public:
    GreenScalingFactorQuad() {} // useless
    GreenScalingFactorQuad( point_t const & u0 , point_t const & u1 , point_t const & u2 , point_t const & u3 ) :
        v0( u0 ) , v1( u1 ) , v2( u2 ) , v3( u3 ) , _s0( 1.0 ) , _s1( 1.0 ) , _s2( 1.0 ) , _s3( 1.0 ) {}


    void computeScalingFactor( point_t const & u0 , point_t const & u1 , point_t const & u2 , point_t const & u3 )
    {
        _s0 = confFactor( v1-v0,v3-v0 , u1-u0 , u3-u0 );
        _s1 = confFactor( v2-v1,v0-v1 , u2-u1 , u0-u1 );
        _s2 = confFactor( v3-v2,v1-v2 , u3-u2 , u1-u2 );
        _s3 = confFactor( v0-v3,v2-v3 , u0-u3 , u2-u3 );
    }
    void computeScalingFactorBasedOnAreas( point_t const & u0 , point_t const & u1 , point_t const & u2 , point_t const & u3 )
    {
        _s0 = areaRatioSqrt( v1-v0,v3-v0 , u1-u0 , u3-u0 );
        _s1 = areaRatioSqrt( v2-v1,v0-v1 , u2-u1 , u0-u1 );
        _s2 = areaRatioSqrt( v3-v2,v1-v2 , u3-u2 , u1-u2 );
        _s3 = areaRatioSqrt( v0-v3,v2-v3 , u0-u3 , u2-u3 );
    }

    T scalingFactor(unsigned int i) const {
        if(i==0)return _s0;
        if(i==1)return _s1;
        if(i==2)return _s2;
        if(i==3)return _s3;
        assert(0);
        return 0.0;
    }
};






namespace URAGO
{
// return the signed solid angle of the tetrahedron supported by the three vectors a,b, and c:
template< class point_t >
double get_signed_solid_angle(point_t const & a, point_t const & b, point_t const & c) {
    typedef typename point_t::type_t    T;
    T det = point_t::dot(a,point_t::cross(b,c));
    if( fabs(det) < 0.0000000001 ) // then you're on the limit case where you cover half the sphere
        return 2.0 * M_PI; // that is particularly shitty, because the sign is difficult to estimate...

    T al = a.norm(),   bl = b.norm(),   cl = c.norm();

    T div = al*bl*cl + point_t::dot(a,b)*cl + point_t::dot(a,c)*bl + point_t::dot(b,c)*al;
    T at = atan2( fabs(det) , div );
    if(at < 0) at += M_PI; // If det>0 && div<0 atan2 returns < 0, so add pi.
    T omega = 2.0 * at;

    if(det > 0.0) return omega;
    return -omega;
}






template< class float_t , class point_t >
void computePhiAndPsiForOneTriangle( point_t const & eta ,
                                     point_t * tri_vertices , // an array of 3 points
                                     float_t * phi , // an array of 3 floats
                                     float_t & psi) {
    typedef typename point_t::type_t    T;
    point_t Nt = point_t::cross( tri_vertices[1]-tri_vertices[0], tri_vertices[2]-tri_vertices[0]);
    T NtNorm = Nt.norm();
    T At = NtNorm / 2.0;
    Nt /= NtNorm;

    psi = 0.0;
    for( unsigned int v = 0 ; v < 3 ; ++v ) phi[v] = 0.0;

    point_t e[3];    T e_norm[3];   point_t e_normalized[3];    T R[3];    point_t d[3];    T d_norm[3];     T C[3];     point_t J[3];
    for( unsigned int v = 0 ; v < 3 ; ++v ) e[v] = tri_vertices[v] - eta;
    for( unsigned int v = 0 ; v < 3 ; ++v ) e_norm[v] =  e[v].norm();
    for( unsigned int v = 0 ; v < 3 ; ++v ) e_normalized[v] = e[v] / e_norm[v];

    T signed_solid_angle = get_signed_solid_angle (e_normalized[0], e_normalized[1], e_normalized[2]) / (4.f * M_PI);
    T signed_volume = point_t::dot( point_t::cross(e[0],e[1]) , e[2] ) / 6.0;

    for( unsigned int v = 0 ; v < 3 ; ++v ) R[v] = e_norm[(v+1)%3] + e_norm[(v+2)%3];
    for( unsigned int v = 0 ; v < 3 ; ++v ) d[v] = tri_vertices[(v+1)%3] - tri_vertices[(v+2)%3];
    for( unsigned int v = 0 ; v < 3 ; ++v ) d_norm[v] =  d[v].norm();
    for( unsigned int v = 0 ; v < 3 ; ++v ) C[v] = log( (R[v] + d_norm[v]) / (R[v] - d_norm[v]) ) / (4.0 * M_PI * d_norm[v]);

    point_t Pt( - signed_solid_angle * Nt );
    for( unsigned int v = 0 ; v < 3 ; ++v ) Pt += point_t::cross( Nt , C[v]*d[v] );
    for( unsigned int v = 0 ; v < 3 ; ++v ) J[v] = point_t::cross( e[(v+2)%3] , e[(v+1)%3] );

    psi = - 3.0 * signed_solid_angle * signed_volume / At ;
    for( unsigned int v = 0 ; v < 3 ; ++v ) psi -= C[v]* point_t::dot(J[v],Nt);
    for( unsigned int v = 0 ; v < 3 ; ++v ) phi[ v ] += point_t::dot(Pt , J[v]) / (2.0 * At);
}








template< class int_t , class float_t , class point_t >
void computeCoordinates(
        point_t const & eta ,
        std::vector< std::vector< int_t > > const & cage_triangles ,
        std::vector< point_t > const & cage_vertices ,
        std::vector< float_t > & _VC_phi ,
        std::vector< float_t > & _TC_psi) {
    typedef typename point_t::type_t    T;

    _VC_phi.clear();  _VC_phi.resize( cage_vertices.size() , 0.0 );
    _TC_psi.clear();  _TC_psi.resize( cage_triangles.size() , 0.0 );

    // iterate over the triangles:
    for( unsigned int t = 0 ; t < cage_triangles.size() ; ++t ) {
        std::vector<int_t> const & tri = cage_triangles[t];
        point_t tri_verts[3] = {cage_vertices[tri[0]] , cage_vertices[tri[1]] , cage_vertices[tri[2]]};
        T psi;  T phi[3];

        computePhiAndPsiForOneTriangle( eta , tri_verts , phi , psi );

        _TC_psi[t] = psi;
        for( unsigned int v = 0 ; v < 3 ; ++v ) _VC_phi[ tri[v] ] += phi[v];
    }
}




}


}



#endif // GC_H
