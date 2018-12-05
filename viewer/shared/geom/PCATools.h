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
#ifndef PCATOOLS_H
#define PCATOOLS_H

#include <vector>
#include <cassert>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>


template< typename point_t , typename normal_t >
class basic_symetric_matrix3f
{
public:
    double vals[9];

    basic_symetric_matrix3f()
    {
        vals[0] = 0.f;
        vals[1] = 0.f;
        vals[2] = 0.f;
        vals[4] = 0.f;
        vals[5] = 0.f;
        vals[8] = 0.f;
    }
    void clear()
    {
        vals[0] = 0.f;
        vals[1] = 0.f;
        vals[2] = 0.f;
        vals[4] = 0.f;
        vals[5] = 0.f;
        vals[8] = 0.f;
    }

    void diagonalize( float & lambdaX_ , float & lambdaY_ , float & lambdaZ_ , normal_t & RepX_ , normal_t & RepY_ , normal_t & RepZ_)
    {
        vals[3] = vals[1];
        vals[6] = vals[2];
        vals[7] = vals[5];

        gsl_matrix_view m = gsl_matrix_view_array (vals, 3, 3);

        gsl_vector *eval = gsl_vector_alloc (3);
        gsl_matrix *evec = gsl_matrix_alloc (3, 3);
        gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
        gsl_eigen_symmv (&m.matrix, eval, evec, w);
        gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

        lambdaX_ = gsl_vector_get (eval, 0);
        gsl_vector_view evec_i = gsl_matrix_column (evec, 0);
        RepX_ = normal_t(
                    gsl_vector_get (&evec_i.vector, 0),
                    gsl_vector_get (&evec_i.vector, 1),
                    gsl_vector_get (&evec_i.vector, 2)
                    );

        lambdaY_ = gsl_vector_get (eval, 1);
        evec_i = gsl_matrix_column (evec, 1);
        RepY_ = normal_t(
                    gsl_vector_get (&evec_i.vector, 0),
                    gsl_vector_get (&evec_i.vector, 1),
                    gsl_vector_get (&evec_i.vector, 2)
                    );

        lambdaZ_ = gsl_vector_get (eval, 2);
        evec_i = gsl_matrix_column (evec, 2);
        RepZ_ = normal_t(
                    gsl_vector_get (&evec_i.vector, 0),
                    gsl_vector_get (&evec_i.vector, 1),
                    gsl_vector_get (&evec_i.vector, 2)
                    );

        // free stuff:
        gsl_eigen_symmv_free(w);
        gsl_vector_free(eval);
        gsl_matrix_free(evec);
    }

    void addPoint( const point_t & pi , float weight = 1.f )
    {
        vals[0] += weight*pi[0]*pi[0];
        vals[1] += weight*pi[0]*pi[1];
        vals[2] += weight*pi[0]*pi[2];
        vals[4] += weight*pi[1]*pi[1];
        vals[5] += weight*pi[1]*pi[2];
        vals[8] += weight*pi[2]*pi[2];
    }

    void removePoint( const point_t & pi , float weight = 1.f )
    {
        vals[0] -= weight*pi[0]*pi[0];
        vals[1] -= weight*pi[0]*pi[1];
        vals[2] -= weight*pi[0]*pi[2];
        vals[4] -= weight*pi[1]*pi[1];
        vals[5] -= weight*pi[1]*pi[2];
        vals[8] -= weight*pi[2]*pi[2];
    }

    void addNormal( const normal_t & pi )
    {
        vals[0] += pi[0]*pi[0];
        vals[1] += pi[0]*pi[1];
        vals[2] += pi[0]*pi[2];
        vals[4] += pi[1]*pi[1];
        vals[5] += pi[1]*pi[2];
        vals[8] += pi[2]*pi[2];
    }

    void removeNormal( const normal_t & pi )
    {
        vals[0] -= pi[0]*pi[0];
        vals[1] -= pi[0]*pi[1];
        vals[2] -= pi[0]*pi[2];
        vals[4] -= pi[1]*pi[1];
        vals[5] -= pi[1]*pi[2];
        vals[8] -= pi[2]*pi[2];
    }

    float operator()( int i , int j ) const
    {
        return vals[i + 3*j];
    }

    void set( float xx , float yy , float zz , float xy , float yz , float xz )
    {
        vals[0] = xx;
        vals[1] = xy;
        vals[2] = xz;
        vals[4] = yy;
        vals[5] = yz;
        vals[8] = zz;
    }

    void add( float xx , float yy , float zz , float xy , float yz , float xz )
    {
        vals[0] += xx;
        vals[1] += xy;
        vals[2] += xz;
        vals[4] += yy;
        vals[5] += yz;
        vals[8] += zz;
    }

    void operator = ( const basic_symetric_matrix3f< point_t , normal_t > & c )
    {
        vals[0] = c.vals[0];
        vals[1] = c.vals[1];
        vals[2] = c.vals[2];
        vals[4] = c.vals[4];
        vals[5] = c.vals[5];
        vals[8] = c.vals[8];
    }

    void operator += ( const basic_symetric_matrix3f< point_t , normal_t > & c )
    {
        vals[0] += c.vals[0];
        vals[1] += c.vals[1];
        vals[2] += c.vals[2];
        vals[4] += c.vals[4];
        vals[5] += c.vals[5];
        vals[8] += c.vals[8];
    }
};





namespace PCATools{

template< typename point_t , typename normal_t >
class PCASolver3f
{
protected:
    point_t Origine_;
    normal_t RepX_ , RepY_ , RepZ_;
    float lambdaX_ , lambdaY_ , lambdaZ_;

    basic_symetric_matrix3f< point_t , normal_t > cov;

    bool _order_and_normalize;

public:
    PCASolver3f()
    {
        _order_and_normalize = true;
        Origine_ = point_t(0,0,0);
        cov.clear();
    }
    PCASolver3f( bool _o_a_n )
    {
        _order_and_normalize = _o_a_n;
        Origine_ = point_t(0,0,0);
        cov.clear();
    }
    ~PCASolver3f(){}

    void clear()
    {
        Origine_ = point_t(0,0,0);
        cov.clear();
    }

    void order_and_normalize( bool _o_a_n )
    {
        _order_and_normalize = _o_a_n;
    }

    void order_basis()
    {
        float _l;
        normal_t _r;

        // WE ORDER THE BASIS THAT WAY :  lambdaX_ < lambdaY_ < lambdaZ_

        if( lambdaY_ < lambdaX_ )
        {
            _l = lambdaY_;
            _r = RepY_;

            lambdaY_ = lambdaX_;
            RepY_ = -RepX_;

            lambdaX_ = _l;
            RepX_ = _r;
        }

        if( lambdaZ_ < lambdaY_ )
        {
            _l = lambdaY_;
            _r = RepY_;

            lambdaY_ = lambdaZ_;
            RepY_ = -RepZ_;

            lambdaZ_ = _l;
            RepZ_ = _r;
        }
    }

    void diagonalize()
    {
        cov.diagonalize(lambdaX_ , lambdaY_ , lambdaZ_ , RepX_ , RepY_ , RepZ_);
    }

    void compute()
    {
        this->diagonalize();

        if( this->_order_and_normalize )
        {
            RepX_.normalize();
            RepZ_ = cross(RepX_ , RepY_);
            RepZ_.normalize();
            RepY_ = cross(RepZ_ , RepX_);
        }
    }

    void addPoint( const point_t & p , float weight = 1.f )
    {
        const point_t & pi = p - Origine_;

        cov.addPoint( pi , weight );
    }

    void removePoint( const point_t & p , float weight = 1.f )
    {
        const point_t & pi = p - Origine_;

        cov.removePoint( pi , weight );
    }

    void setOrigine( const point_t & o ){ Origine_ = o; }

    normal_t RepX(){ return RepX_; }
    normal_t RepY(){ return RepY_; }
    normal_t RepZ(){ return RepZ_; }
    float lambdaX(){ return lambdaX_; }
    float lambdaY(){ return lambdaY_; }
    float lambdaZ(){ return lambdaZ_; }
};





template< typename point_t , typename normal_t >
class DirectionSolver : public PCASolver3f< point_t , normal_t >
{
private:
    bool _computed;
    float _pertinence;
    float _planar;
    float _rich;
    normal_t _direction;

public:
    DirectionSolver()
    {
        this->_computed = false;
        this->_order_and_normalize = false;
        this->Origine_ = point_t(0,0,0);

        this->cov.clear();
    }
    DirectionSolver( bool _o_a_n )
    {
        this->_computed = false;
        this->_order_and_normalize = false;
        this->Origine_ = point_t(0,0,0);

        this->cov.clear();
    }
    ~DirectionSolver(){}

    void clear()
    {
        this->_computed = false;

        this->cov.clear();
    }

    void order_and_normalize( bool _o_a_n ){  }

    void setOrigine( const point_t & o ){  }


    normal_t XXYYZZ()
    {
        return normal_t(this->cov(0,0) , this->cov(1,1) , this->cov(2,2));
    }
    normal_t XYYZXZ()
    {
        return normal_t(this->cov(0,1) , this->cov(1,2) , this->cov(0,2));
    }

    DirectionSolver( const normal_t & xxyyzz , const normal_t & xyyzxz )
    {
        this->cov.set( xxyyzz[0] , xxyyzz[1] , xxyyzz[2] , xyyzxz[0] , xyyzxz[1] , xyyzxz[2] );
    }


    void addPoint( const point_t & pi )
    {
        assert( false && "Don't use POINTS with DirectionSolver !!  USE NORMALS !!" );
    }
    void addNormal( const normal_t & ni )
    {
        this->_computed = false;

        this->cov.addNormal( ni );
    }

    void removePoint( const point_t & pi )
    {
        assert( false && "Don't use POINTS with DirectionSolver !!  USE NORMALS !!" );
    }
    void removeNormal( const normal_t & ni )
    {
        this->_computed = false;

        this->cov.removeNormal( ni );
    }

    void compute()
    {
        if( this->_computed )
            return;

        this->diagonalize();


        this->order_basis();
        // Now we have lambdaX < lambdaY < lambdaZ :


        this->_planar = ( this->lambdaZ_ ) / ( this->lambdaX_ + this->lambdaY_ + this->lambdaZ_ );
        this->_planar = ( 3.f * (this->_planar) - 1.f ) / 2.f;


        //            this->_rich = ( this->lambdaX_ + this->lambdaY_ ) / ( this->lambdaX_ + this->lambdaY_ + this->lambdaZ_ );
        //            this->_rich = 3.f * _rich / 2.f ;
        this->_rich = 1.f  -  this->_planar;


        this->_pertinence = this->lambdaZ_ / ( this->lambdaZ_ + this->lambdaY_ );
        this->_pertinence = 2 * ( this->_pertinence ) - 1.f;


        this->_direction = this->RepZ_;


        this->_computed = true;
    }



    float planar()
    {
        // planar() belongs to [ 0 ; 1 ] :
        //      _planar == 0.f => not planar at all, inside the DirectionSolvver were put normals covering the unit sphere in a homogeneous way
        //      _planar == 1.f => completely planar, there is one normal inside the DirectionSolver
        this->compute();
        return this->_planar;
    }



    float pertinence()
    {
        this->compute();
        return this->_pertinence;
    }


    float rich()
    {
        this->compute();
        return this->_rich;
    }


    float poor()
    {
        this->compute();
        return 1.f  -  this->_rich;
    }



    normal_t direction()
    {
        // direction is the axis with the higher eigenvalue : that should be the normal that is the most colinear to every normal inside the DirectionSolver
        this->compute();
        return this->_direction;
    }



    void operator << ( const normal_t & p )
    {
        this->addNormal( p );
    }
    void operator += ( const normal_t & p )
    {
        this->addNormal( p );
    }

#if normal_t != point_t
    void operator << ( const point_t & p )
    {
        this->addPoint( p );
    }
    void operator += ( const point_t & p )
    {
        this->addPoint( p );
    }
#endif




    void setSolver( const DirectionSolver & S )
    {
        this->_computed = false;

        this->cov = S.cov;
    }

    void addSolver( const DirectionSolver & S )
    {
        this->_computed = false;

        this->cov += (S.cov);
    }

    void addComponents( const normal_t & xxyyzz , const normal_t & xyyzxz )
    {
        this->_computed = false;

        this->cov.add( xxyyzz[0] , xxyyzz[1] , xxyyzz[2] , xyyzxz[0] , xyyzxz[1] , xyyzxz[2] );
    }


    void operator<< ( const DirectionSolver & S )
    {
        this->addSolver( S );
    }
    void operator+= ( const DirectionSolver & S )
    {
        this->addSolver( S );
    }
    void operator= ( const DirectionSolver & S )
    {
        this->setSolver( S );
    }
};














template< typename point_t , typename normal_t  >
class AdaptivePCASolver : public PCASolver3f< point_t , normal_t >
{
private:
    bool _computed;
    float _sum_wi_pi_pi[3][3];
    point_t _sum_wi_pi;
    float _sum_wi;
    point_t _mean_point;

    float _planar;
    normal_t _plan_normal;

public:

    AdaptivePCASolver() :
        _computed(false),
        _sum_wi_pi(0.f,0.f,0.f),
        _sum_wi(0.f)
    {
        _sum_wi_pi_pi[0][0] = 0.f;
        _sum_wi_pi_pi[0][1] = 0.f;
        _sum_wi_pi_pi[0][2] = 0.f;
        //      _sum_wi_pi_pi[1][0] = 0.f;
        _sum_wi_pi_pi[1][1] = 0.f;
        _sum_wi_pi_pi[1][2] = 0.f;
        //      _sum_wi_pi_pi[2][0] = 0.f;
        //      _sum_wi_pi_pi[2][1] = 0.f;
        _sum_wi_pi_pi[2][2] = 0.f;
    }



    void addPoint( const point_t  & _p , float _weight = 1.f )
    {
        _computed = false;
        _sum_wi += _weight;
        _sum_wi_pi += _weight*_p;

        _sum_wi_pi_pi[0][0] += _weight * _p[0] * _p[0];
        _sum_wi_pi_pi[0][1] += _weight * _p[0] * _p[1];
        _sum_wi_pi_pi[0][2] += _weight * _p[0] * _p[2];
        _sum_wi_pi_pi[1][1] += _weight * _p[1] * _p[1];
        _sum_wi_pi_pi[1][2] += _weight * _p[1] * _p[2];
        _sum_wi_pi_pi[2][2] += _weight * _p[2] * _p[2];
    }





    void compute()
    {
        if( this->_computed )
            return;

        this->_mean_point = this->_sum_wi_pi / this->_sum_wi;

        this->cov.set
                (
                    _sum_wi_pi_pi[0][0] - _sum_wi * _mean_point[0] * _mean_point[0],
                _sum_wi_pi_pi[1][1] - _sum_wi * _mean_point[1] * _mean_point[1],
                _sum_wi_pi_pi[2][2] - _sum_wi * _mean_point[2] * _mean_point[2],
                _sum_wi_pi_pi[0][1] - _sum_wi * _mean_point[0] * _mean_point[1],
                _sum_wi_pi_pi[1][2] - _sum_wi * _mean_point[1] * _mean_point[2],
                _sum_wi_pi_pi[0][2] - _sum_wi * _mean_point[0] * _mean_point[2]
                );

        this->diagonalize();

        this->order_basis();

        this->_planar = this->lambdaX_ / ( this->lambdaX_ + this->lambdaY_ + this->lambdaZ_ );
        this->_planar = 1.f  -  3.f  *  (this->_planar);

        this->_plan_normal = this->RepX_;

        this->_computed = true;
    }




    point_t mean_point()
    {
        this->compute();
        return this->_mean_point;
    }



    float planar()
    {
        this->compute();
        return this->_planar;
    }



    void best_fitting_plan( float & _planar , point_t & _center , normal_t & _normal )
    {
        this->compute();

        _planar = this->_planar;
        _center = this->_mean_point;
        _normal = this->_plan_normal;
    }





    void addSolver( const AdaptivePCASolver & S )
    {
        this->_computed = false;

        this->_sum_wi_pi_pi[0][0] += S._sum_wi_pi_pi[0][0];
        this->_sum_wi_pi_pi[0][1] += S._sum_wi_pi_pi[0][1];
        this->_sum_wi_pi_pi[0][2] += S._sum_wi_pi_pi[0][2];
        //     this->_sum_wi_pi_pi[1][0] += S._sum_wi_pi_pi[1][0];
        this->_sum_wi_pi_pi[1][1] += S._sum_wi_pi_pi[1][1];
        this->_sum_wi_pi_pi[1][2] += S._sum_wi_pi_pi[1][2];
        //     this->_sum_wi_pi_pi[2][0] += S._sum_wi_pi_pi[2][0];
        //     this->_sum_wi_pi_pi[2][1] += S._sum_wi_pi_pi[2][1];
        this->_sum_wi_pi_pi[2][2] += S._sum_wi_pi_pi[2][2];

        this->_sum_wi += S._sum_wi;
        this->_sum_wi_pi += S._sum_wi_pi;
    }


    void setSolver( const AdaptivePCASolver & S )
    {
        this->_computed = false;

        this->_sum_wi_pi_pi[0][0] = S._sum_wi_pi_pi[0][0];
        this->_sum_wi_pi_pi[0][1] = S._sum_wi_pi_pi[0][1];
        this->_sum_wi_pi_pi[0][2] = S._sum_wi_pi_pi[0][2];
        //    this->_sum_wi_pi_pi[1][0] = S._sum_wi_pi_pi[1][0];
        this->_sum_wi_pi_pi[1][1] = S._sum_wi_pi_pi[1][1];
        this->_sum_wi_pi_pi[1][2] = S._sum_wi_pi_pi[1][2];
        //    this->_sum_wi_pi_pi[2][0] = S._sum_wi_pi_pi[2][0];
        //    this->_sum_wi_pi_pi[2][1] = S._sum_wi_pi_pi[2][1];
        this->_sum_wi_pi_pi[2][2] = S._sum_wi_pi_pi[2][2];

        this->_sum_wi = S._sum_wi;
        this->_sum_wi_pi = S._sum_wi_pi;
    }


    void operator << ( const AdaptivePCASolver & S )
    {
        this->addSolver( S );
    }
    void operator += ( const AdaptivePCASolver & S )
    {
        this->addSolver( S );
    }
    void operator = ( const AdaptivePCASolver & S )
    {
        this->setSolver( S );
    }


};



}












#endif // PCATOOLS_H
