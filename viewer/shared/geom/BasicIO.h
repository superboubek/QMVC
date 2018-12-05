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
#ifndef BASICIO_H
#define BASICIO_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <QString>
#include <QStringList>
#include <QFile>
#include <QDataStream>
#include <cfloat>
#include <cmath>
#include <cassert>



namespace StupidPointDontUseThatEverEverEver{
struct Point3f
{
    typedef float type_t;

    Point3f() :x(0),y(0),z(0) {}

    Point3f(float xI, float yI, float zI)
        :x(xI)
        ,y(yI)
        ,z(zI)
    {
        pos[0] = xI;
        pos[1] = yI;
        pos[2] = zI;
    }

    float operator[]( unsigned int i )const
    {
        return pos[i];
    }
    float & operator[]( unsigned int i )
    {
        return pos[i];
    }

    union
    {
        struct
        {
            float x;
            float y;
            float z;
        };
        float pos[3];
    };
};
}




namespace SCALARIO{

template< class T > bool open( const std::string & filename , std::vector< T > & function )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    std::string magic_s;

    myfile >> magic_s;

    if( magic_s != "SCALAR"   &&   magic_s != "SSF" )
    {
        std::cout << magic_s << " != SCALAR | SSF :   We handle ONLY *.sca and *.ssf files." << std::endl;
        myfile.close();
        return false;
    }

    unsigned int ns;
    myfile >> ns;

    function.resize(ns);
    for( unsigned int s = 0 ; s < ns ; ++s )
    {
        T x;
        myfile >> x;
        function[s] = x;
    }
    myfile.close();
    return true;
}

template< class T > bool save( const std::string & filename ,const std::vector< T > & function )
{
    std::ofstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    myfile << "SCALAR" << std::endl;
    myfile << function.size() << std::endl;

    unsigned int ns = function.size();
    for( unsigned int s = 0 ; s < ns ; ++s )
    {
        myfile << function[s] << std::endl;
    }
    myfile.close();
    return true;
}

}


namespace SOFFIO{

template< class point_t > bool open( const std::string & filename , std::vector< point_t > & vertices )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    std::string magic_s;

    myfile >> magic_s;

    if( magic_s != "SOFF" )
    {
        std::cout << magic_s << " != SOFF :   We handle ONLY *.soff files." << std::endl;
        myfile.close();
        return false;
    }

    unsigned int nv , ne;
    myfile >> nv >> ne;

    vertices.resize(nv);
    for( unsigned int v = 0 ; v < nv ; ++v )
    {
        typename point_t::type_t x , y , z;
        myfile >> x >> y >> z;
        vertices[v] = point_t(x,y,z);
    }
    myfile.close();
    return true;
}

template< class point_t > bool open( const std::string & filename , std::vector< point_t > & vertices , std::vector< int > & edges )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    std::string magic_s;

    myfile >> magic_s;

    if( magic_s != "SOFF" )
    {
        std::cout << magic_s << " != SOFF :   We handle ONLY *.soff files." << std::endl;
        myfile.close();
        return false;
    }

    unsigned int nv , ne;
    myfile >> nv >> ne;

    vertices.resize(nv);
    for( unsigned int v = 0 ; v < nv ; ++v )
    {
        typename point_t::type_t x , y , z;
        myfile >> x >> y >> z;
        vertices[v] = point_t(x,y,z);
    }
    edges.resize(2*ne);
    for( unsigned int e = 0 ; e < 2*ne ; ++e )
    {
        int ei;
        myfile >> ei;
        edges[e] = ei;
    }
    myfile.close();
    return true;
}

}


namespace OFFIO{

template< class point_t , class type_t > bool open( const std::string & filename ,
                                                    std::vector< point_t > & vertices ,
                                                    std::vector< std::vector< type_t > > & faces,
                                                    bool convertToTriangles = true,
                                                    bool randomize = false )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    std::string magic_s;

    myfile >> magic_s;

    if( magic_s != "OFF" )
    {
        std::cout << magic_s << " != OFF :   We handle ONLY *.off files." << std::endl;
        myfile.close();
        return false;
    }

    int n_vertices , n_faces , dummy_int;
    myfile >> n_vertices >> n_faces >> dummy_int;

    vertices.resize(n_vertices);

    for( int v = 0 ; v < n_vertices ; ++v )
    {
        typename point_t::type_t x , y , z;
        myfile >> x >> y >> z;
        if( std::isnan(x) )
            x = typename point_t::type_t(0.0);
        if( std::isnan(y) )
            y = typename point_t::type_t(0.0);
        if( std::isnan(z) )
            z = typename point_t::type_t(0.0);
        vertices[v] = point_t( x , y , z );
    }


    for( int f = 0 ; f < n_faces ; ++f )
    {
        int n_vertices_on_face;
        myfile >> n_vertices_on_face;
        if( n_vertices_on_face == 3 )
        {
            type_t _v1 , _v2 , _v3;
            std::vector< type_t > _v;
            myfile >> _v1 >> _v2 >> _v3;
            _v.push_back( _v1 );
            _v.push_back( _v2 );
            _v.push_back( _v3 );
            faces.push_back( _v );
        }
        else if( n_vertices_on_face > 3 )
        {
            std::vector< type_t > vhandles;
            vhandles.resize(n_vertices_on_face);
            for( int i=0 ; i < n_vertices_on_face ; ++i )
                myfile >> vhandles[i];

            if( convertToTriangles )
            {
                unsigned int k=(randomize)?(rand()%vhandles.size()):0;
                for (unsigned int i=0;i<vhandles.size()-2;++i)
                {
                    std::vector< type_t > tri(3);
                    tri[0]=vhandles[(k+0)%vhandles.size()];
                    tri[1]=vhandles[(k+i+1)%vhandles.size()];
                    tri[2]=vhandles[(k+i+2)%vhandles.size()];
                    faces.push_back(tri);
                }
            }
            else
            {
                faces.push_back(vhandles);
            }
        }
        else
        {
            std::cout << "OFFIO::open error : Face number " << f << " has " << n_vertices_on_face << " vertices" << std::endl;
            myfile.close();
            return false;
        }
    }

    myfile.close();
    return true;
}


template< class point_t , class type_t > bool save( const std::string & filename , std::vector< point_t > & vertices , std::vector< std::vector< type_t > > & faces )
{
    std::ofstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    myfile << "OFF" << std::endl;

    unsigned int n_vertices = vertices.size() , n_faces = faces.size();
    myfile << n_vertices << " " << n_faces << " 0" << std::endl;

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        myfile << vertices[v][0] << " " << vertices[v][1] << " " << vertices[v][2] << std::endl;
    }
    for( unsigned int f = 0 ; f < n_faces ; ++f )
    {
        myfile << faces[f].size();
        for( unsigned int vof = 0 ; vof < faces[f].size() ; ++vof )
            myfile << " " << faces[f][vof];
        myfile << std::endl;
    }
    myfile.close();
    return true;
}
template< class point_t , class face_t > bool saveFromMeshFormat( const std::string & filename , const point_t * vertices , unsigned int n_vertices , const face_t * faces , unsigned int n_faces )
{
    std::ofstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    myfile << "OFF" << std::endl;

    myfile << n_vertices << " " << n_faces << " 0" << std::endl;

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        myfile << vertices[v][0] << " " << vertices[v][1] << " " << vertices[v][2] << std::endl;
    }
    for( unsigned int f = 0 ; f < n_faces ; ++f )
    {
        myfile << faces[f].size();
        for( unsigned int vof = 0 ; vof < faces[f].size() ; ++vof )
            myfile << " " << faces[f][vof];
        myfile << std::endl;
    }
    myfile.close();
    return true;
}


template< class point_t > bool open( const std::string & filename , std::vector< point_t > & vertices )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    std::string magic_s;

    myfile >> magic_s;

    if( magic_s != "OFF" )
    {
        std::cout << magic_s << " != OFF :   We handle ONLY *.off files." << std::endl;
        myfile.close();
        return false;
    }

    int n_vertices , n_faces , dummy_int;
    myfile >> n_vertices >> n_faces >> dummy_int;

    vertices.resize(n_vertices);

    point_t bb( FLT_MAX , FLT_MAX , FLT_MAX );
    point_t BB( -FLT_MAX , -FLT_MAX , -FLT_MAX );

    for( int v = 0 ; v < n_vertices ; ++v )
    {
        typename point_t::type_t x , y , z;
        myfile >> x >> y >> z;
        vertices[v] = point_t( x , y , z );
        bb = point_t::min( bb , point_t( x , y , z ) );
        BB = point_t::max( BB , point_t( x , y , z ) );
    }

    myfile.close();

    return true;
}




}





namespace PTSIO{
template< class point_t > bool open(
        const std::string & filename,
        std::vector<point_t> & vertices)
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    while( myfile.good() )
    {
        std::string line;
        getline (myfile,line);
        QString QTLine = QString::fromStdString( line );
        QStringList lineElements = QTLine.split(" ");

        if( lineElements.size() != 3 )
            break;

        vertices.push_back(point_t( lineElements[0].toDouble() , lineElements[1].toDouble() , lineElements[2].toDouble() ));

    }
    std::cout << "THE FILE WAS READ SUCCESSFULLY" << std::endl;
    myfile.close();
    return true;
}
}




namespace OBJIO{



template< class point_t , class int_type_t > bool open(
        const std::string & filename,
        std::vector<point_t> & vertices,
        std::vector< std::vector< int_type_t > > & faces,
        bool convertToTriangles = true,
        bool randomize = false,
        bool convertEdgesToDegenerateTriangles = true)
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    vertices.clear();
    faces.clear();

    while( myfile.good() )
    {
        std::string line;
        getline (myfile,line);
        QString QTLine = QString::fromStdString( line );
        QRegExp reg("\\s+");
        QStringList lineElements = QTLine.split(reg);

        if(  lineElements.size() > 0  )
        {
            QString elementType = lineElements[0];
            // vertex
            if ( elementType == QString("v") )
            {
                vertices.push_back(point_t( lineElements[1].toDouble() , lineElements[2].toDouble() , lineElements[3].toDouble() ));
            }
            // face
            else if ( elementType == QString("f") )
            {
                std::vector< int_type_t > vhandles;
                for( int i = 1 ; i < lineElements.size() ; ++i )
                {
                    QStringList FaceElements = lineElements[i].split("/", QString::SkipEmptyParts);
                    if( FaceElements.size() > 0 )
                        vhandles.push_back( (int_type_t)( (abs(FaceElements[0].toInt()) - 1) ) );
                }

                if (vhandles.size()>3)
                {
                    if( convertToTriangles )
                    {
                        //model is not triangulated, so let us do this on the fly...
                        //to have a more uniform mesh, we add randomization
                        unsigned int k=(randomize)?(rand()%vhandles.size()):0;
                        for (unsigned int i=0;i<vhandles.size()-2;++i)
                        {
                            std::vector< int_type_t > tri(3);
                            tri[0]= (int_type_t) vhandles[(k+0)%vhandles.size()];
                            tri[1]= (int_type_t) vhandles[(k+i+1)%vhandles.size()];
                            tri[2]= (int_type_t) vhandles[(k+i+2)%vhandles.size()];
                            faces.push_back(tri);
                        }
                    }
                    else
                        faces.push_back(vhandles);
                }
                else if (vhandles.size()==3)
                {
                    faces.push_back(vhandles);
                }
                else if (vhandles.size()==2)
                {
                    if( convertEdgesToDegenerateTriangles )
                    {
                        printf("Unexpected number of face vertices (2). Converting edge to degenerate triangle");
                        vhandles.push_back(vhandles[1]);
                        faces.push_back(vhandles);
                    }
                    else
                        printf("Unexpected number of face vertices (2). Ignoring face");
                }
            }
        }
    }
    myfile.close();
    return true;
}



template< class point_t , class int_type_t > bool open(
        const std::string & filename,
        std::vector<point_t> & vertices,
        std::vector< std::vector< int_type_t > > & faces,
        std::vector< std::pair< int_type_t , int_type_t > > & additional_edges,
        bool convertToTriangles = true,
        bool randomize = false)
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    vertices.clear();
    faces.clear();

    while( myfile.good() )
    {
        std::string line;
        getline (myfile,line);
        QString QTLine = QString::fromStdString( line );
        QRegExp reg("\\s+");
        QStringList lineElements = QTLine.split(reg);

        if(  lineElements.size() > 0  )
        {
            QString elementType = lineElements[0];
            // vertex
            if ( elementType == QString("v") )
            {
                vertices.push_back(point_t( lineElements[1].toDouble() , lineElements[2].toDouble() , lineElements[3].toDouble() ));
            }
            // face
            else if ( elementType == QString("f") )
            {
                std::vector< int_type_t > vhandles;
                for( int i = 1 ; i < lineElements.size() ; ++i )
                {
                    QStringList FaceElements = lineElements[i].split("/", QString::SkipEmptyParts);
                    if( FaceElements.size() > 0 )
                        vhandles.push_back( (int_type_t)( (abs(FaceElements[0].toInt()) - 1) ) );
                }

                if (vhandles.size()>3)
                {
                    if( convertToTriangles )
                    {
                        //model is not triangulated, so let us do this on the fly...
                        //to have a more uniform mesh, we add randomization
                        unsigned int k=(randomize)?(rand()%vhandles.size()):0;
                        for (unsigned int i=0;i<vhandles.size()-2;++i)
                        {
                            std::vector< int_type_t > tri(3);
                            tri[0]= (int_type_t) vhandles[(k+0)%vhandles.size()];
                            tri[1]= (int_type_t) vhandles[(k+i+1)%vhandles.size()];
                            tri[2]= (int_type_t) vhandles[(k+i+2)%vhandles.size()];
                            faces.push_back(tri);
                        }
                    }
                    else
                        faces.push_back(vhandles);
                }
                else if (vhandles.size()==3)
                {
                    faces.push_back(vhandles);
                }
                else if (vhandles.size()==2)
                {
                    if( true )
                    {
                        // printf("Unexpected number of face vertices (2). Converting edge to degenerate triangle");
                        additional_edges.push_back( std::pair<int_type_t,int_type_t>(vhandles[0],vhandles[1]) );
                    }
                    else
                        printf("Unexpected number of face vertices (2). Ignoring face");
                }
            }
            else if ( elementType == QString("l") )
            {
                std::vector< int_type_t > vhandles;
                for( int i = 1 ; i < lineElements.size() ; ++i )
                {
                    QStringList FaceElements = lineElements[i].split("/", QString::SkipEmptyParts);
                    if( FaceElements.size() > 0 )
                        vhandles.push_back( (int_type_t)( (abs(FaceElements[0].toInt()) - 1) ) );
                }
                if (vhandles.size()==2)
                {
                    if( true )
                    {
                        // printf("Unexpected number of face vertices (2). Converting edge to degenerate triangle");
                        additional_edges.push_back( std::pair<int_type_t,int_type_t>(vhandles[0],vhandles[1]) );
                    }
                    else
                        printf("Unexpected number of face vertices (2). Ignoring face");
                }
            }
        }
    }
    myfile.close();
    return true;
}






template< class point_t > bool open(
        const std::string & filename,
        std::vector<point_t> & vertices)
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    vertices.clear();

    while( myfile.good() )
    {
        std::string line;
        getline (myfile,line);
        QString QTLine = QString::fromStdString( line );
        QRegExp reg("\\s+");
        QStringList lineElements = QTLine.split(reg);
        //        QStringList lineElementsBeforeCleaning = QTLine.split(" ", QString::SkipEmptyParts);
        //        QStringList lineElements;
        //        for( int e = 0 ; e < lineElementsBeforeCleaning.size() ; ++e )
        //        {
        //            if( lineElementsBeforeCleaning[e] != "" && lineElementsBeforeCleaning[e] != " " )
        //                lineElements.push_back( lineElementsBeforeCleaning[e] );
        //        }
        if(  lineElements.size() > 0  )
        {
            QString elementType = lineElements[0];
            // vertex
            if ( elementType == QString("v") )
            {
                vertices.push_back(point_t( lineElements[1].toDouble() , lineElements[2].toDouble() , lineElements[3].toDouble() ));
            }
        }
    }
    myfile.close();
    return true;
}








template< class point_t , typename int_t > bool save(
        std::string const & filename,
        std::vector<point_t> & verticesP,
        std::vector< std::vector<int_t> > & facesP )
{
    std::ofstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    myfile << "#OBJ" << std::endl;
    myfile << (verticesP.size()) << " " << (facesP.size()) << " 0" << std::endl;

    for( unsigned int v = 0 ; v < verticesP.size() ; ++v )
    {
        myfile << "v " << (verticesP[v][0]) << " " << (verticesP[v][1]) << " " << (verticesP[v][2]) << std::endl;
    }

    for( unsigned int t = 0 ; t < facesP.size() ; ++t )
    {
        unsigned int nv = facesP[t].size();
        myfile << "f";
        for( unsigned int vof = 0 ; vof < nv ; ++vof )
            myfile << " " << (facesP[t][vof] + 1);
        myfile << std::endl;
    }

    myfile.close();
    return true;
}





template< class point_t , class face_t , class edge_t > bool saveFromMeshFormat(
        const std::string & filename ,
        const point_t * vertices ,
        unsigned int n_vertices ,
        const face_t * faces ,
        unsigned int n_faces ,
        const edge_t * edges ,
        unsigned int n_edges )
{
    std::ofstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    for( unsigned int v = 0 ; v < n_vertices ; ++v )
    {
        myfile << "v " << vertices[v][0] << " " << vertices[v][1] << " " << vertices[v][2] << std::endl;
    }
    for( unsigned int f = 0 ; f < n_faces ; ++f )
    {
        myfile << "f";
        for( unsigned int vof = 0 ; vof < faces[f].size() ; ++vof )
            myfile << " " << (faces[f][vof]+1);
        myfile << std::endl;
    }
    for( unsigned int e = 0 ; e < n_edges ; ++e )
    {
        myfile << "l " << (edges[e][0]+1) << " " << (edges[e][1]+1) << std::endl;
    }
    myfile.close();
    return true;
}




}



namespace MeshIO{
template< class point_t > bool open( const std::string & filename , std::vector< point_t > & vertices , std::vector< int > & tetras )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    std::string dummy_s;

    unsigned int meshVersion , dataDimension , nVertices , nTriangles , nTetras;

    myfile >> dummy_s >> meshVersion;
    myfile >> dummy_s >> dataDimension;

    if( dataDimension != 3 )
    {
        std::cout << "dataDimension != 3 , " << filename << " is not a tetrahedral mesh file !" << std::endl;
        myfile.close();
        return false;
    }

    myfile >> dummy_s >> nVertices;
    vertices.resize( nVertices );
    for( unsigned int v = 0 ; v < nVertices ; ++v )
    {
        double x , y , z , dummy_i;
        myfile >> x >> y >> z >> dummy_i;
        vertices[v] = point_t(x,y,z);
    }

    myfile >> dummy_s >> nTriangles;
    for( unsigned int v = 0 ; v < nTriangles ; ++v )
    {
        int x , y , z , dummy_i;
        myfile >> x >> y >> z >> dummy_i;
    }

    myfile >> dummy_s >> nTetras;
    tetras.resize(4 * nTetras );
    for( unsigned int t = 0 ; t < nTetras ; ++t )
    {
        int w , x , y , z , dummy_i;
        myfile >> w >> x >> y >> z >> dummy_i;
        tetras[ 4*t ] = w-1;
        tetras[ 4*t + 1 ] = x-1;
        tetras[ 4*t + 2 ] = y-1;
        tetras[ 4*t + 3 ] = z-1;
    }
    myfile.close();

    return true;
}



template< class point_t > bool openTrisAndTets( const std::string & filename , std::vector< point_t > & vertices , std::vector< int > & tris , std::vector< int > & tetras )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    std::string dummy_s;

    unsigned int meshVersion , dataDimension , nVertices , nTriangles , nTetras;

    myfile >> dummy_s >> meshVersion;
    myfile >> dummy_s >> dataDimension;

    if( dataDimension != 3 )
    {
        std::cout << "dataDimension != 3 , " << filename << " is not a tetrahedral mesh file !" << std::endl;
        myfile.close();
        return false;
    }

    myfile >> dummy_s >> nVertices;
    vertices.resize( nVertices );
    for( unsigned int v = 0 ; v < nVertices ; ++v )
    {
        double x , y , z , dummy_i;
        myfile >> x >> y >> z >> dummy_i;
        vertices[v] = point_t(x,y,z);
    }

    myfile >> dummy_s >> nTriangles;
    tris.resize(3 * nTriangles );
    for( unsigned int t = 0 ; t < nTriangles ; ++t )
    {
        int x , y , z , dummy_i;
        myfile >> x >> y >> z >> dummy_i;
        tris[ 3*t + 0 ] = x-1;
        tris[ 3*t + 1 ] = y-1;
        tris[ 3*t + 2 ] = z-1;
    }

    myfile >> dummy_s >> nTetras;
    tetras.resize(4 * nTetras );
    for( unsigned int t = 0 ; t < nTetras ; ++t )
    {
        int w , x , y , z , dummy_i;
        myfile >> w >> x >> y >> z >> dummy_i;
        tetras[ 4*t ] = w-1;
        tetras[ 4*t + 1 ] = x-1;
        tetras[ 4*t + 2 ] = y-1;
        tetras[ 4*t + 3 ] = z-1;
    }
    myfile.close();

    return true;
}








template< class point_t , class int_t > bool openTris( const std::string & filename , std::vector< point_t > & vertices , std::vector< std::vector< int_t > > & tris )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    std::string dummy_s;

    unsigned int meshVersion , dataDimension , nVertices , nTriangles;
    //    unsigned int nTetras;

    myfile >> dummy_s >> meshVersion;
    myfile >> dummy_s >> dataDimension;

    if( dataDimension != 3 )
    {
        std::cout << "dataDimension != 3 , " << filename << " is not a tetrahedral mesh file !" << std::endl;
        myfile.close();
        return false;
    }

    myfile >> dummy_s >> nVertices;
    std::cout << "allocate " << nVertices << " vertices" << std::endl;
    vertices.resize( nVertices );
    for( unsigned int v = 0 ; v < nVertices ; ++v )
    {
        double x , y , z , dummy_i;
        myfile >> x >> y >> z >> dummy_i;
        vertices[v] = point_t(x,y,z);
    }

    myfile >> dummy_s >> nTriangles;
    std::cout << "allocate " << nTriangles << " tris" << std::endl;
    tris.resize( nTriangles );
    for( unsigned int t = 0 ; t < nTriangles ; ++t )
    {
        tris[t].resize(3);
        int_t x , y , z , dummy_i;
        myfile >> x >> y >> z >> dummy_i;
        tris[ t ][ 0 ] = x-1;
        tris[ t ][ 1 ] = y-1;
        tris[ t ][ 2 ] = z-1;
    }

    myfile.close();

    return true;
}

}



namespace BasicConvertor{

static
inline
bool OBJ2OFF(std::string const & OBJfilename, std::string const & OFFfilename)
{
    std::vector< StupidPointDontUseThatEverEverEver::Point3f > verticesP;
    std::vector< std::vector<int> > facesP;
    if( OBJIO::open(OBJfilename,verticesP,facesP,true , true) )
        return OFFIO::save(OFFfilename,verticesP,facesP);
    return false;
}

static
inline
bool OBJ2OFF(QString  fileName)
{
    bool madeIt = true;
    //    for(int i = 0 ; i < fileNames.size() ; ++i ) {
    //        QString fileName = fileNames[i];  // got a file name
    fileName.truncate( fileName.size() - 4 );
    madeIt = madeIt && OBJ2OFF( (fileName + QString(".obj")).toStdString() , (fileName + QString(".off")).toStdString() );
    //    }
    return madeIt;
}

static
inline
bool OBJ2OFF(QStringList const & fileNames)
{
    bool madeIt = true;
    for(int i = 0 ; i < fileNames.size() ; ++i ) {
        QString fileName = fileNames[i];  // got a file name
        fileName.truncate( fileName.size() - 4 );
        madeIt = madeIt && OBJ2OFF( (fileName + QString(".obj")).toStdString() , (fileName + QString(".off")).toStdString() );
    }
    return madeIt;
}



static
inline
bool OFF2OBJ(std::string const & OFFfilename, std::string const & OBJfilename)
{
    std::vector< StupidPointDontUseThatEverEverEver::Point3f > verticesP;
    std::vector< std::vector<int> > facesP;
    if( OFFIO::open(OFFfilename,verticesP,facesP,false) )
        return OBJIO::save(OBJfilename,verticesP,facesP);
}

static
inline
bool OFF2OBJ(QString  fileName)
{
    bool madeIt = true;
    //    for(int i = 0 ; i < fileNames.size() ; ++i ) {
    //        QString fileName = fileNames[i];  // got a file name
    fileName.truncate( fileName.size() - 4 );
    madeIt = madeIt && OFF2OBJ( (fileName + QString(".off")).toStdString() , (fileName + QString(".obj")).toStdString() );
    //    }
    return madeIt;
}
static
inline
bool OFF2OBJ(QStringList const & fileNames)
{
    bool madeIt = true;
    for(int i = 0 ; i < fileNames.size() ; ++i ) {
        QString fileName = fileNames[i];  // got a file name
        fileName.truncate( fileName.size() - 4 );
        madeIt = madeIt && OFF2OBJ( (fileName + QString(".off")).toStdString() , (fileName + QString(".obj")).toStdString() );
    }
    return madeIt;
}


}


#endif // BASICIO_H
