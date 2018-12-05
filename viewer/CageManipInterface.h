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
#ifndef CAGEMANIPINTERFACE_H
#define CAGEMANIPINTERFACE_H

#include <vector>
#include <queue>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <cmath>

#include <omp.h>

// shared files from jmShared:
#include "geom/BasicIO.h"
#include "geom/PCATools.h"
#include "qt/Manipulator.h"

// MAIN FILES FOR CAGE COORDINATES:
#include "../coordinates/gc/gc.h"
#include "../coordinates/mvc/mvc.h"
#include "../coordinates/qmvc/qmvc.h"
#include "../coordinates/smvc/smvc.h"
#include "../coordinates/mec/mec.h"

#include <QTimer>


#define ALLOW_TRI_MVC
#define ALLOW_QMVC
#define ALLOW_QMVC_MEC // DO NOT COMPUTE QMVC_MEC IF YOU DO NOT COMPUTE QMVC !!!!
#define ALLOW_SMVC
#define ALLOW_GC



enum MeshModificationMode {VERTEX_UPDATE_INTERACTIVE , VERTEX_UPDATE_REALTIME};
enum MeshNormalsUpdateMode {NORMAL_UPDATE_NONE, NORMAL_UPDATE_INTERACTIVE , NORMAL_UPDATE_REALTIME};
enum MeshModificationMethod {
#ifdef ALLOW_TRI_MVC
    CAGECOORDS_MVC,
#endif
#ifdef ALLOW_QMVC
    CAGECOORDS_QMVC,
#endif
#ifdef ALLOW_QMVC_MEC
    CAGECOORDS_QMVC_MEC,
#endif
#ifdef ALLOW_SMVC
    CAGECOORDS_SMVC,
#endif
#ifdef ALLOW_GC
    CAGECOORDS_GREEN_Urago,
#endif
    CAGECOORDS_NUMBER_OF_COORDINATE_SYSTEMS
};











template< class point_t >
class CMInterface
{
    // MESH :
    std::vector< point_t > mesh_vertices;
    std::vector< std::vector< int > > mesh_triangles;

    // MODIFIED MESH :
    std::vector< point_t > mesh_modified_vertices;
    std::vector< point_t > mesh_modified_normals;

    // CAGE :
    std::vector< point_t > cage_vertices;
    std::vector< std::vector< int > > cage_triangles;
    std::vector< point_t > cage_triangle_normals;

    // CAGE VISU :
    std::vector< std::vector< int > > cage_visu_triangles;
    std::vector< std::vector< int > > cage_visu_quads;
    std::vector< std::pair< int,int > > cage_quads_cuts;

    // CAGE MANIP :
    std::vector< bool > cage_selected_vertices;

    // Deformation parameters :
    MeshModificationMode deformationMode;
    MeshModificationMethod deformationMethod;
    MeshNormalsUpdateMode normalsUpdateMode;





    // GREEN COORDS :
    std::vector< GreenCoordinates3D::GreenScalingFactor< point_t > > GREENscalingFactors;

    // GREEN COORDS Urago FORMULATION :
    std::vector< std::vector< double > > GREENphiCoordinates_Urago;
    std::vector< std::vector< double > > GREENpsiCoordinates_Urago;

    // MVC COORDS :
    std::vector< std::vector< double > > MVCCoordinates;

    // QMVC :
    std::vector< std::vector< double > > QMVCCoordinates;
    std::vector< std::vector< double > > QMVCCoordinates_MEC;
    std::vector< point_t > nanVertices;

    // SMVC :
    std::vector< std::vector< double > > SMVCCoordinates;





public:

    //  ACCESS the input :
    inline std::vector< point_t > & get_mesh_vertices() { return mesh_vertices ; }
    inline std::vector< point_t > const & get_mesh_vertices() const { return mesh_vertices ; }

    inline std::vector< point_t > & get_mesh_vertex_normals() { return mesh_modified_normals ; }
    inline std::vector< point_t > const & get_mesh_vertex_normals() const { return mesh_modified_normals ; }
    
    inline std::vector< std::vector< int > > & get_mesh_triangles() { return mesh_triangles ; }
    inline std::vector< std::vector< int > > const & get_mesh_triangles() const { return mesh_triangles ; }
    
    inline std::vector< point_t > & get_cage_vertices() { return cage_vertices ; }
    inline std::vector< point_t > const & get_cage_vertices() const { return cage_vertices ; }
    
    inline std::vector< std::vector< int > > & get_cage_visu_triangles() { return cage_visu_triangles ; }
    inline std::vector< std::vector< int > > const & get_cage_visu_triangles() const { return cage_visu_triangles ; }

    inline std::vector< std::vector< int > > & get_cage_visu_quads() { return cage_visu_quads ; }
    inline std::vector< std::vector< int > > const & get_cage_visu_quads() const { return cage_visu_quads ; }
    inline std::vector< std::vector< int > > & get_cage_quads_cuts() { return cage_quads_cuts ; }
    inline std::vector< std::vector< int > > const & get_cage_quads_cuts() const { return cage_quads_cuts ; }

    
    // ACCESS the output :
    inline std::vector< point_t > const & get_mesh_modified_vertices() const { return mesh_modified_vertices ; }
    
    inline std::vector< bool > & get_cage_selected_vertices () { return cage_selected_vertices; }
    inline const std::vector< bool > & get_cage_selected_vertices () const { return cage_selected_vertices; }





    CMInterface()
    {
        deformationMode = VERTEX_UPDATE_REALTIME;
        normalsUpdateMode = NORMAL_UPDATE_REALTIME;
    }
    void setMode( int m )
    {
        deformationMode = MeshModificationMode(m);
    }
    void setNormalUpdateMode( int m )
    {
        normalsUpdateMode = MeshNormalsUpdateMode(m);
    }
    void setMethod( int m )
    {
        deformationMethod = MeshModificationMethod(m);
    }


    void clear()
    {
        mesh_vertices.clear();
        mesh_triangles.clear();

        mesh_modified_vertices.clear();
        mesh_modified_normals.clear();

        cage_vertices.clear();
        cage_triangles.clear();
        cage_triangle_normals.clear();

        cage_selected_vertices.clear();

        cage_visu_triangles.clear();
        cage_visu_quads.clear();
        cage_quads_cuts.clear();

        GREENscalingFactors.clear();
        GREENphiCoordinates_Urago.clear();
        GREENpsiCoordinates_Urago.clear();

        MVCCoordinates.clear();

        QMVCCoordinates.clear();

        SMVCCoordinates.clear();
    }



    void reInitialize()
    {
        mesh_modified_normals.clear();

        mesh_modified_vertices.clear();
        for( unsigned int v = 0 ; v < mesh_vertices.size() ; ++v )
            mesh_modified_vertices.push_back( mesh_vertices[ v ] );

        cage_triangle_normals.clear();

        cage_selected_vertices.clear();
        cage_selected_vertices.resize( mesh_vertices.size() , false);

        cage_triangles.clear();

        GREENscalingFactors.clear();
        GREENphiCoordinates_Urago.clear();
        GREENpsiCoordinates_Urago.clear();

        MVCCoordinates.clear();

        QMVCCoordinates.clear();

        SMVCCoordinates.clear();


        // mesh normals :
        update_mesh_vertex_normals();

        // convert cage TriQuad (cage_visu_quads and cage_visu_triangles) to triangles (cage_triangles) :
        for( unsigned int t = 0 ; t < cage_visu_triangles.size() ; ++t )
        {
            int _v1 = cage_visu_triangles[t][0] , _v2 = cage_visu_triangles[t][1] , _v3 = cage_visu_triangles[t][2];
            std::vector< int > _v;

            _v.push_back( _v1 );
            _v.push_back( _v2 );
            _v.push_back( _v3 );

            cage_triangles.push_back( _v );

            GREENscalingFactors.push_back( GreenCoordinates3D::GreenScalingFactor< point_t >( cage_vertices[_v2] - cage_vertices[_v1] , cage_vertices[_v3] - cage_vertices[_v1] ) );
        }
        for( unsigned int q = 0 ; q < cage_visu_quads.size() ; ++q )
        {
            int _v1 = cage_visu_quads[q][0] , _v2 = cage_visu_quads[q][1] , _v3 = cage_visu_quads[q][2] , _v4 = cage_visu_quads[q][3];
            std::vector< int > _v;

            cage_quads_cuts.push_back( std::pair<int,int>(_v1,_v3) );

            _v.push_back( _v1 );
            _v.push_back( _v2 );
            _v.push_back( _v3 );

            cage_triangles.push_back( _v );

            GREENscalingFactors.push_back( GreenCoordinates3D::GreenScalingFactor< point_t >( cage_vertices[_v2] - cage_vertices[_v1] , cage_vertices[_v3] - cage_vertices[_v1] ) );

            _v.clear();
            _v.push_back( _v1 );
            _v.push_back( _v3 );
            _v.push_back( _v4 );

            cage_triangles.push_back( _v );

            GREENscalingFactors.push_back( GreenCoordinates3D::GreenScalingFactor< point_t >( cage_vertices[_v3] - cage_vertices[_v1] , cage_vertices[_v4] - cage_vertices[_v1] ) );
        }

        // cage normals :
        update_cage_triangle_normals();

        // compute cage coordinates :
        computeCoordinates();
    }





    void update_mesh_vertex_normals()
    {
        mesh_modified_normals.clear();
        mesh_modified_normals.resize( mesh_modified_vertices.size() , point_t(0,0,0) );

#pragma omp parallel for
        for( unsigned int t = 0 ; t < mesh_triangles.size(); ++t )
        {
            point_t const & tri_normal = ( mesh_modified_vertices[ mesh_triangles[t][1] ] - mesh_modified_vertices[ mesh_triangles[t][0] ] ) %
                    ( mesh_modified_vertices[ mesh_triangles[t][2] ] - mesh_modified_vertices[ mesh_triangles[t][0] ] );

            mesh_modified_normals[ mesh_triangles[t][0] ] += tri_normal;
            mesh_modified_normals[ mesh_triangles[t][1] ] += tri_normal;
            mesh_modified_normals[ mesh_triangles[t][2] ] += tri_normal;
        }

#pragma omp parallel for
        for( unsigned int v = 0 ; v < mesh_vertices.size() ; ++v )
        {
            mesh_modified_normals[ v ].normalize();
        }
    }



    void update_cage_triangle_normals()
    {
        cage_triangle_normals.resize( cage_triangles.size() );
#pragma omp parallel for
        for( unsigned int t = 0 ; t < cage_triangles.size(); ++t )
        {
            cage_triangle_normals[ t ] = ( cage_vertices[ cage_triangles[t][1] ] - cage_vertices[ cage_triangles[t][0] ] ) %
                    ( cage_vertices[ cage_triangles[t][2] ] - cage_vertices[ cage_triangles[t][0] ] );

            cage_triangle_normals[ t ].normalize();
        }
    }


    void update_cage_triangle_GREENscalingFactors()
    {
#pragma omp parallel for
        for( unsigned int t = 0 ; t < cage_triangles.size() ; ++t )
        {
            GREENscalingFactors[ t ].computeScalingFactor( cage_vertices[ cage_triangles[t][1] ] - cage_vertices[ cage_triangles[t][0] ] ,
                    cage_vertices[ cage_triangles[t][2] ] - cage_vertices[ cage_triangles[t][0] ] );

            //   GREENscalingFactors[ t ].computeScalingFactorBasedOnAreas( cage_vertices[ cage_triangles[t][1] ] - cage_vertices[ cage_triangles[t][0] ] ,
            //           cage_vertices[ cage_triangles[t][2] ] - cage_vertices[ cage_triangles[t][0] ] );
        }
    }




    void saveDeformedCage( std::string const & filename )
    {
        std::ofstream myfile;
        myfile.open(filename.c_str());
        if (!myfile.is_open())
        {
            std::cout << filename << " cannot be opened" << std::endl;
            return;
        }

        myfile << "OFF" << std::endl;
        myfile << (cage_vertices.size()) << " " << (cage_visu_quads.size() + cage_visu_triangles.size()) << " 0" << std::endl;

        for( unsigned int v = 0 ; v < cage_vertices.size() ; ++v )
        {
            myfile << (cage_vertices[v]) << std::endl;
        }

        for( unsigned int t = 0 ; t < cage_visu_triangles.size() ; ++t )
        {
            myfile << "3 " << (cage_visu_triangles[t][0]) << " " << (cage_visu_triangles[t][1]) << " " << (cage_visu_triangles[t][2]) << std::endl;
        }
        for( unsigned int t = 0 ; t < cage_visu_quads.size() ; ++t )
        {
            myfile << "4 " << (cage_visu_quads[t][0]) << " " << (cage_visu_quads[t][1]) << " " << (cage_visu_quads[t][2]) << " " << (cage_visu_quads[t][3]) << std::endl;
        }

        myfile.close();
    }


    void saveDeformedModel( std::string const & filename )
    {
        std::ofstream myfile;
        myfile.open(filename.c_str());
        if (!myfile.is_open())
        {
            std::cout << filename << " cannot be opened" << std::endl;
            return;
        }

        myfile << "OFF" << std::endl;
        myfile << (mesh_modified_vertices.size()) << " " << (mesh_triangles.size()) << " 0" << std::endl;

        for( unsigned int v = 0 ; v < mesh_modified_vertices.size() ; ++v )
        {
            myfile << (mesh_modified_vertices[v]) << std::endl;
        }

        for( unsigned int t = 0 ; t < mesh_triangles.size() ; ++t )
        {
            myfile << "3 " << (mesh_triangles[t][0]) << " " << (mesh_triangles[t][1]) << " " << (mesh_triangles[t][2]) << std::endl;
        }

        myfile.close();
    }



    void open_OFF_mesh( std::string const & filename  ,  point_t & bb  ,  point_t & BB )
    {
        OFFIO::open(filename , mesh_vertices , mesh_triangles , true);
        mesh_modified_vertices = mesh_vertices;
        bb = point_t( FLT_MAX , FLT_MAX , FLT_MAX );
        BB = point_t( -FLT_MAX , -FLT_MAX , -FLT_MAX );
        for( unsigned int v = 0 ; v < mesh_vertices.size() ; ++v )
        {
            point_t const & vpos = mesh_vertices[v];
            bb = point_t::min( bb , vpos );
            BB = point_t::max( BB , vpos );
        }
        update_mesh_vertex_normals();
    }
    void open_OBJ_mesh( std::string const & filename  ,  point_t & bb  ,  point_t & BB )
    {
        OBJIO::open(filename , mesh_vertices , mesh_triangles , true);

        mesh_modified_vertices = mesh_vertices;
        bb = point_t( FLT_MAX , FLT_MAX , FLT_MAX );
        BB = point_t( -FLT_MAX , -FLT_MAX , -FLT_MAX );
        for( unsigned int v = 0 ; v < mesh_vertices.size() ; ++v )
        {
            point_t const & vpos = mesh_vertices[v];
            bb = point_t::min( bb , vpos );
            BB = point_t::max( BB , vpos );
        }
        update_mesh_vertex_normals();

    }

    void initialize_cage_structure(std::vector< std::vector< unsigned int > > const & cfaces) {
        cage_selected_vertices.clear();
        cage_selected_vertices.resize(cage_vertices.size() , false);

        cage_triangles.clear();
        cage_visu_triangles.clear();
        cage_visu_quads.clear();
        GREENscalingFactors.clear();
        for( unsigned int f = 0 ; f < cfaces.size() ; ++f )
        {
            unsigned int n_vertices_on_face = cfaces[f].size();
            if( n_vertices_on_face == 3 )
            {
                int _v1 = cfaces[f][0] , _v2 = cfaces[f][1] , _v3 = cfaces[f][2];
                std::vector< int > _v;
                _v.push_back( _v1 );
                _v.push_back( _v2 );
                _v.push_back( _v3 );

                cage_triangles.push_back( _v );

                GREENscalingFactors.push_back( GreenCoordinates3D::GreenScalingFactor< point_t >( cage_vertices[_v2] - cage_vertices[_v1] , cage_vertices[_v3] - cage_vertices[_v1] ) );

                // VISU :
                cage_visu_triangles.push_back( _v );
            }
            else if( n_vertices_on_face == 4 )
            {
                int _v1 = cfaces[f][0] , _v2 = cfaces[f][1] , _v3 = cfaces[f][2] , _v4 = cfaces[f][3];
                std::vector< int > _v;
                _v.push_back( _v1 );
                _v.push_back( _v2 );
                _v.push_back( _v3 );

                cage_triangles.push_back( _v );

                GREENscalingFactors.push_back( GreenCoordinates3D::GreenScalingFactor< point_t >( cage_vertices[_v2] - cage_vertices[_v1] , cage_vertices[_v3] - cage_vertices[_v1] ) );

                _v.clear();
                _v.push_back( _v1 );
                _v.push_back( _v3 );
                _v.push_back( _v4 );

                cage_triangles.push_back( _v );

                GREENscalingFactors.push_back( GreenCoordinates3D::GreenScalingFactor< point_t >( cage_vertices[_v3] - cage_vertices[_v1] , cage_vertices[_v4] - cage_vertices[_v1] ) );

                // VISU :
                _v.clear();
                _v.push_back( _v1 );
                _v.push_back( _v2 );
                _v.push_back( _v3 );
                _v.push_back( _v4 );
                cage_visu_quads.push_back( _v );
            }
            else
            {
                std::cout << "We handle ONLY *.off files with 3 or 4 vertices per face" << std::endl;
                exit(1);
            }
        }
    }

    void open_OFF_cage( std::string const & filename )
    {
        std::vector< std::vector< unsigned int > > cfaces;
        OFFIO::open(filename , cage_vertices , cfaces , false);
        initialize_cage_structure(cfaces);
        update_cage_triangle_normals();
        computeCoordinates();
    }
    void open_OBJ_cage( std::string const & filename )
    {
        std::vector< std::vector< unsigned int > > cfaces;
        OBJIO::open(filename , cage_vertices , cfaces , false);
        initialize_cage_structure(cfaces);
        update_cage_triangle_normals();
        computeCoordinates();
    }

    void setCageDeformedVertices( std::vector< point_t > const & newCageVerts ) {
        if( cage_vertices.size() == newCageVerts.size() ) {
            cage_vertices = newCageVerts;

            {
                update_cage_triangle_normals();

                if( deformationMode == VERTEX_UPDATE_REALTIME )
                {
                    updateModel();
                    if( normalsUpdateMode == NORMAL_UPDATE_REALTIME ) {
                        update_mesh_vertex_normals();
                    }
                }
            }
        }
    }


    // This function gives you the index for the cage face you clicked on :  TODO
    int cageClickedFace( int x , int y , bool & isAQuad )
    {
        return -1;
    }


    //  You  don't  really  need  that :
    void drawCage( unsigned int nSubdivisionsForQuads = 5 )
    {
        glBegin( GL_TRIANGLES );

        // Draw triangles :
        for (unsigned int t = 0 ; t < cage_visu_triangles.size() ; ++t )
        {
            for( unsigned int v = 0 ; v < 3 ; ++v )
                glVertex3f( cage_vertices[ cage_visu_triangles[t][v] ][0],
                        cage_vertices[ cage_visu_triangles[t][v] ][1],
                        cage_vertices[ cage_visu_triangles[t][v] ][2]);
        }
        glEnd();
        // Draw quads :
        glBegin( GL_QUADS );
        for (unsigned int q = 0 ; q < cage_visu_quads.size() ; ++q )
        {
            if(nSubdivisionsForQuads > 0) {
                for( unsigned int k = 0 ; k < nSubdivisionsForQuads ; ++k ) {
                    for( unsigned int kv = 0 ; kv < nSubdivisionsForQuads ; ++kv ) {
                        double u = (double)(k) / (double)(nSubdivisionsForQuads) ;
                        double v = (double)(kv) / (double)(nSubdivisionsForQuads) ;
                        point_t vert = (1-u)*(1-v) * cage_vertices[ cage_visu_quads[q][0] ]
                                + u*(1-v) * cage_vertices[ cage_visu_quads[q][1] ]
                                + u*v * cage_vertices[ cage_visu_quads[q][2] ]
                                + (1-u)*v * cage_vertices[ cage_visu_quads[q][3] ];
                        glVertex3f( vert[0] , vert[1] , vert[2] );
                        u = (double)(k+1) / (double)(nSubdivisionsForQuads) ;
                        v = (double)(kv) / (double)(nSubdivisionsForQuads) ;
                        vert = (1-u)*(1-v) * cage_vertices[ cage_visu_quads[q][0] ]
                                + u*(1-v) * cage_vertices[ cage_visu_quads[q][1] ]
                                + u*v * cage_vertices[ cage_visu_quads[q][2] ]
                                + (1-u)*v * cage_vertices[ cage_visu_quads[q][3] ];
                        glVertex3f( vert[0] , vert[1] , vert[2] );
                        u = (double)(k+1) / (double)(nSubdivisionsForQuads) ;
                        v = (double)(kv+1) / (double)(nSubdivisionsForQuads) ;
                        vert = (1-u)*(1-v) * cage_vertices[ cage_visu_quads[q][0] ]
                                + u*(1-v) * cage_vertices[ cage_visu_quads[q][1] ]
                                + u*v * cage_vertices[ cage_visu_quads[q][2] ]
                                + (1-u)*v * cage_vertices[ cage_visu_quads[q][3] ];
                        glVertex3f( vert[0] , vert[1] , vert[2] );
                        u = (double)(k) / (double)(nSubdivisionsForQuads) ;
                        v = (double)(kv+1) / (double)(nSubdivisionsForQuads) ;
                        vert = (1-u)*(1-v) * cage_vertices[ cage_visu_quads[q][0] ]
                                + u*(1-v) * cage_vertices[ cage_visu_quads[q][1] ]
                                + u*v * cage_vertices[ cage_visu_quads[q][2] ]
                                + (1-u)*v * cage_vertices[ cage_visu_quads[q][3] ];
                        glVertex3f( vert[0] , vert[1] , vert[2] );
                    }
                }
            }
            else
                for( unsigned int v = 0 ; v < 4 ; ++v )
                    glVertex3f( cage_vertices[ cage_visu_quads[q][v] ][0],
                            cage_vertices[ cage_visu_quads[q][v] ][1],
                            cage_vertices[ cage_visu_quads[q][v] ][2]);
        }
        glEnd();
    }
    void drawcage_quads_cuts(  )
    {
        glBegin( GL_LINES );

        // Draw cage_quads_cuts :
        for (unsigned int c = 0 ; c < cage_quads_cuts.size() ; ++c )
        {
            glVertex3f( cage_vertices[ cage_quads_cuts[c].first ][0],
                    cage_vertices[ cage_quads_cuts[c].first ][1],
                    cage_vertices[ cage_quads_cuts[c].first ][2]);
            glVertex3f( cage_vertices[ cage_quads_cuts[c].second ][0],
                    cage_vertices[ cage_quads_cuts[c].second ][1],
                    cage_vertices[ cage_quads_cuts[c].second ][2]);
        }
        glEnd();
    }



    //  You  don't  really  need  that :
    void drawBaseModel()
    {
        glEnable( GL_LIGHTING );
        glPolygonMode( GL_FRONT_AND_BACK , GL_FILL );
        glBegin( GL_TRIANGLES );
        for( unsigned int t = 0 ; t < mesh_triangles.size() ; ++t )
        {
            for( unsigned int v = 0 ; v < 3 ; ++v )
            {
                glNormal3f( mesh_modified_normals[ mesh_triangles[t][v] ][0],
                        mesh_modified_normals[ mesh_triangles[t][v] ][1],
                        mesh_modified_normals[ mesh_triangles[t][v] ][2]);
                glVertex( mesh_vertices[ mesh_triangles[t][v] ] );
            }
        }
        glEnd();
    }


    //  You  don't  really  need  that :
    void drawModifiedModel(  )
    {
        glBegin( GL_TRIANGLES );
        for( unsigned int t = 0 ; t < mesh_triangles.size() ; ++t )
        {
            for( unsigned int v = 0 ; v < 3 ; ++v )
            {
                glNormal3f( mesh_modified_normals[ mesh_triangles[t][v] ][0],
                        mesh_modified_normals[ mesh_triangles[t][v] ][1],
                        mesh_modified_normals[ mesh_triangles[t][v] ][2]);
                glVertex3f( mesh_modified_vertices[ mesh_triangles[t][v] ][0],
                        mesh_modified_vertices[ mesh_triangles[t][v] ][1],
                        mesh_modified_vertices[ mesh_triangles[t][v] ][2]);
            }
        }
        glEnd();
    }
    template< class ScalarTextureHandler_t >
    void drawModifiedModelWithTexture(ScalarTextureHandler_t * scalarTexture , int coordinate_id_to_show = -1 )
    {
        glBegin( GL_TRIANGLES );
        for( unsigned int t = 0 ; t < mesh_triangles.size() ; ++t )
        {
            for( unsigned int v = 0 ; v < 3 ; ++v )
            {
                double value = getCoordinateValue(mesh_triangles[t][v],coordinate_id_to_show);
                scalarTexture->setVal(value);
                glNormal3f( mesh_modified_normals[ mesh_triangles[t][v] ][0],
                        mesh_modified_normals[ mesh_triangles[t][v] ][1],
                        mesh_modified_normals[ mesh_triangles[t][v] ][2]);
                glVertex3f( mesh_modified_vertices[ mesh_triangles[t][v] ][0],
                        mesh_modified_vertices[ mesh_triangles[t][v] ][1],
                        mesh_modified_vertices[ mesh_triangles[t][v] ][2]);
            }
        }
        glEnd();
    }


    double getCoordinateValue( unsigned int v , int coordinate_id_to_show ) const {
#ifdef ALLOW_GC
        if( deformationMethod == CAGECOORDS_GREEN_Urago ) {
            if( coordinate_id_to_show < 0 ) return 0.0;
            if( coordinate_id_to_show >= (int)(cage_vertices.size() + cage_triangles.size()) ) return 0.0;
            if( coordinate_id_to_show < (int)(cage_vertices.size()) ) return GREENphiCoordinates_Urago[v][coordinate_id_to_show];
            return GREENpsiCoordinates_Urago[v][coordinate_id_to_show - cage_vertices.size()];
        }
#endif

#ifdef ALLOW_TRI_MVC
        if( deformationMethod == CAGECOORDS_MVC ) {
            if( coordinate_id_to_show < 0 ) return 0.0;
            if( coordinate_id_to_show >= (int)cage_vertices.size() ) return 0.0;
            return MVCCoordinates[v][coordinate_id_to_show];
        }
#endif
#ifdef ALLOW_SMVC
        if( deformationMethod == CAGECOORDS_SMVC ) {
            if( coordinate_id_to_show < 0 ) return 0.0;
            if( coordinate_id_to_show >= (int)cage_vertices.size() ) return 0.0;
            return SMVCCoordinates[v][coordinate_id_to_show];
        }
#endif
#ifdef ALLOW_QMVC
        if( deformationMethod == CAGECOORDS_QMVC ) {
            if( coordinate_id_to_show < 0 ) return 0.0;
            if( coordinate_id_to_show >= (int)cage_vertices.size() ) return 0.0;
            return QMVCCoordinates[v][coordinate_id_to_show];
        }
#endif
#ifdef ALLOW_QMVC_MEC
        if( deformationMethod == CAGECOORDS_QMVC_MEC ) {
            if( coordinate_id_to_show < 0 ) return 0.0;
            if( coordinate_id_to_show >= (int)cage_vertices.size() ) return 0.0;
            return QMVCCoordinates_MEC[v][coordinate_id_to_show];
        }
#endif

        return 0.0;
    }


    //  You  don't  really  need  that :
    void drawCageSelectedVertices()
    {
        glBegin( GL_POINTS );
        for( unsigned int v = 0 ; v < cage_vertices.size() ; ++v )
        {
            if( cage_selected_vertices[ v ] )
                glVertex3f( cage_vertices[ v ][0],
                        cage_vertices[ v ][1],
                        cage_vertices[ v ][2]);
        }
        glEnd();
    }

    void drawNANVertices()
    {
        glBegin( GL_POINTS );
        for( unsigned int v = 0 ; v < nanVertices.size() ; ++v )
        {
                glVertex3f( nanVertices[ v ][0],
                        nanVertices[ v ][1],
                        nanVertices[ v ][2]);
        }
        glEnd();
    }










    void computeGreenCoordinates_Urago()
    {
        QTime timer;
        timer.start();


        GREENphiCoordinates_Urago.clear();
        GREENpsiCoordinates_Urago.clear();
        GREENphiCoordinates_Urago.resize( mesh_vertices.size() );
        GREENpsiCoordinates_Urago.resize( mesh_vertices.size() );

        for (unsigned int i =0; i < mesh_vertices.size (); i++) {
            GREENphiCoordinates_Urago[i].clear();
            GREENpsiCoordinates_Urago[i].clear();
            GREENphiCoordinates_Urago[i].resize( cage_vertices.size() , 0.f );
            GREENpsiCoordinates_Urago[i].resize( cage_triangles.size() , 0.f );
        }
        timer.restart();
#pragma omp parallel for
        for( unsigned int p_idx = 0 ; p_idx < mesh_vertices.size() ; ++p_idx )
        {
            GreenCoordinates3D::URAGO::computeCoordinates(
                        mesh_vertices[ p_idx ] ,
                        cage_triangles ,
                        cage_vertices ,
                        GREENphiCoordinates_Urago[ p_idx ],
                        GREENpsiCoordinates_Urago[ p_idx ]);
        }

        std::cout << "Computed Green coordinates in " << timer.elapsed() << " ms (CPU version)" << std::endl;
    }






    void computeMVCCoordinates()
    {
        QTime timer;
        timer.start();
        MVCCoordinates.clear();
        MVCCoordinates.resize( mesh_vertices.size() );

#pragma omp parallel for 
        for( unsigned int p_idx = 0 ; p_idx < mesh_vertices.size() ; ++p_idx )
        {
            MVC3D::computeCoordinates(
                        mesh_vertices[ p_idx ] ,
                        cage_triangles ,
                        cage_vertices ,
                        cage_triangle_normals ,
                        MVCCoordinates[ p_idx ]);
        }
        std::cout << "Computed mean value coordinates in " << timer.elapsed() << " ms (CPU version)" << std::endl;
    }



    void computeSMVCCoordinates()
    {
        QTime timer;
        timer.start();
        SMVCCoordinates.clear();
        SMVCCoordinates.resize( mesh_vertices.size() );

        std::vector< std::vector< int > > cage_faces;
        for(unsigned int f = 0 ; f < cage_visu_triangles.size() ; ++f)
            cage_faces.push_back( cage_visu_triangles[f] );
        for(unsigned int f = 0 ; f < cage_visu_quads.size() ; ++f)
            cage_faces.push_back( cage_visu_quads[f] );

#pragma omp parallel for
        for( unsigned int p_idx = 0 ; p_idx < mesh_vertices.size() ; ++p_idx )
        {
            SMVC3D::computeCoordinates(
                        mesh_vertices[ p_idx ] ,
                        cage_faces ,
                        cage_vertices ,
                        SMVCCoordinates[ p_idx ]);
        }
        std::cout << "Computed Spherical Mean Value Coordinates in " << timer.elapsed() << " ms (CPU version)" << std::endl;
    }



    void computeQMVCCoordinates()
    {
        QMVCCoordinates.clear();
        QMVCCoordinates.resize( mesh_vertices.size() );

        QTime timer;
        timer.start();
#pragma omp parallel for
        for( unsigned int p_idx = 0 ; p_idx < mesh_vertices.size() ; ++p_idx )
        {
            QMVC3D::computeCoordinates(
                        mesh_vertices[ p_idx ] ,
                        cage_visu_triangles ,
                        cage_visu_quads ,
                        cage_vertices ,
                        QMVCCoordinates[ p_idx ]);
        }
        std::cout << "Computed Quad mean value coordinates in " << timer.elapsed() << " ms (CPU version)" << std::endl;

        // YOU NEED TO COMMENT THE OMP PRAGMA IF YOU WANT TO HAVE A PROPER ESTIMATION OF THE COMPUTATION TIME NEEDED PER COORDINATE !!!!!!!!
        //        double elapsed = (double)(timer.elapsed());
        //        std::cout << "Computed "<< mesh_vertices.size() << " x " << cage_vertices.size() << " = " << mesh_vertices.size() * cage_vertices.size() << " coordinates in " << elapsed << " ms : " << (elapsed)/(mesh_vertices.size()*cage_vertices.size())
        //                  << "ms per coordinate and " << (elapsed)/(mesh_vertices.size()) << " ms per vertex for a cage made of " << cage_vertices.size() << " vertices, " <<
        //                     cage_visu_triangles.size() << " triangles and " << cage_visu_quads.size() << " quads" << std::endl;
    }


    void computeCoordinates()
    {
#ifdef ALLOW_TRI_MVC
        computeMVCCoordinates();
#endif
#ifdef ALLOW_QMVC
        computeQMVCCoordinates();
        #ifdef ALLOW_QMVC_MEC
                computeQMVCMECCoordinates();
        #endif
#endif
#ifdef ALLOW_SMVC
        computeSMVCCoordinates();
#endif
#ifdef ALLOW_GC
        computeGreenCoordinates_Urago();
#endif
    }






    void computeQMVCMECCoordinates() {
        unsigned int nV = QMVCCoordinates.size();
        QMVCCoordinates_MEC.resize(nV);
//        double coordinateShift = 0.05;
        QTime timer;
        timer.start();
#pragma omp parallel for
        for(unsigned int v = 0 ; v < nV ; ++v) {
            std::vector< double > clampCoords( QMVCCoordinates[v].size() );
            for( unsigned int i = 0 ; i < QMVCCoordinates[v].size() ; ++i ) {
//                clampCoords[i] = std::max<double>( 0.0 , QMVCCoordinates[v][i] - coordinateShift );
                clampCoords[i] = std::max<double>( 0.0 , QMVCCoordinates[v][i] );
            }
            bool isNAN = MEC::computeCoordinates<mat33d,point3d,double,double>( mesh_vertices[v] , cage_vertices , clampCoords , QMVCCoordinates_MEC[v] );
            if(isNAN)
                nanVertices.push_back(mesh_vertices[v]);
        }
        std::cout << "Converted Quad mean value coordinates with MEC in " << timer.elapsed() << " ms (CPU version)" << std::endl;
    }









    // You need to have the correct QGLContext activated for that !!!!!!!!
    // This function select the cage vertices drawn inside the QRect "zone"
    void select( QRectF const & zone )
    {
        float modelview[16];
        glGetFloatv(GL_MODELVIEW_MATRIX , modelview);
        float projection[16];
        glGetFloatv(GL_PROJECTION_MATRIX , projection);

        for( unsigned int v = 0 ; v < cage_vertices.size() ; ++v )
        {
            point_t const & p = cage_vertices[ v ];

            float x = modelview[0] * p[0] + modelview[4] * p[1] + modelview[8] * p[2] + modelview[12];
            float y = modelview[1] * p[0] + modelview[5] * p[1] + modelview[9] * p[2] + modelview[13];
            float z = modelview[2] * p[0] + modelview[6] * p[1] + modelview[10] * p[2] + modelview[14];
            float w = modelview[3] * p[0] + modelview[7] * p[1] + modelview[11] * p[2] + modelview[15];
            x /= w;
            y /= w;
            z /= w;
            w = 1.f;

            float xx = projection[0] * x + projection[4] * y + projection[8] * z + projection[12] * w;
            float yy = projection[1] * x + projection[5] * y + projection[9] * z + projection[13] * w;
            float ww = projection[3] * x + projection[7] * y + projection[11] * z + projection[15] * w;
            xx /= ww;
            yy /= ww;

            xx = ( xx + 1.f ) / 2.f;
            yy = 1.f - ( yy + 1.f ) / 2.f;

            if( zone.contains( xx , yy ) )
                cage_selected_vertices[ v ] = true;
        }
    }

    // You need to have the correct QGLContext activated for that !!!!!!!!
    // This function unselect the cage vertices drawn inside the QRect "zone"
    void unselect( QRectF const & zone )
    {
        float modelview[16];
        glGetFloatv(GL_MODELVIEW_MATRIX , modelview);
        float projection[16];
        glGetFloatv(GL_PROJECTION_MATRIX , projection);

        for( unsigned int v = 0 ; v < cage_vertices.size() ; ++v )
        {
            point_t const & p = cage_vertices[ v ];

            float x = modelview[0] * p[0] + modelview[4] * p[1] + modelview[8] * p[2] + modelview[12];
            float y = modelview[1] * p[0] + modelview[5] * p[1] + modelview[9] * p[2] + modelview[13];
            float z = modelview[2] * p[0] + modelview[6] * p[1] + modelview[10] * p[2] + modelview[14];
            float w = modelview[3] * p[0] + modelview[7] * p[1] + modelview[11] * p[2] + modelview[15];
            x /= w;
            y /= w;
            z /= w;
            w = 1.f;

            float xx = projection[0] * x + projection[4] * y + projection[8] * z + projection[12] * w;
            float yy = projection[1] * x + projection[5] * y + projection[9] * z + projection[13] * w;
            float ww = projection[3] * x + projection[7] * y + projection[11] * z + projection[15] * w;
            xx /= ww;
            yy /= ww;

            xx = ( xx + 1.f ) / 2.f;
            yy = 1.f - ( yy + 1.f ) / 2.f;

            if( zone.contains( xx , yy ) )
                cage_selected_vertices[ v ] = false;
        }
    }





    // Once you're OK with the selection you made , call this function to compute the manipulator , and it will activate it :
    void computeManipulatorForSelection( SimpleManipulator * manipulator )
    {
        manipulator->resetScales();
        manipulator->clear();

        int nb=0;
        qglviewer::Vec oo( 0.f , 0.f , 0.f );
        for( unsigned int v = 0 ; v < cage_vertices.size() ; ++v )
        {
            if( cage_selected_vertices[v] )
            {
                point_t const & p = cage_vertices[v];
                oo += qglviewer::Vec( p[0] , p[1] , p[2] );
                ++nb;
            }
        }
        oo /= nb;

        PCATools::PCASolver3f< qglviewer::Vec , qglviewer::Vec > solver;
        solver.setOrigine( oo );

        for( unsigned int v = 0 ; v < cage_vertices.size() ; ++v )
        {
            if( cage_selected_vertices[v] )
            {
                point_t const & p = cage_vertices[v];
                solver.addPoint( qglviewer::Vec( p[0] , p[1] , p[2] ) );
            }
        }

        solver.compute();

        manipulator->setOrigine( oo );
        manipulator->setRepX( solver.RepX() );
        manipulator->setRepY( solver.RepY() );
        manipulator->setRepZ( solver.RepZ() );

        for( unsigned int v = 0 ; v < cage_vertices.size() ; ++v )
        {
            if( cage_selected_vertices[v] )
            {
                /// ! Changed ! ///
                //cage_selected_vertices[v] = false;
                point_t const & p = cage_vertices[v];
                manipulator->addPoint( v , qglviewer::Vec( p[0] , p[1] , p[2] ) );
            }
        }
        manipulator->activate();
    }




    // When you move your manipulator , it sends you a SIGNAL.
    // When it happens, call that function with the manipulator as the parameter, it will update everything :
    void cageChanged( SimpleManipulator * manipulator )
    {
        unsigned int n_points = manipulator->n_points();
        qglviewer::Vec p;
        int idx;
        for( unsigned int i = 0 ; i < n_points ; ++i )
        {
            manipulator->getTransformedPoint( i , idx , p );

            cage_vertices[ idx ] = point_t( p[0] , p[1] , p[2] );
        }
        update_cage_triangle_normals();

        if( deformationMode == VERTEX_UPDATE_REALTIME )
        {
            updateModel();
            if( normalsUpdateMode == NORMAL_UPDATE_REALTIME ) {
                update_mesh_vertex_normals();
            }
        }
    }


    // When you release the mouse after moving the manipulator, it sends you a SIGNAL.
    // When it happens, call that function with the manipulator as the parameter, it will update everything :
    void manipulatorReleased()
    {
        if( deformationMode == VERTEX_UPDATE_INTERACTIVE )
        {
            updateModel();
        }
        if( normalsUpdateMode == NORMAL_UPDATE_INTERACTIVE ) {
            update_mesh_vertex_normals();
        }
    }



    // You should never have to call that function on your own, but why not :
    void updateModel()
    {
#ifdef ALLOW_TRI_MVC
        if( deformationMethod == CAGECOORDS_MVC )
        {
#pragma omp parallel for
            for( unsigned int v = 0 ; v < mesh_modified_vertices.size() ; ++v )
            {
                point_t pos(0,0,0);
                for( unsigned int vc = 0 ; vc < MVCCoordinates[ v ].size() ; ++vc )
                    pos += MVCCoordinates[ v ][ vc ] * cage_vertices[ vc ];

                mesh_modified_vertices[ v ] = pos;
            }
            return;
        }
#endif
#ifdef ALLOW_QMVC
        if( deformationMethod == CAGECOORDS_QMVC )
        {
#pragma omp parallel for
            for( unsigned int v = 0 ; v < mesh_modified_vertices.size() ; ++v )
            {
                point_t pos(0,0,0);
                for( unsigned int vc = 0 ; vc < QMVCCoordinates[ v ].size() ; ++vc )
                    pos += QMVCCoordinates[ v ][ vc ] * cage_vertices[ vc ];

                mesh_modified_vertices[ v ] = pos;
            }
            return;
        }
#endif
#ifdef ALLOW_QMVC_MEC
        if( deformationMethod == CAGECOORDS_QMVC_MEC )
        {
#pragma omp parallel for
            for( unsigned int v = 0 ; v < mesh_modified_vertices.size() ; ++v )
            {
                point_t pos(0,0,0);
                for( unsigned int vc = 0 ; vc < QMVCCoordinates_MEC[ v ].size() ; ++vc )
                    pos += QMVCCoordinates_MEC[ v ][ vc ] * cage_vertices[ vc ];

                mesh_modified_vertices[ v ] = pos;
            }
            return;
        }
#endif
#ifdef ALLOW_SMVC
        if( deformationMethod == CAGECOORDS_SMVC )
        {
#pragma omp parallel for
            for( unsigned int v = 0 ; v < mesh_modified_vertices.size() ; ++v )
            {
                point_t pos(0,0,0);
                for( unsigned int vc = 0 ; vc < SMVCCoordinates[ v ].size() ; ++vc )
                    pos += SMVCCoordinates[ v ][ vc ] * cage_vertices[ vc ];

                mesh_modified_vertices[ v ] = pos;
            }
            return;
        }
#endif
#ifdef ALLOW_GC
        if( deformationMethod == CAGECOORDS_GREEN_Urago )
        {
            update_cage_triangle_GREENscalingFactors();

#pragma omp parallel for
            for( unsigned int v = 0 ; v < mesh_modified_vertices.size() ; ++v )
            {
                point_t pos(0,0,0);
                for( unsigned int vc = 0 ; vc < cage_vertices.size() ; ++vc )
                    pos += GREENphiCoordinates_Urago[ v ][ vc ] * cage_vertices[vc];
                for( unsigned int tc = 0 ; tc < cage_triangles.size() ; ++tc ){
                    pos += GREENpsiCoordinates_Urago[ v ][ tc ] * GREENscalingFactors[ tc ].scalingFactor() * cage_triangle_normals[tc];
                }

                mesh_modified_vertices[ v ] = pos;
            }
            return;
        }
#endif
    }
};







#endif // CAGEMANIPINTERFACE_H
