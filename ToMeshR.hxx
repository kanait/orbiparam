////////////////////////////////////////////////////////////////////
//
// $Id: ToMeshR.hxx 2021/06/13 22:13:49 kanai Exp $
//
// Copyright (c) 2021 Takashi Kanai
// Released under the MIT license
//
////////////////////////////////////////////////////////////////////

#ifndef _TOMESHR_HXX
#define _TOMESHR_HXX 1

#include "MyMesh.hxx"
#include "MeshR.hxx"

class ToMeshR {

public:
  // triangles only
  void apply( MyMesh& mymesh, MeshR& meshr ) {

    int n_v = mymesh.n_vertices();
    int n_f = mymesh.n_faces();
    if ( n_v )  meshr.reservePoints( n_v );
    if ( n_f )  meshr.reserveIndices( n_f );
    if ( mymesh.has_vertex_texcoords2D() )
      meshr.reserveTexcoords( n_v, 2 );

    MyMesh::VertexIter v_it, v_end(mymesh.vertices_end());
    int v_id = 0;
    for (v_it=mymesh.vertices_begin(); v_it!=v_end; ++v_it)
      {
        MyMesh::VertexHandle vh = *v_it;
        MyMesh::Point p = mymesh.point(vh);
        meshr.setPoint( v_id, (float)p[0], (float)p[1], (float)p[2] );
        v_id += nXYZ;
      }

    MyMesh::FaceIter f_it, f_end(mymesh.faces_end());
    int f_id = 0;
    for (f_it=mymesh.faces_begin(); f_it!=f_end; ++f_it)
      {
        MyMesh::FaceVertexIter fv_it;
        for ( fv_it = mymesh.fv_iter( *f_it ); fv_it.is_valid(); ++fv_it )
          {
            MyMesh::VertexHandle fvh = *fv_it;
            meshr.setIndex( f_id, fvh.idx() ); ++f_id;
          }
      }

    if ( mymesh.has_vertex_texcoords2D() )
      {
        for (v_it=mymesh.vertices_begin(); v_it!=v_end; ++v_it)
          {
            MyMesh::VertexHandle vh = *v_it;
            MyMesh::TexCoord2D t = mymesh.texcoord2D( vh );
            // cout << "id " << vh.idx() << " " << t[0] << " " << t[1] << endl;
            meshr.setTexcoord( 2*vh.idx(), t[0], t[1] );
          }
      }
  };

  void apply_param( MyMesh& mymesh, MeshR& meshr ) {

    int n_v = mymesh.n_vertices();
    int n_f = mymesh.n_faces();
    if ( !n_v || !n_f ) return;
    if ( !(mymesh.has_vertex_texcoords2D()) ) return;

    meshr.reservePoints( n_v );
    meshr.reserveIndices( n_f );

    // create vertex normals
    if ( !(mymesh.has_vertex_normals()) )
      {
        mymesh.request_face_normals();
        mymesh.request_vertex_normals();
        mymesh.update_normals();
        mymesh.update_vertex_normals();
      }

    meshr.reserveNormals( n_v );

    if ( mymesh.has_vertex_colors() )
      meshr.reserveColors( n_v );

    MyMesh::VertexIter v_it, v_end(mymesh.vertices_end());
    int v_id = 0;
    int n_id = 0;
    int c_id = 0;
    for (v_it=mymesh.vertices_begin(); v_it!=v_end; ++v_it)
      {
        MyMesh::VertexHandle vh = *v_it;
        MyMesh::TexCoord2D t = mymesh.texcoord2D( vh );
        meshr.setPoint( v_id, (float)t[0], (float)t[1], 0.0f );
        v_id += nXYZ;
        MyMesh::Normal n = mymesh.normal( vh );
        // cout << "n_id " << n_id << " " << n << endl;
        meshr.setNormal( n_id, (float)n[0], (float)n[1], (float)n[2] );
        n_id += nXYZ;
        if ( mymesh.has_vertex_colors() )
          {
            MyMesh::Color c = mymesh.color( vh );
            meshr.setColor( c_id, (unsigned char) c[0], (unsigned char) c[1], (unsigned char) c[2] );
            // cout << "c_id " << c_id << " c " << (int) c[0] << " " << (int) c[1] << " " << (int) c[2] << endl;
            c_id += nXYZ;
          }
      }

    MyMesh::FaceIter f_it, f_end(mymesh.faces_end());
    int f_id = 0;
    for (f_it=mymesh.faces_begin(); f_it!=f_end; ++f_it)
      {
        MyMesh::FaceVertexIter fv_it;
        for ( fv_it = mymesh.fv_iter( *f_it ); fv_it.is_valid(); ++fv_it )
          {
            MyMesh::VertexHandle fvh = *fv_it;
            meshr.setIndex( f_id, fvh.idx() ); ++f_id;
          }
      }

  };
};

#endif // _TOMESHR_HXX
