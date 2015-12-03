#ifndef _TOMESHR_HXX
#define _TOMESHR_HXX

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
};

#endif // _TOMESHR_HXX
