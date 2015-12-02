#ifndef _ORBIFOLD_HXX
#define _ORBIFOLD_HXX

#include "MyMesh.hxx"
#include "ShortestPathDijkstra.hxx"

#define OF_INTERNAL 0
#define OF_BOUNDARY 1
#define OF_CONESINGULARITY_PI  2
#define OF_CONESINGULARITY_PI2 3

class Orbifold {
public:

  Orbifold() {};
  Orbifold( MyMesh& mymesh ) {
    init( mymesh );
  };
  ~Orbifold() {
    for ( int i = 0; i < path_.size(); ++i ) path_[i].clear();
    path_.clear();
  };
  void init( MyMesh& mymesh ) {
    setMyMesh( mymesh );
    vertex_type_.resize( mymesh.n_vertices() );
  };
  void setMyMesh( MyMesh& mymesh ) { mymesh_ = &mymesh; };
  MyMesh& mymesh() { return *mymesh_; };

  // void addCenter( MyMesh::VertexHandle& a, MyMesh::VertexHandle& b ) {
  //   center_.insert( make_pair( a, b ) );
  // };
  // MyMesh::VertexHandle center(MyMesh::VertexHandle& a) { return center_[a]; };

  std::vector<std::vector<MyMesh::VertexHandle> >& path() { return path_; };
  std::vector<unsigned int>& vertex_type() { return vertex_type_; };

  void setCSVertices( std::vector<MyMesh::VertexHandle>& cs_vertices ) {
    for ( int i = 0; i < cs_vertices.size(); ++i )
      {
        // cout << i << endl;
        cs_vertex_.push_back( cs_vertices[i] );
      }
  };
  std::vector<MyMesh::VertexHandle>& cs_vertices() { return cs_vertex_; };

  void setCSVertex( MyMesh::VertexHandle& vh ) {
    cs_vertex_.push_back( vh );
  };

  bool calcBoundaries() {
    if ( cs_vertex_.size() != 3 )
      {
        cerr << "the number of cone singularites is not three. " << endl;
        return false;
      }

    // four paths
    path_.resize(4);

    ShortestPathDijkstra sp(mymesh());

    // calc path no.1
    MyMesh::VertexHandle sv = cs_vertex_[0];
    MyMesh::VertexHandle ev = cs_vertex_[1];
    sp.apply( sv, ev, path_[0] );

    // copy path 0 to 3 in the reverse order
    for ( int i = path_[0].size() - 1; i >= 0; --i )
      {
        path_[3].push_back( path_[0][i] );
        int id = path_[0][i].idx();
        if ( i == 0 )
          vertex_type_[id] = OF_CONESINGULARITY_PI;
        else if ( i == path_[0].size() - 1 )
          vertex_type_[id] = OF_CONESINGULARITY_PI2;
        else
          vertex_type_[id] = OF_BOUNDARY;
      }

    // calc path no.2
    sv = cs_vertex_[1];
    ev = cs_vertex_[2];
    sp.apply( sv, ev, path_[1] );
 
    // copy path 1 to 2 in the reverse order
    for ( int i = path_[1].size() - 1; i >= 0; --i )
      {
        path_[2].push_back( path_[1][i] );
        int id = path_[1][i].idx();
        if ( i == 0 )
          vertex_type_[id] = OF_CONESINGULARITY_PI2;
        else if ( i == path_[1].size() - 1 )
          vertex_type_[id] = OF_CONESINGULARITY_PI;
        else
          vertex_type_[id] = OF_BOUNDARY;
      }

#if 0    
    for ( int i = 0; i < path_.size(); ++i )
      {
        cout << "path " << i+1 << endl;
        for ( int j = 0; j < path_[i].size(); ++j )
          {
            int id = path_[i][j].idx();
            cout << "\t v " << id << " type " << vertex_type_[id] << endl;
          }
      }
#endif
    return true;
  };


private:

  MyMesh* mymesh_;

  // // center OF_CONESINGULARITY_PI2 vertex for boundary vertices (in eq.9b)
  // std::map<MyMesh::VertexHandle,MyMesh::VertexHandle> center_;

  // first and last vertices are cone singularities
  std::vector<std::vector<MyMesh::VertexHandle> > path_;

  // 0, 2: OF_CONESINGULARITY_PI, 1: OF_CONESINGULARITY_PI2
  std::vector<MyMesh::VertexHandle> cs_vertex_;

  // vertex type (INTERNAL, BOUNDARY, OF_CONESINGULARITY_PI or OF_CONESINGULARITY_PI2)
  std::vector<unsigned int> vertex_type_;
};

#endif // _ORBIFOLD_HXX
