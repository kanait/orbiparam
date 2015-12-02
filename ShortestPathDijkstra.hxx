#ifndef _SHORTESTPATHDIJKSTRA_HXX
#define _SHORTESTPATHDIJKSTRA_HXX

#include <iostream>
#include <limits>
using namespace std;

#include "MyMesh.hxx"
#include "VertexPQ.hxx"

class ShortestPathDijkstra {

public:
  ShortestPathDijkstra() {};
  ShortestPathDijkstra( MyMesh& mymesh ) {
    setMyMesh( mymesh );

    init();
  };
  ~ShortestPathDijkstra() {
    mymesh().remove_property(elen);
  };

  void init( MyMesh& mymesh ) {
    setMyMesh( mymesh );
    init();
  };

  void init() {
    mymesh().add_property(elen);
    computeEdgeLength();
  };

  void apply_init() {
    verr_.resize( mymesh().n_vertices() );
    prev_vt_.resize( mymesh().n_vertices() );
    for ( int i = 0; i < mymesh().n_vertices(); ++i )
      {
        verr_[i] = std::numeric_limits<double>::max();
      }
  };

  void apply_finish() {
    verr_.clear();
    prev_vt_.clear();
  }

  void setMyMesh( MyMesh& mymesh ) { mymesh_ = &mymesh; };
  MyMesh& mymesh() { return *mymesh_; };

  void computeEdgeLength() {
    MyMesh& mesh = mymesh();
    MyMesh::EdgeIter e_it, e_end(mesh.edges_end());
    for (e_it=mesh.edges_begin(); e_it!=e_end; ++e_it)
      {
        MyMesh::HalfedgeHandle h = mesh.halfedge_handle( *e_it, 0 );
        MyMesh::VertexHandle v0 = mesh.to_vertex_handle( h );
        MyMesh::VertexHandle v1 = mesh.from_vertex_handle( h );
        MyMesh::Point p0 = mesh.point( v0 );
        MyMesh::Point p1 = mesh.point( v1 );
        Elen& el = mesh.property( elen, *e_it );
        el.setLen( (p1 - p0).length() );
      }
  };

  void apply( MyMesh::VertexHandle& sv, MyMesh::VertexHandle& ev,
              std::vector<MyMesh::VertexHandle>& path ) {

    apply_init();

    MyMesh& mesh = mymesh();

    MyMesh::VertexIHalfedgeIter vh_it;
    for ( vh_it = mesh.vih_iter( sv ); vh_it.is_valid(); ++vh_it )
        {
          MyMesh::VertexHandle ov = mesh.from_vertex_handle( *vh_it );
          Elen& el = mesh.property( elen, mesh.edge_handle( *vh_it ) );
          // push pair(el, ov)
          double len = el.len();
          verr_[ov.idx()] = len;
          prev_vt_[ov.idx()] = sv;
          pq_.push( len, ov );
        }

    //cout << "pq size " << pq_.size() << endl;

    while ( !(pq_.empty()) )
      {
        //cout << "err " << pq_.top_err() << " v " << pq_.top_vh().idx() << endl;
        double min_err = pq_.top_err();
        MyMesh::VertexHandle min_vt = pq_.top_vh();
        pq_.pop();

        MyMesh::VertexIHalfedgeIter nvh_it;
        for ( nvh_it = mesh.vih_iter( min_vt ); nvh_it.is_valid(); ++nvh_it )
          {
            MyMesh::VertexHandle nvt = mesh.from_vertex_handle( *nvh_it );
            Elen& nel = mesh.property( elen, mesh.edge_handle( *nvh_it ) );
            //cout << "\t " << verr_[nvt.idx()] << " " << min_err + nel.len() << endl;
            if ( verr_[nvt.idx()] > min_err + nel.len() )
              {
                double len = min_err + nel.len();
                verr_[nvt.idx()] = len;
                prev_vt_[nvt.idx()] = min_vt;
                pq_.push( len, nvt );
              }
          }
      }

    // create path
    std::vector<MyMesh::VertexHandle> path0;
    MyMesh::VertexHandle vt = ev;
    while ( vt != sv )
      {
        path0.push_back( vt );
        vt = prev_vt_[ vt.idx() ];
      }
    path0.push_back( vt );

    // copy path in reverse order
    for ( int i = path0.size()-1; i >= 0; --i )
      path.push_back(path0[i]);

    apply_finish();

  };

private:

  MyMesh* mymesh_;
  VertexPQ pq_;
  std::vector<double> verr_;
  std::vector<MyMesh::VertexHandle> prev_vt_;

};

#endif // _SHORTESTPATHDIJKSTRA_HXX
