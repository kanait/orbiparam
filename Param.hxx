#ifndef _PARAM_HXX
#define _PARAM_HXX

#include "MyMesh.hxx"

#define EIGEN_NO_DEBUG

#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>
typedef Eigen::Triplet<double> T;

#define RECTANGLE 4

#define SPARSELU 1
#define BICGSTAB 2
int sol = SPARSELU;

#define COTW 1
#define MVW  2
int wei = MVW;

#define SAVE_VERTEX 1
#define SAVE_TEXCOORD 2
int sav = SAVE_VERTEX;

#if 0
// Boundary vertex is fixed or not
class Fixed {
public:
  Fixed(): isFixed_(false) {};
  ~Fixed(){};
  bool isFixed() const { return isFixed_; };
  void setIsFixed( bool f ) { isFixed_ = f; };

private:
  bool isFixed_;
};

OpenMesh::VPropHandleT<Fixed> ffs;
#endif

#define FIXED_BOUNDARY 1
#define NATURAL_BOUNDARY 2
int bou = FIXED_BOUNDARY;

// 0: vertex i, 1-N: its neighbor vertices j
class VVWeights {

public:

  VVWeights(){};
  ~VVWeights(){ w_.clear(); };

  double w( int i ) { return w_[i]; };
  int size() const { return w_.size(); };
  std::vector<double>& ws() { return w_; };
  void addWeight( double d ) { w_.push_back(d); };
  // for orbifold
  // void addWeight( double d, int pid ) {
  //   w_.push_back(d);
  //   param_id_.push_back( pid );
  // };
  void addParamID( MyMesh::VertexHandle& vh, int id ) {
    param_id_.insert( make_pair( vh, id ) );
  };
  int param_id( MyMesh::VertexHandle& vh ) { return param_id_[vh]; };
  bool setWeight( int i, double d ) {
    if ( i >= (int) w_.size() ) return false;
    w_[i] = d;
    return true;
  };
  
  void Print() {
    cout << "vvw" << endl;
    cout << "weight size " << size() << endl;
    for ( int i = 0; i < w_.size(); ++i )
      {
        cout << w_[i] << endl;
      }
    cout << "param_id " << endl;
    for ( std::map<MyMesh::VertexHandle,int>::iterator pi = param_id_.begin();
          pi != param_id_.end(); ++pi )
      {
        cout << "vt " << pi->first.idx() << " id " << pi->second << endl;
      }
  };

private:

  std::vector<double> w_;
  // parameter id for orbifold
  // for OF_BOUNDARY vertex, values are 1 (true), or -1 (false)
  std::map<MyMesh::VertexHandle,int> param_id_;
  
  

};

//////////////////////////////////////////////////////////////////////////////

#include "Orbifold.hxx"

class Param {

public:
  
  double cotw( MyMesh::Point& p0, MyMesh::Point& p1, MyMesh::Point& p2 ) {
    OMVector3d c1( p1 - p0 );
    OMVector3d c2( p2 - p0 );
    //double sin = OpenMesh::cross( c1, c2 ).length();
    double sin = c1.cross( c2 ).norm();
    //double angle = std::fabs( std::atan2(sin, OpenMesh::dot(c1, c2)) );
    double angle = std::fabs( std::atan2(sin, c1.dot(c2)) );

    // std::cout << "\tcot " << 1.0 / std::tan( angle ) << " tan " << std::tan( M_PI_2 - angle ) << std::endl;
    // if ( 1.0 / std::tan( angle ) - std::tan( M_PI_2 - angle ) > 1.0e-05 )
    //   std::cout << "diff!" << std::endl;
    return std::tan( M_PI_2 - angle );
  };

  double tan2w( MyMesh::Point& p0, MyMesh::Point& p1, MyMesh::Point& p2 ) {
    OMVector3d c1( p1 - p0 );
    OMVector3d c2( p2 - p0 );
    //double sin = OpenMesh::cross( c1, c2 ).length();
    double sin = c1.cross( c2 ).norm();
    //double angle = std::fabs( std::atan2(sin, OpenMesh::dot(c1, c2)) );
    double angle = std::fabs( std::atan2(sin, c1.dot(c2)) );

    // std::cout << "\tcot " << 1.0 / std::tan( angle ) << " tan " << std::tan( M_PI_2 - angle ) << std::endl;
    // if ( 1.0 / std::tan( angle ) - std::tan( M_PI_2 - angle ) > 1.0e-05 )
    //   std::cout << "diff!" << std::endl;
    return std::tan( angle / 2.0 );
  };

  //
  // A weight has to be doubled for natural boundary as described in [Karni 05]
  //
  double computeCotangentWeight( MyMesh& mesh, MyMesh::VertexIHalfedgeIter& vih_it ) {
    MyMesh::HalfedgeHandle iheh( *vih_it ); // vih_it.handle();
    // if ( mesh.from_vertex_handle( iheh ) == vh ) std::cout << "i from" << std::endl;
    // else if ( mesh.to_vertex_handle( iheh ) == vh ) std::cout << "i to" << std::endl;
    //MyMesh::HalfedgeHandle next_iheh = iheh.next();
    MyMesh::HalfedgeHandle next_iheh = mesh.next_halfedge_handle( *vih_it ); 
    // if ( mesh.from_vertex_handle( next_iheh ) == vh ) std::cout << "next i from" << std::endl;
    // else if ( mesh.to_vertex_handle( next_iheh ) == vh ) std::cout << "next i to" << std::endl

    // MyMesh::Point p0, p1, p2;
    auto p0 = mesh.point( mesh.to_vertex_handle( next_iheh ) );
    auto p1 = mesh.point( mesh.from_vertex_handle( iheh ) );
    auto p2 = mesh.point( mesh.to_vertex_handle( iheh ) );
    double al = cotw( p0, p1, p2 );

    // left face
    MyMesh::HalfedgeHandle mheh = mesh.opposite_halfedge_handle( iheh );
    // if ( mesh.from_vertex_handle( mheh ) == vh ) std::cout << "m from" << std::endl;
    // else if ( mesh.to_vertex_handle( mheh ) == vh ) std::cout << "m to" << std::endl;
    MyMesh::HalfedgeHandle next_mheh = mesh.next_halfedge_handle( mheh );

    p0 = mesh.point( mesh.to_vertex_handle( next_mheh ) );
    p1 = mesh.point( mesh.from_vertex_handle( mheh ) );
    p2 = mesh.point( mesh.to_vertex_handle( mheh ) );
    double be = cotw( p0, p1, p2 );

    //return (al + be) / 2.0;
    return (al + be);
  };

  double computeMeanValueWeight( MyMesh& mesh, MyMesh::VertexIHalfedgeIter& vih_it ) {
    //MyMesh::Point p0, p1, p2;

    // gamma
    MyMesh::HalfedgeHandle iheh( *vih_it ); // vih_it.handle();
    MyMesh::HalfedgeHandle next_iheh = mesh.next_halfedge_handle( iheh );
    auto p0 = mesh.point( mesh.from_vertex_handle( next_iheh ) );
    auto p1 = mesh.point( mesh.to_vertex_handle( next_iheh ) );
    auto p2 = mesh.point( mesh.from_vertex_handle( iheh ) );
    double ga = tan2w( p0, p1, p2 );
#if 0  
    std::cout << "\t gamma p0 " << mesh.from_vertex_handle( next_iheh ).idx() << " "
              << " p1 "  << mesh.to_vertex_handle( next_iheh ).idx() << " "
              << " p2 "  << mesh.from_vertex_handle( iheh ).idx() << std::endl;
#endif  
  
    // delta
    MyMesh::HalfedgeHandle mheh = mesh.opposite_halfedge_handle( iheh );
    MyMesh::HalfedgeHandle next_mheh = mesh.next_halfedge_handle( mheh );
    p0 = mesh.point( mesh.from_vertex_handle( mheh ) ); // i
    p1 = mesh.point( mesh.from_vertex_handle( next_mheh ) ); //j
    p2 = mesh.point( mesh.to_vertex_handle( next_mheh ) );
    double de = tan2w( p0, p1, p2 );
#if 0
    std::cout << "\t delta p0 " << mesh.from_vertex_handle( mheh ).idx() << " "
              << " p1 "  << mesh.from_vertex_handle( next_mheh ).idx() << " "
              << " p2 "  << mesh.to_vertex_handle( next_mheh ).idx() << std::endl;
#endif

    // return (ga + de) / (p0 - p1).length();
    return (ga + de) / (p0 - p1).norm();
  };

  void constructBoundary( MyMesh& mesh, std::vector<MyMesh::VertexHandle>& bverts ) {

    // find a boundary vertex
    MyMesh::VertexHandle vh0;
    for ( auto vh : mesh.vertices() )
      {
        if ( mesh.is_boundary( vh ) )
          {
            vh0 = vh; // v_it.handle();
            //std::cout << "boundary vertex " << vh0.idx() << std::endl;
            break;
          }
      }

#if 0
    MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
    for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
      {
        MyMesh::VertexHandle vh = *v_it;
        if ( mesh.is_boundary( vh ) )
          {
            vh0 = vh; // v_it.handle();
            //std::cout << "boundary vertex " << vh0.idx() << std::endl;
            break;
          }
      }
#endif

    std::vector<MyMesh::VertexHandle> bv_tmp;

    // construct boundary vertices
    MyMesh::VertexHandle vh = vh0, prev_vh;
    do {
      bv_tmp.push_back( vh );
      MyMesh::VertexVertexIter vv_it;
      for ( vv_it = mesh.vv_iter( vh ); vv_it.is_valid(); ++vv_it )
        {
          MyMesh::VertexHandle vc = *vv_it; // vv_it.handle();
          if ( mesh.is_boundary( vc ) && (prev_vh != vc) )
            {
              //std::cout << "\tfound! " << vc.idx() << std::endl;
              prev_vh = vh;
              vh = vc;
              break;
            }
        }
    } while ( vh != vh0 );
    std::cout << "#boundary vertices: " << bv_tmp.size() << std::endl;

#if 1
    // reverse direction
    bverts.push_back( bv_tmp[0] );
    for ( int i = bv_tmp.size()-1; i > 0; --i )
      bverts.push_back( bv_tmp[i] );
#endif
#if 0
    for ( int i = 0; i < bv_tmp.size(); ++i )
      bverts.push_back( bv_tmp[i] );
#endif
  };

  void computeBoundaryMapping( MyMesh& mesh,
                               std::vector<double>& paramx,  std::vector<double>& paramy ) {
    // property handle to store a flag
    mesh.add_property(fffs);
  
    // construct boundary vertices
    std::vector<MyMesh::VertexHandle> bverts;
    constructBoundary( mesh, bverts );

    // corner coordinates
    std::vector<double> cornerx(RECTANGLE);
    std::vector<double> cornery(RECTANGLE);
    std::vector<MyMesh::VertexHandle> cornerv;
    // fix two vertices of the coordinates (0,0) or (1,1)
    cornerx[0] = 0.0; cornerx[1] = 1.0; cornerx[2] = 1.0; cornerx[3] = 0.0;
    cornery[0] = 0.0; cornery[1] = 0.0; cornery[2] = 1.0; cornery[3] = 1.0;

    // compute boundary length
    double blength = 0.0;
    int bvn = bverts.size();
    for ( int i = 0; i < bvn-1; ++i )
      {
        // blength += (mesh.point( bverts[i] ) - mesh.point( bverts[i+1] )).length();
        blength += (mesh.point( bverts[i] ) - mesh.point( bverts[i+1] )).norm();
      }
    // blength += (mesh.point( bverts[bvn-1] ) - mesh.point( bverts[0] )).length();
    blength += (mesh.point( bverts[bvn-1] ) - mesh.point( bverts[0] )).norm();
    std::cout << "boundary length " << blength << std::endl;

    // determine corner vertices
    double l4 = blength / (double) RECTANGLE;
    double tl = 0.0;
    cornerv.push_back( bverts[0] );
    paramx[ bverts[0].idx() ] = cornerx[0];
    paramy[ bverts[0].idx() ] = cornery[0];
    // (0,0)
    Fixed& ff = mesh.property( fffs, bverts[0] );
    ff.setIsFixed( true );

    int count = 1;
    for ( int i = 0; i < bvn-1; ++i )
      {
        // tl += (mesh.point( bverts[i] ) - mesh.point( bverts[i+1] )).length();
        tl += (mesh.point( bverts[i] ) - mesh.point( bverts[i+1] )).norm();
        if ( (tl > l4) && (count < RECTANGLE) )
          {
            cornerv.push_back( bverts[i] );
            paramx[ bverts[i].idx() ] = cornerx[count];
            paramy[ bverts[i].idx() ] = cornery[count];
            if ( count == 2 ) // (1,1)
              {
                Fixed& ff = mesh.property( fffs, bverts[i] );
                ff.setIsFixed( true );
              }
            tl = 0.0;
            ++count;
          }
      }

    std::cout << "corner vertices: ";
    for ( int i = 0; i < 4; ++i ) std::cout << cornerv[i].idx() << " ";
    std::cout << std::endl;

    // compute inner boundary coordinates
    std::vector<MyMesh::VertexHandle> innerv;
    count = 1;
    for ( int i = 0; i < bvn; ++i )
      {
        innerv.push_back( bverts[i] );
        if ( cornerv[count] == bverts[i] )
          {
            // compute the length of poly line
            double pll = 0.0;
            for ( int j = 0; j < (int) innerv.size()-1; ++j )
              {
                // pll += (mesh.point( innerv[j] ) - mesh.point( bverts[j+1] )).length();
                pll += (mesh.point( innerv[j] ) - mesh.point( bverts[j+1] )).norm();
              }
            // compute the ratio of length and compute parameter
            double cll = 0.0;
            // std::cout << "count " << count << std::endl;
            for ( int j = 0; j < (int) innerv.size()-1; ++j )
              {
                double t = cll / pll;
                double px = (1 - t) * cornerx[count-1] + t * cornerx[count];
                double py = (1 - t) * cornery[count-1] + t * cornery[count];
                paramx[ innerv[j].idx() ] = px;
                paramy[ innerv[j].idx() ] = py;
                // std::cout << "\t t " << t << " param " << px << " " << py << std::endl;
                // cll += (mesh.point( innerv[j] ) - mesh.point( bverts[j+1] )).length();
                cll += (mesh.point( innerv[j] ) - mesh.point( bverts[j+1] )).norm();
              }

            // next poly line
            ++count;
            innerv.clear();
            innerv.push_back( bverts[i] );
          }
      }
    //
    // last poly line
    //
    innerv.push_back( bverts[0] );
    double pll = 0.0;
    for ( int j = 0; j < (int) innerv.size()-1; ++j )
      {
        // pll += (mesh.point( innerv[j] ) - mesh.point( bverts[j+1] )).length();
        pll += (mesh.point( innerv[j] ) - mesh.point( bverts[j+1] )).norm();
      }
    // compute the ratio of length and compute parameter
    double cll = 0.0;
    // std::cout << "count " << count << std::endl;
    for ( int j = 0; j < (int) innerv.size()-1; ++j )
      {
        double t = cll / pll;
        double px = (1 - t) * cornerx[RECTANGLE-1] + t * cornerx[0];
        double py = (1 - t) * cornery[RECTANGLE-1] + t * cornery[0];
        paramx[ innerv[j].idx() ] = px;
        paramy[ innerv[j].idx() ] = py;
        // std::cout << "\t t " << t << " param " << px << " " << py << std::endl;
        // cll += (mesh.point( innerv[j] ) - mesh.point( bverts[j+1] )).length();
        cll += (mesh.point( innerv[j] ) - mesh.point( bverts[j+1] )).norm();
      }
  };

  ////////////////////////////////////////////////////////////////////////////////////////

  void applyParam_FixedBoundary( MyMesh& mesh,
                                 std::vector<double>& paramx, std::vector<double>& paramy,
                                 int sol_, int wei_ ) {
    int n_vt = mesh.n_vertices();

    //
    // boundary mapping
    // compute 2D parameters for boundary
    //
    computeBoundaryMapping( mesh, paramx, paramy );

    // property handle to store weights
    OpenMesh::VPropHandleT<VVWeights> vvws;
    mesh.add_property(vvws);

    //
    // compute weights
    //
    MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
    for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
      {
        // for check
        //std::cout << i << " " << v_it.handle().idx() << std::endl;
        //if ( i != v_it.handle().idx() ) std::cout << "ng" << std::endl;

        VVWeights& vvw = mesh.property(vvws,*v_it);
        double wd = 0.0;
        vvw.addWeight( wd );
        //MyMesh::VertexHandle vh = v_it.handle();
        //std::cout << "v " << i << std::endl;

#if 0
        MyMesh::VertexHandle vh = *v_it;
        std::cout << "i = " << vh.idx() << std::endl;
#endif
        MyMesh::VertexIHalfedgeIter vih_it;
        for (vih_it=mesh.vih_iter( *v_it ); vih_it.is_valid(); ++vih_it)
          {
            double wdc;
            if ( wei_ == COTW )
              wdc = computeCotangentWeight( mesh, vih_it );
            else if ( wei_ == MVW )
              wdc = computeMeanValueWeight( mesh, vih_it );
            vvw.addWeight( wdc );
            wd -= wdc;
          }

        vvw.setWeight( 0, wd );

#if 0
        std::vector<double>& ws = vvw.ws();
        for ( int j = 0; j < ws.size(); ++j )
          std::cout << "\t " << ws[j] << std::endl;
        std::cout << std::endl;
#endif
      }

    //
    // compute inner parameters
    //
    // build up sparse matrix A
    std::cout << "Setup sparse matrix ..." << std::endl;
    std::vector<Eigen::Triplet<double> > tripletList;
    for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
      {
        MyMesh::VertexHandle vh = *v_it;
        int i = vh.idx();
        if ( !(mesh.is_boundary(*v_it)) )
          {
            VVWeights& vvw = mesh.property(vvws,*v_it);
            tripletList.push_back( Eigen::Triplet<double>(i, i, vvw.w(0)) );
            int k;
            MyMesh::VertexVertexIter vv_it;
            for (k=1, vv_it=mesh.vv_iter( *v_it ); vv_it.is_valid(); ++vv_it, ++k)
              {
                MyMesh::VertexHandle vvh = *vv_it;
                int j = vvh.idx();
                tripletList.push_back( Eigen::Triplet<double>(i, j, vvw.w(k)) );
              }
          }
        else
          tripletList.push_back( Eigen::Triplet<double>(i, i, 1.0) );
      }

    Eigen::SparseMatrix<double> spmat( n_vt, n_vt );
    spmat.setFromTriplets( tripletList.begin(), tripletList.end() );
    spmat.makeCompressed();

    // setup vector b
    Eigen::VectorXd xx(n_vt), xy(n_vt), bx(n_vt), by(n_vt);
    for ( int i = 0; i < n_vt; ++i )
      {
        bx[i] = paramx[i];
        by[i] = paramy[i];
      }

    // to symmetric positive-definite matrix
    // toSpd( spmat, bverts, bx, by );

    // solve x
    if ( sol_ == BICGSTAB )
      {
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver(spmat);
        std::cout << "solve xcoords ..." << std::endl;
        xx = solver.solve(bx);
        std::cout << "solve ycoords ..." << std::endl;
        xy = solver.solve(by);
      }
    else if ( sol_ == SPARSELU )
      {
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;

        solver.analyzePattern(spmat); 
        // Compute the numerical factorization 
        solver.factorize(spmat); 
        //Use the factors to solve the linear system 
        std::cout << "solve xcoords ..." << std::endl;
        xx = solver.solve(bx);
        std::cout << "solve ycoords ..." << std::endl;
        xy = solver.solve(by);
      }
    for ( int i = 0; i < n_vt; ++i )
      {
        paramx[i] = xx[i];
        paramy[i] = xy[i];
      }

#if 0
    for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
      {
        MyMesh::VertexHandle vh = *v_it;
        int i = vh.idx();
        Fixed& ff = mesh.property( fffs, vh );
        if ( ff.isFixed() )
          {
            std::cout << "i " << i << " param " << paramx[i] << " " << paramy[i] << std::endl;
          }
      }
#endif
  
    mesh.remove_property(vvws);

    std::cout << "done." << std::endl;
  };

  void applyParam_NaturalBoundary( MyMesh& mesh,
                                   std::vector<double>& paramx, std::vector<double>& paramy,
                                   int sol_, int wei_ ) {
    const int n_vt = mesh.n_vertices();

#if 0
    //
    // boundary mapping
    //
    // compute 2D parameters for boundary
    computeBoundaryMapping( mesh, paramx, paramy );
#endif

    MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
#if 1
    for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
      {
        MyMesh::VertexHandle vh = *v_it;
        //if ( !(mesh.is_boundary(vh)) ) 
        Fixed& ff = mesh.property( fffs, vh );
        if ( !(ff.isFixed()) ) // not fixed
          {
            paramx[vh.idx()] = paramy[vh.idx()] = 0.0;
          }
      }
#endif

    // property handle to store weights
    OpenMesh::VPropHandleT<VVWeights> vvws;
    mesh.add_property(vvws);

    //
    // compute weights
    //
    //MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
    for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
      {
        // for check
        //std::cout << i << " " << v_it.handle().idx() << std::endl;
        //if ( i != v_it.handle().idx() ) std::cout << "ng" << std::endl;

        VVWeights& vvw = mesh.property(vvws,*v_it);
        double wd = 0.0;
        vvw.addWeight( wd );
        //MyMesh::VertexHandle vh = v_it.handle();
        //std::cout << "v " << i << std::endl;

        MyMesh::VertexIHalfedgeIter vih_it;
        for (vih_it=mesh.vih_iter( *v_it ); vih_it.is_valid(); ++vih_it)
          {
            double wdc;
            if ( wei_ == COTW )
              wdc = computeCotangentWeight( mesh, vih_it );
            else if ( wei_ == MVW )
              wdc = computeMeanValueWeight( mesh, vih_it );
            vvw.addWeight( wdc );
            wd -= wdc;
          }

        vvw.setWeight( 0, wd );

#if 0
        std::vector<double>& ws = vvw.ws();
        for ( int j = 0; j < ws.size(); ++j )
          std::cout << "\t " << ws[j] << std::endl;
        std::cout << std::endl;
#endif
      }

    //
    // compute inner parameters
    //
    // build up sparse matrix A
    std::cout << "Setup sparse matrix ..." << std::endl;
    std::vector<Eigen::Triplet<double> > tripletList;
    for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
      {
        MyMesh::VertexHandle vh = *v_it;
        int i = vh.idx();
        Fixed& ff = mesh.property( fffs, vh );
        if ( !(ff.isFixed()) ) // not fixed
          {
            VVWeights& vvw = mesh.property(vvws,*v_it);

            // set for x coord
            tripletList.push_back( Eigen::Triplet<double>(i, i, vvw.w(0)) );
            // set for y coord
            tripletList.push_back( Eigen::Triplet<double>(i+n_vt, i+n_vt, vvw.w(0)) );
            int k;
            MyMesh::VertexVertexIter vv_it;
            for (k=1, vv_it=mesh.vv_iter( *v_it ); vv_it.is_valid(); ++vv_it, ++k)
              {
                MyMesh::VertexHandle vvh = *vv_it;
                int j = vvh.idx();
                // set for x coord
                tripletList.push_back( Eigen::Triplet<double>(i, j, vvw.w(k)) );
                // set for y coord
                tripletList.push_back( Eigen::Triplet<double>(i+n_vt, j+n_vt, vvw.w(k)) );
              }

            // boundary vertices
            if ( mesh.is_boundary(*v_it) )
              {
#if 0
                std::cout << "id " << i << std::endl;
                for ( int j = 0; j < vvw.size(); ++j )
                  {
                    std::cout << "\t" << vvw.w(j) << std::endl;
                  }
#endif
                // get neighbor boundary vertices
                MyMesh::VertexVertexIter vv_it;
                vv_it=mesh.vv_iter( *v_it );
                MyMesh::VertexHandle sbvh = *vv_it, ebvh;
                while ( vv_it.is_valid() )
                  {
                    ebvh = *vv_it;
                    ++vv_it;
                  }
                int i1 = sbvh.idx();
                int i2 = ebvh.idx();

                if ( wei == COTW )
                  {
                    // std::cout << "i1 " << i1 << " " << mesh.is_boundary( sbvh ) << std::endl;
                    // std::cout << "i2 " << i2 << " " << mesh.is_boundary( ebvh ) << std::endl;

                    // for x coords
                    tripletList.push_back( Eigen::Triplet<double>(i, i2+n_vt, 1.0) );
                    tripletList.push_back( Eigen::Triplet<double>(i, i1+n_vt, -1.0) );
                    // tripletList.push_back( Eigen::Triplet<double>(i, i2+n_vt, -1.0) );
                    // tripletList.push_back( Eigen::Triplet<double>(i, i1+n_vt, 1.0) );

                    // for y coords
                    tripletList.push_back( Eigen::Triplet<double>(i+n_vt, i1, 1.0) );
                    tripletList.push_back( Eigen::Triplet<double>(i+n_vt, i2, -1.0) );
                    // tripletList.push_back( Eigen::Triplet<double>(i+n_vt, i1, -1.0) );
                    // tripletList.push_back( Eigen::Triplet<double>(i+n_vt, i2, 1.0) );
                  }
                else if ( wei == MVW )
                  {
                    MyMesh::Point p0 = mesh.point( *v_it );
                    MyMesh::Point p1 = mesh.point( sbvh );
                    MyMesh::Point p2 = mesh.point( ebvh );
                    // double ir1 = 1.0/((p1 - p0).length());
                    double ir1 = 1.0/((p1 - p0).norm());
                    // double ir2 = 1.0/((p2 - p0).length());
                    double ir2 = 1.0/((p2 - p0).norm());
                    // std::cout << ir1 << " " << ir2 << std::endl;
                    // for x coords
                    tripletList.push_back( Eigen::Triplet<double>(i, i2+n_vt, ir2) );
                    tripletList.push_back( Eigen::Triplet<double>(i, i1+n_vt, -ir1) );
                    tripletList.push_back( Eigen::Triplet<double>(i, i+n_vt, ir1 - ir2) );
                    // for y coords
                    tripletList.push_back( Eigen::Triplet<double>(i+n_vt, i1, ir1) );
                    tripletList.push_back( Eigen::Triplet<double>(i+n_vt, i2, -ir2) );
                    tripletList.push_back( Eigen::Triplet<double>(i+n_vt, i, ir2 - ir1) );
                  }
              }
          }
        else // fixed
          {
            // std::cout << "i = " << i << " param " << paramx[i] << " " << paramy[i] << std::endl;
            tripletList.push_back( Eigen::Triplet<double>(i, i, 1.0) );
            tripletList.push_back( Eigen::Triplet<double>(i+n_vt, i+n_vt, 1.0) );
          }
      }

    Eigen::SparseMatrix<double> spmat( 2*n_vt, 2*n_vt );
    spmat.setFromTriplets( tripletList.begin(), tripletList.end() );
    spmat.makeCompressed();

    // setup vector b
    //Eigen::VectorXd xx(n_vt), xy(n_vt), bx(n_vt), by(n_vt);
    Eigen::VectorXd xx(2*n_vt), bx(2*n_vt);
    for ( int i = 0; i < n_vt; ++i )
      {
        // set for x coord
        bx[i] = paramx[i];
        // set for y coord
        bx[i+n_vt] = paramy[i];
      }

    // solve x
    if ( sol_ == BICGSTAB )
      {
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver(spmat);
        std::cout << "solve x and y coords ..." << std::endl;
        xx = solver.solve(bx);
        // std::cout << "solve ycoords ..." << std::endl;
        // xy = solver.solve(by);
      }
    else if ( sol_ == SPARSELU )
      {
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;

        solver.analyzePattern(spmat); 
        // Compute the numerical factorization 
        solver.factorize(spmat); 
        //Use the factors to solve the linear system 
        std::cout << "solve x and y coords ..." << std::endl;
        xx = solver.solve(bx);
        // std::cout << "solve ycoords ..." << std::endl;
        // xy = solver.solve(by);
      }
    for ( int i = 0; i < n_vt; ++i )
      {
        paramx[i] = xx[i];
        paramy[i] = xx[i+n_vt];
      }

#if 0
    for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
      {
        MyMesh::VertexHandle vh = *v_it;
        int i = vh.idx();
        Fixed& ff = mesh.property( fffs, vh );
        if ( ff.isFixed() )
          {
            std::cout << "i " << i << " param " << paramx[i] << " " << paramy[i] << std::endl;
          }
      }
#endif

    mesh.remove_property(vvws);

    std::cout << "done." << std::endl;
  };

  MyMesh::VertexIHalfedgeCWIter start_vih_cwiter( MyMesh::VertexHandle& vh,
                                                  MyMesh::VertexHandle& pvh,
                                                  MyMesh& mesh ) {
    int c = 0;
    MyMesh::VertexIHalfedgeCWIter vih_cwit = mesh.vih_cwiter( vh );
    while (1)
      {
        if ( mesh.from_vertex_handle( *vih_cwit ) == pvh )
          return vih_cwit;
        ++vih_cwit; ++c;
      }
  };

  MyMesh::VertexOHalfedgeCWIter start_voh_cwiter( MyMesh::VertexHandle& vh,
                                                  MyMesh::VertexHandle& pvh,
                                                  MyMesh& mesh ) {
    int c = 0;
    MyMesh::VertexOHalfedgeCWIter voh_cwit = mesh.voh_cwiter( vh );
    while (1)
      {
        if ( mesh.to_vertex_handle( *voh_cwit ) == pvh )
          return voh_cwit;
        ++voh_cwit; ++c;
      }
  };

  MyMesh::VertexOHalfedgeCCWIter start_voh_ccwiter( MyMesh::VertexHandle& vh,
                                                    MyMesh::VertexHandle& pvh,
                                                    MyMesh& mesh ) {
    int c = 0;
    MyMesh::VertexOHalfedgeCCWIter voh_ccwit = mesh.voh_ccwiter( vh );
    while (1)
      {
        if ( mesh.to_vertex_handle( *voh_ccwit ) == pvh )
          return voh_ccwit;
        ++voh_ccwit; ++c;
      }
  };

  MyMesh::VertexVertexIter start_vv_iter( MyMesh::VertexVertexIter vv_it,
                                          MyMesh::VertexHandle& pvt ) {
    int c = 0;
    while (1)
      {
        if ( *vv_it == pvt )
          {
            // cout << "count " << c << endl;
            return vv_it;
          }
        ++vv_it; ++c;
      }
  };

  //
  // set param_id to vvw of 1-ring neighbor vertices for boundary vertices
  //
  void setParamID_I( MyMesh::VertexHandle& vh,
                     VVWeights& vvw,
                     std::map<MyMesh::VertexHandle,int>& alt_param_id,
                     MyMesh& mesh,
                     std::vector<unsigned int>& vertex_type,
                     bool reverse ) {
    // if vh is OF_INTERNAL, apply add param_id to only 1-ring neighbor points
    // with OF_BOUNDARY or OF_CONESINGULARITY_PI2.
    // cout << "internal vt " << vh.idx() << " 1-ring neighbor start" << endl;
    for ( MyMesh::VertexVertexIter vv_it = mesh.vv_iter( vh ); vv_it.is_valid(); ++vv_it )
      {
        MyMesh::VertexHandle vvh = *vv_it;
        int id = vvh.idx();
        int rev_id = alt_param_id[vvh];
        // cout << "\tvertex " << id << endl;
        if ( (vertex_type[ id ] == OF_BOUNDARY) ||
             (vertex_type[ id ] == OF_CONESINGULARITY_PI2)  )
          {
            // if a 1-ring neighbor vertex is OF_BOUNDARY or OF_CONESINGULARITY_PI2,
            // set param_id to its original id or alt_param_id
            // TODO: 11/29 ここのところもう一度見直す
            if ( reverse == false )
              {
                // 11/29 下と交換
                // vvw.addParamID( vvh, id );
                vvw.addParamID( vvh, rev_id );
                // cout << "\tvt " << id << " type " << vertex_type[id] << " reverse " << reverse << " param " << rev_id << endl;
              }
            else
              {
                // 11/29 上と交換
                // vvw.addParamID( vvh, rev_id );
                vvw.addParamID( vvh, id );
                // cout << "\tvt " << id << " type " << vertex_type[id] << " reverse " << reverse << " param " << id << endl;
              }
          }
      }
    // cout << "internal 1-ring neighbor end" << endl;
  };

  void setParamID_B( MyMesh::VertexHandle& vh,
                     MyMesh::VertexHandle& pvh,
                     MyMesh::VertexHandle& nvh,
                     VVWeights& vvw,
                     std::map<MyMesh::VertexHandle,int>& alt_param_id,
                     MyMesh& mesh,
                     std::vector<unsigned int>& vertex_type,
                     bool reverse ) {
    // if vh is OF_BOUNDARY or OF_CONESINGULARITY_PI2, apply add param_id
    /// to all 1-ring neighbor points.
    // TODO: reverse と rev の関係について見直す
    // cout << "boundary 1-ring neighbor start" << endl;
    bool rev = reverse;
    MyMesh::VertexVertexIter vv_it = start_vv_iter( mesh.vv_iter( vh ), pvh );
    MyMesh::VertexHandle svh = *vv_it;
    // int i;
    // for ( i = 0, vv_it = mesh.vv_iter( vh ); vv_it.is_valid(); ++vv_it, ++i ) {
    do {
      MyMesh::VertexHandle vvh = *vv_it;
      int id = vvh.idx();

      // cout << "vertex " << id << endl;

      // set param_id to 1.0 (reverse=true) or -1.0 (reverse=false)
      if ( rev == false )
        {
          // 11/30 下と交換
          // vvw.addParamID( vvh, -1 );
          vvw.addParamID( vvh, 1 );
          // cout << "\tvt " << id << " type " << vertex_type[id] << " reverse " << rev << " param " << 1.0 << endl;
        }
      else // rev == true
        {
          // 11/30 上と交換
          // vvw.addParamID( vvh, 1 );
          vvw.addParamID( vvh, -1 );
          // cout << "\tvt " << id << " type " << vertex_type[id] << " reverse " << rev << " param " << -1.0 << endl;
        }
      ++vv_it;
      if ( !(vv_it.is_valid()) ) vv_it = mesh.vv_iter( vh ); // reset vv_it to start
      if ( *vv_it == nvh )
        {
          rev = ( rev == true ) ? false : true; // parameter is changed from nvt
        }
    } while ( *vv_it != svh );
    // cout << "boundary 1-ring neighbor end" << endl;
  };

  int applyParam_Orbifold( MyMesh& mesh,
                           Orbifold& orbi,
                           // std::vector<MyMesh::VertexHandle>& cs_vertices,
                           std::vector<double>& paramx, std::vector<double>& paramy,
                           int sol_, int wei_ ) {

    cout << "compute orbifold parameterization ... " << endl;

    int n_vt = mesh.n_vertices();

    // set orbifold boundary
    orbi.calcBoundaries();

    // cone singularity vertices
    std::vector<MyMesh::VertexHandle>& cs_vertices = orbi.cs_vertices();

    // parameters setup
    std::vector<std::vector<MyMesh::VertexHandle> >& path = orbi.path();

    // the number of parameter vertices
    // - mesh.n_vertices() + path[2].size()-2 + path[3].size()-2 + 1 (cs_vertex[3])
    // - two parameters in each vertex of type OF_BOUNDARY or OF_CONESINGULARITY_PI2
    static int n_param = mesh.n_vertices() + path[2].size()-2 + path[3].size()-2 + 1;
    // cout << "vertices " << mesh.n_vertices() << " parameters " << n_param << endl;
    paramx.resize( n_param );
    paramy.resize( n_param );

    // cout << "0----------------------------------------------------------" << endl;
    // set a param_id pair and a center for OF_BOUNDARY vertex
    std::map<MyMesh::VertexHandle,int> alt_param_id;
    int param_id = mesh.n_vertices(); // alt_param_id starts from mesh.n_vertices()
    for ( int i = 0; i < path.size() / 2; ++i ) // search for only two paths
      {
        // cone singularities are not concerned
        for ( int j = 1; j < path[i].size() - 1; ++j )
          {
            MyMesh::VertexHandle vt = path[i][j];
            alt_param_id.insert( make_pair(vt, param_id) );
            // cout << "\t\t " << vt.idx() << " " << alt_param_id[vt] << endl;
            ++param_id;
          }
      }
    // set a param_id pair for  OF_CONESINGULARITY_PI2 vertex
    // alt parameter is last one (n_param - 1)
    MyMesh::VertexHandle cs_vt = cs_vertices[1];
    alt_param_id.insert( make_pair(cs_vt, param_id) );
    ++param_id;

    // cout << "parameters " << n_param << " counts " << param_id << endl;

    // property handle to store weights
    OpenMesh::VPropHandleT<VVWeights> vvws;
    mesh.add_property(vvws);

    std::vector<unsigned int>& vertex_type = orbi.vertex_type();

    // cout << "1----------------------------------------------------------" << endl;
    // make consistent path from path no.0 and no.1
    std::vector<MyMesh::VertexHandle> path01;
    for ( int i = 0; i < path.size() / 2; ++i )
      for ( int j = 0; j < path[i].size(); ++j )
        if ( (i != 1) || (j != 0) ) path01.push_back( path[i][j] );

    // determine param_id for OF_BOUNDARY and OF_CONESINGULARITY_PI2 vertices of 1-ring neighbors
    for ( int i = 1; i < path01.size() - 1; ++i )
      {
    // for ( int i = 0; i < path.size()/2; ++i )
    //   {
    //     for ( int j = 1; j < path[i].size(); ++j )
    //       {
        // MyMesh::VertexHandle vh = path[i][j];    // current boundary vertex
        // MyMesh::VertexHandle pvh = path[i][j-1]; // prev boundary vertex
        // MyMesh::VertexHandle nvh = path[i][j+1]; // next boundary vertex
        MyMesh::VertexHandle vh = path01[i];    // current boundary vertex
        MyMesh::VertexHandle pvh = path01[i-1]; // prev boundary vertex
        MyMesh::VertexHandle nvh = path01[i+1]; // next boundary vertex

        // bool reverse = ( i > 1 ) ? true : false;
        // Since i < 2, reverse is always false.
        bool reverse = false;

        // processing OF_BOUNDARY point
        // cout << "(a) boundary vt " << vh.idx() << endl;
        if ( vertex_type[vh.idx()] == OF_BOUNDARY )
          {
            VVWeights& vvw = mesh.property(vvws, vh);
            setParamID_B( vh, pvh, nvh, vvw, alt_param_id, mesh, vertex_type, reverse );
            // vvw.Print();
          }

        // 1-ring neighbor vertices
        // start vv_it is set to prev boundary vertex
        MyMesh::VertexVertexIter vv_it = start_vv_iter( mesh.vv_iter( vh ), pvh );
        // cout << "vt " << vt.idx() << " pvt " << pvt.idx() << " start " << svt.idx() << endl;
        MyMesh::VertexHandle svh = *vv_it;
        // cout << "svh " << svh.idx() << endl;
        do {
          // processing 1-ring neighbors of OF_BOUNDARY point
          MyMesh::VertexHandle vvh = *vv_it;
          // if ( (vvh != pvh) && (vvh != nvh) )
          if ( vertex_type[vvh.idx()] == OF_INTERNAL )
            {
              // cout << "(b) should be internal vt " << vvh.idx() << endl;
              VVWeights& vvvw = mesh.property(vvws, vvh);
              setParamID_I( vvh, vvvw, alt_param_id, mesh, vertex_type, reverse );
            }

          ++vv_it;
          if ( !(vv_it.is_valid()) ) vv_it = mesh.vv_iter( vh ); // reset vv_it to start
          if ( *vv_it == nvh )
            {
              // cout << "evh " << nvh.idx() << endl;
              reverse = ( reverse == true ) ? false : true; // parameter is changed from nvt
            }
        } while ( *vv_it != svh );
        // }
      }
    // cout << "2----------------------------------------------------------" << endl;

    // set cone sigularities fixed
    //std::vector<MyMesh::VertexHandle>& cs_vertices = orbi.cs_vertices();
    int i0 = cs_vertices[0].idx();
    int i1 = cs_vertices[1].idx();
    int i2 = cs_vertices[2].idx();
    int i3 = alt_param_id[cs_vertices[1]];
    paramx[i0] = 0.0; paramy[i0] = 0.0;
    paramx[i1] = 1.0; paramy[i1] = 0.0;
    paramx[i2] = 1.0; paramy[i2] = 1.0;
    paramx[i3] = 0.0; paramy[i3] = 1.0;

    // cout << "cs " << cs_vertices[0].idx() << " " << cs_vertices[1].idx() << " "
    //      << cs_vertices[2].idx() << " " << n_param - 1 << endl;

    //
    // compute weights
    //
    MyMesh::VertexIter v_it, v_end(mesh.vertices_end());
    for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
      {
        // for check
        //std::cout << i << " " << v_it.handle().idx() << std::endl;
        //if ( i != v_it.handle().idx() ) std::cout << "ng" << std::endl;

        VVWeights& vvw = mesh.property(vvws,*v_it);
        double wd = 0.0;
        vvw.addWeight( wd );
        //MyMesh::VertexHandle vh = v_it.handle();
        //std::cout << "v " << i << std::endl;

        MyMesh::VertexIHalfedgeIter vih_it;
        for (vih_it=mesh.vih_iter( *v_it ); vih_it.is_valid(); ++vih_it)
          {
            double wdc;
            if ( wei_ == COTW )
              wdc = computeCotangentWeight( mesh, vih_it );
            else if ( wei_ == MVW )
              wdc = computeMeanValueWeight( mesh, vih_it );
            vvw.addWeight( wdc );
            wd -= wdc;
          }

        vvw.setWeight( 0, wd );

#if 0
        std::vector<double>& ws = vvw.ws();
        for ( int j = 0; j < ws.size(); ++j )
          std::cout << "\t " << ws[j] << std::endl;
        std::cout << std::endl;
#endif
      }

    //
    // compute inner parameters
    //
    // build up sparse matrix A
    std::cout << "Setup sparse matrix ..." << std::endl;
    std::vector<Eigen::Triplet<double> > tripletList;
    for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
      {
        MyMesh::VertexHandle vh = *v_it;
        int i = vh.idx();
        if ( vertex_type[i] == OF_INTERNAL ) // internal vertex
          {
            VVWeights& vvw = mesh.property(vvws,*v_it);

            // set for x coord
            tripletList.push_back( Eigen::Triplet<double>(i, i, vvw.w(0)) );
            // set for y coord
            tripletList.push_back( Eigen::Triplet<double>(i+n_param, i+n_param, vvw.w(0)) );

            int k;
            MyMesh::VertexVertexIter vv_it;
            // cout << "internal vt " << i << endl;
            for (k=1, vv_it=mesh.vv_iter( *v_it ); vv_it.is_valid(); ++vv_it, ++k)
              {
                MyMesh::VertexHandle vvh = *vv_it;
                int id = vvh.idx();
                int j;
                if ( (vertex_type[id] == OF_BOUNDARY) ||
                     (vertex_type[id] == OF_CONESINGULARITY_PI2))
                  j = vvw.param_id( vvh );
                else
                  j = id;

                // if ( (vertex_type[id] == OF_BOUNDARY) ||
                //      (vertex_type[id] == OF_CONESINGULARITY_PI2))
                // cout << "\tvt " << id << " type " << vertex_type[id] << " param id " << j << " w " << vvw.w(k) << endl;

                // set for x coord
                tripletList.push_back( Eigen::Triplet<double>(i, j, vvw.w(k)) );
                // set for y coord
                tripletList.push_back( Eigen::Triplet<double>(i+n_param, j+n_param, vvw.w(k)) );
              }
          }
        else if ( vertex_type[i] == OF_CONESINGULARITY_PI ) // fixed
          {
            // original parameter id for cs_vertices[0] or cs_vertices[2]
            // set for x coord
            tripletList.push_back( Eigen::Triplet<double>(i, i, 1.0) );
            // set for y coord
            tripletList.push_back( Eigen::Triplet<double>(i+n_param, i+n_param, 1.0) );
          }
        else if ( vertex_type[i] == OF_CONESINGULARITY_PI2 ) // fixed
          {
            // original parameter id for cs_vertices[1]
            // set for x coord
            tripletList.push_back( Eigen::Triplet<double>(i, i, 1.0) );
            // set for y coord
            tripletList.push_back( Eigen::Triplet<double>(i+n_param, i+n_param, 1.0) );
            
            // altanative parameter id for cs_vertices[3]
            int id = alt_param_id[vh];
            // set for x coord
            tripletList.push_back( Eigen::Triplet<double>(id, id, 1.0) );
            // set for y coord
            tripletList.push_back( Eigen::Triplet<double>(id+n_param, id+n_param, 1.0) );
          }
        else // vertex_type[i] == OF_BOUNDARY
          {
#if 0
            cout << "boundary vt " << i << endl;

            VVWeights& vvw = mesh.property(vvws, *v_it);
            vvw.Print();

            //tripletList.push_back( Eigen::Triplet<double>(i, i, vvw.w(0)) );
            double sum_wei = 0.0;  // weight for non-reverse 
            double sum_wei_rev = 0.0; // weight for reverse
            int rev_i = alt_param_id[ vh ];

            //
            // eq.9a
            //
            // R(90)
            int k;
            MyMesh::VertexVertexIter vv_it;
            for (k=1, vv_it=mesh.vv_iter( *v_it ); vv_it.is_valid(); ++vv_it, ++k)
              {
                MyMesh::VertexHandle vvh = *vv_it;
                int j = vvh.idx();
                int rev = vvw.param_id( vvh );
                // cout << "vt id " << j << " reverse " << rev << endl;
                // classify weights to two according to the reverse flag
                if ( rev == -1 ) // non-reverse
                  {
                    sum_wei -= vvw.w(k);
                    // set for x coord
                    tripletList.push_back( Eigen::Triplet<double>(i, j, vvw.w(k)) );
                    // set for y coord
                    tripletList.push_back( Eigen::Triplet<double>(i+n_param, j+n_param, vvw.w(k)) );
                  }
                else // reverse
                  {
                    sum_wei_rev -= vvw.w(k);
                    // - wij * y_j
                    // set for x coord
                    int jj;
                    if ( (vertex_type[j] == OF_BOUNDARY) ||
                         (vertex_type[j] == OF_CONESINGULARITY_PI2) )
                      jj = alt_param_id[vvh];
                    else // internal
                      jj = j;
                    tripletList.push_back( Eigen::Triplet<double>(i, jj+n_param, -vvw.w(k)) );
                    // + wij * x_j
                    // set for y coord
                    tripletList.push_back( Eigen::Triplet<double>(i+n_param, jj, vvw.w(k)) );
                  }
              }
            // non-reverse weight
            // + sum wij * x_i
            // set for x coord
            tripletList.push_back( Eigen::Triplet<double>(i, i, sum_wei) );
            // + sum wij * y_i
            // set for y coord
            tripletList.push_back( Eigen::Triplet<double>(i+n_param, i+n_param, sum_wei) );

            // reverse weight
            // + sum wij * y_i
            // set for x coord
            tripletList.push_back( Eigen::Triplet<double>(i, rev_i+n_param, -sum_wei_rev) );
            // - sum wij * x_i
            // set for y coord
            tripletList.push_back( Eigen::Triplet<double>(i+n_param, rev_i, sum_wei_rev) );

            cout << "i " << i << " rev_i " << rev_i << " wei " << sum_wei << " wei-rev " << sum_wei_rev << " wei (original) " << vvw.w(0) << endl;

#endif
          }
      }

    //
    // build up sparse matrix A (cont.)
    // eq.9a
    //
    // k = 0: path no.0, k = 1: path no.2
    for ( int k = 0; k < 2; ++k )
      {
        int l = (k == 0) ? 0 : 2;
        for ( int j = 1; j < path[l].size() - 1; ++j )
          {
            MyMesh::VertexHandle vh = path[l][j];
            int i = vh.idx();
            int rev_i = alt_param_id[vh];

            // cout << "boundary vt " << i << endl;

            VVWeights& vvw = mesh.property(vvws, vh);
            // vvw.Print();

            double sum_wei = 0.0;  // weight for non-reverse 
            double sum_wei_rev = 0.0; // weight for reverse

            int m;
            MyMesh::VertexVertexIter vv_it;
            for (m=1, vv_it=mesh.vv_iter( vh ); vv_it.is_valid(); ++vv_it, ++m)
              {
                MyMesh::VertexHandle vvh = *vv_it;
                int j = vvh.idx();
                int rev = vvw.param_id( vvh );

                // classify weights to two according to the reverse flag
                if ( rev == -1 ) // non-reverse
                  {
                    sum_wei -= vvw.w(m);
                    // set for x coord
                    tripletList.push_back( Eigen::Triplet<double>(i, j, vvw.w(m)) );
                    // set for y coord
                    tripletList.push_back( Eigen::Triplet<double>(i+n_param, j+n_param, vvw.w(m)) );
                  }
                else // reverse
                  {
                    sum_wei_rev -= vvw.w(m);
                    int jj;
                    if ( (vertex_type[j] == OF_BOUNDARY) ||
                         (vertex_type[j] == OF_CONESINGULARITY_PI2) )
                      jj = alt_param_id[vvh];
                    else // internal
                      jj = j;

                    if ( k == 0 )
                      {
                        // R(-90)
                        // - wij * y_j
                        // set for x coord
                        tripletList.push_back( Eigen::Triplet<double>(i, jj+n_param, vvw.w(m)) );
                        // + wij * x_j
                        // set for y coord
                        tripletList.push_back( Eigen::Triplet<double>(i+n_param, jj, -vvw.w(m)) );
                      }
                    else
                      {
                        // R(90)
                        // - wij * y_j
                        // set for x coord
                        tripletList.push_back( Eigen::Triplet<double>(i, jj+n_param, -vvw.w(m)) );
                        // + wij * x_j
                        // set for y coord
                        tripletList.push_back( Eigen::Triplet<double>(i+n_param, jj, vvw.w(m)) );
                      }
                  }
              }
            // non-reverse weight
            // + sum wij * x_i
            // set for x coord
            tripletList.push_back( Eigen::Triplet<double>(i, i, sum_wei) );
            // + sum wij * y_i
            // set for y coord
            tripletList.push_back( Eigen::Triplet<double>(i+n_param, i+n_param, sum_wei) );

            // reverse weight
            if ( k == 0 )
              {
                // R(-90)
                // + sum wij * y_i
                // set for x coord
                tripletList.push_back( Eigen::Triplet<double>(i, rev_i+n_param, sum_wei_rev) );
                // - sum wij * x_i
                // set for y coord
                tripletList.push_back( Eigen::Triplet<double>(i+n_param, rev_i, -sum_wei_rev) );
              }
            else
              {
                // R(90)
                // + sum wij * y_i
                // set for x coord
                tripletList.push_back( Eigen::Triplet<double>(i, rev_i+n_param, -sum_wei_rev) );
                // - sum wij * x_i
                // set for y coord
                tripletList.push_back( Eigen::Triplet<double>(i+n_param, rev_i, sum_wei_rev) );
              }

            // cout << "i " << i << " rev_i " << rev_i << " wei " << sum_wei << " wei-rev " << sum_wei_rev << " wei (original) " << vvw.w(0) << endl;

          }
      }

    //
    // build up sparse matrix A (cont.)
    // eq.9b
    //
    // k = 0: path no.0, k = 1: path no.2
    for ( int k = 0; k < 2; ++k )
      {
        // int cen_param_id = (k == 0) ? cs_vertices[1].idx() : alt_param_id[cs_vertices[1]];
        int cen_param_id = (k == 0) ? cs_vertices[0].idx() : cs_vertices[2].idx();
        int l = (k == 0) ? 0 : 2;
        for ( int j = 1; j < path[l].size() - 1; ++j )
          {
            MyMesh::VertexHandle vh = path[l][j];
            int i = vh.idx();
            int rev_i = alt_param_id[vh];

            if ( k == 0 ) // l = 0
              {
                // R(-90)
                // x
                tripletList.push_back( Eigen::Triplet<double>(rev_i, rev_i, 1.0) );
                tripletList.push_back( Eigen::Triplet<double>(rev_i, cen_param_id, -1.0) );
                tripletList.push_back( Eigen::Triplet<double>(rev_i, i+n_param, 1.0) );
                tripletList.push_back( Eigen::Triplet<double>(rev_i, cen_param_id+n_param, -1.0) );
                // y
                tripletList.push_back( Eigen::Triplet<double>(rev_i+n_param, rev_i+n_param, -1.0) );
                tripletList.push_back( Eigen::Triplet<double>(rev_i+n_param, cen_param_id+n_param, 1.0) );
                tripletList.push_back( Eigen::Triplet<double>(rev_i+n_param, i, 1.0) );
                tripletList.push_back( Eigen::Triplet<double>(rev_i+n_param, cen_param_id, -1.0) );
              }
            else // l = 2
              {
                // R(90)
                // x
                tripletList.push_back( Eigen::Triplet<double>(rev_i, rev_i, -1.0) );
                tripletList.push_back( Eigen::Triplet<double>(rev_i, cen_param_id, 1.0) );
                tripletList.push_back( Eigen::Triplet<double>(rev_i, i+n_param, 1.0) );
                tripletList.push_back( Eigen::Triplet<double>(rev_i, cen_param_id+n_param, -1.0) );
                // y
                tripletList.push_back( Eigen::Triplet<double>(rev_i+n_param, rev_i+n_param, 1.0) );
                tripletList.push_back( Eigen::Triplet<double>(rev_i+n_param, cen_param_id+n_param, -1.0) );
                tripletList.push_back( Eigen::Triplet<double>(rev_i+n_param, i, 1.0) );
                tripletList.push_back( Eigen::Triplet<double>(rev_i+n_param, cen_param_id, -1.0) );
              }
          }
      }

    Eigen::SparseMatrix<double> spmat( 2*n_param, 2*n_param );
    spmat.setFromTriplets( tripletList.begin(), tripletList.end() );
    spmat.makeCompressed();

    // cout << "matrix A" << endl;
    // cout << spmat << endl;

    // cout << "before param " << endl;
    // for ( int i = 0; i < paramx.size(); ++i )
    //   {
    //     cout << i << " " << paramx[i] << " " << paramy[i] << endl;
    //   }
    
    // setup vector b
    //Eigen::VectorXd xx(n_param), xy(n_param), bx(n_param), by(n_param);
    Eigen::VectorXd xx(2*n_param), bx(2*n_param);
    for ( int i = 0; i < n_param; ++i )
      {
        // set for x coord
        bx[i] = paramx[i];
        // set for y coord
        bx[i+n_param] = paramy[i];
      }

    // solve x
    if ( sol_ == BICGSTAB )
      {
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver(spmat);
        std::cout << "solve x and y coords ..." << std::endl;
        xx = solver.solve(bx);
        // std::cout << "solve ycoords ..." << std::endl;
        // xy = solver.solve(by);
      }
    else if ( sol_ == SPARSELU )
      {
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;

        solver.analyzePattern(spmat); 
        // Compute the numerical factorization 
        solver.factorize(spmat); 
        //Use the factors to solve the linear system 
        std::cout << "solve x and y coords ..." << std::endl;
        xx = solver.solve(bx);
        // std::cout << "solve ycoords ..." << std::endl;
        // xy = solver.solve(by);
      }

    for ( int i = 0; i < n_param; ++i )
      {
        paramx[i] = xx[i];
        paramy[i] = xx[i+n_param];
      }

#if 0    
    cout << "computed param " << endl;
    for ( int i = 0; i < paramx.size(); ++i )
      {
        cout << i << " " << paramx[i] << " " << paramy[i] << endl;
      }
#endif

    mesh.remove_property(vvws);

    //
    // "cut" mesh along the path
    //

#if 0
    // make consistent path from path no.0 and no.1
    std::vector<MyMesh::VertexHandle> path01;
    for ( int i = 0; i < path.size() / 2; ++i )
      for ( int j = 0; j < path[i].size(); ++j )
        if ( (i != 1) || (j != 0) ) path01.push_back( path[i][j] );
#endif

    // extract right side faces of a path
    std::vector<MyMesh::FaceHandle> right_faces;
    // for OF_BOUNDARY and OF_CONESINGULARITY_PI2 vertices,
    // not for access OF_CONESINGULARITY_PI vertices
    for ( int i = 1; i < path01.size() - 1; ++i )
      {
        MyMesh::VertexHandle vh = path01[i];    // current boundary vertex
        MyMesh::VertexHandle pvh = path01[i-1]; // prev boundary vertex
        MyMesh::VertexHandle nvh = path01[i+1]; // next boundary vertex

        // right faces
        // cout << "vt " << vh.idx() << " cw in halfedge circulator " << endl;
        MyMesh::VertexIHalfedgeCWIter vih_cwit = start_vih_cwiter( vh, pvh, mesh );
        while ( mesh.from_vertex_handle( *vih_cwit ) != nvh )
          {
            MyMesh::VertexHandle ivh = mesh.from_vertex_handle( *vih_cwit );
            MyMesh::FaceHandle fh = mesh.face_handle( *vih_cwit );
            right_faces.push_back( fh );
            // cout << "\tivt " << ivh.idx() << " fc " << fh.idx() << endl;
            ++vih_cwit;
            if ( !(vih_cwit->is_valid()) ) vih_cwit = mesh.vih_cwiter( vh );
          }

        // erase duplicate elements
        std::sort(right_faces.begin(), right_faces.end());
        right_faces.erase(std::unique(right_faces.begin(), right_faces.end()), right_faces.end());

        // left faces
        // cout << "vt " << vh.idx() << " ccw out halfedge circulator " << endl;
        MyMesh::VertexOHalfedgeCCWIter voh_ccwit = start_voh_ccwiter( vh, pvh, mesh );
        while ( mesh.to_vertex_handle( *voh_ccwit ) != nvh )
          {
            MyMesh::VertexHandle ovh = mesh.to_vertex_handle( *voh_ccwit );
            MyMesh::FaceHandle fh = mesh.face_handle( *voh_ccwit );
            // cout << "\tovt " << ovh.idx() << " fc " << fh.idx() << endl;
            ++voh_ccwit;
            if ( !(voh_ccwit->is_valid()) ) voh_ccwit = mesh.voh_ccwiter( vh );
          }
      }

    // create new vertex handles for path vertices with having
    // the same positions to the original vertices
    // OF_BOUNDARY vertices
    for ( int i = 0; i < path.size() / 2; ++i ) // search for only two paths
      {
        for ( int j = 1; j < path[i].size() - 1; ++j )
          {
            MyMesh::VertexHandle vh = path[i][j];    // current boundary vertex
            // new vertex
            MyMesh::VertexHandle new_vh = mesh.add_vertex( mesh.point( vh ) );
            // color
            if ( mesh.has_vertex_colors() )
              {
                MyMesh::Color c = mesh.color(vh);
                mesh.set_color( new_vh, c );
              }
            // cout << "new bd vt " << new_vh.idx() << " param_id " << alt_param_id[vh] << endl;
          }
      }
    // OF_CONESINGULARITY_PI2 vertex
    cs_vt = cs_vertices[1];
    // new vertex
    MyMesh::VertexHandle new_cs_vh = mesh.add_vertex( mesh.point(cs_vt) );
    // color
    if ( mesh.has_vertex_colors() )
      {
        MyMesh::Color cc = mesh.color(cs_vt);
        mesh.set_color( new_cs_vh, cc );
      }
    // cout << "new cs vt " << new_cs_vh.idx() << " param_id " << alt_param_id[cs_vt] << endl;

    // create vertex handles for new faces
    std::vector<std::vector<MyMesh::VertexHandle> > faces_vhandles( right_faces.size() );
    for ( int i = 0; i < right_faces.size(); ++i )
      {
        // cout << "i = " << i << endl;
        for ( MyMesh::FaceVertexIter fv_it = mesh.fv_iter( right_faces[i] );
              fv_it.is_valid(); ++fv_it )
          {
            MyMesh::VertexHandle fvh = *fv_it; // vv_it.handle();
            if ( (vertex_type[fvh.idx()] == OF_BOUNDARY) ||
                 (vertex_type[fvh.idx()] == OF_CONESINGULARITY_PI2) )
              {
                int new_id = alt_param_id[fvh];
                MyMesh::VertexHandle nfvh = mesh.vertex_handle(new_id);
                // cout << "new_id " << new_id << " created " << nfvh.idx() << endl;
                faces_vhandles[i].push_back( nfvh );
              }
            else // internal
              {
                faces_vhandles[i].push_back( fvh );
              }
          }
      }

    // change path vertices to new ones
    for ( int i = path.size() / 2; i < path.size(); ++i ) // search for only two paths
      {
        for ( int j = 0; j < path[i].size(); ++j )
          {
            if ( ((i == path.size() / 2) && (j != 0)) ||
                 ((i == path.size() - 1) && (j != path[i].size()-1)) )
              {
                MyMesh::VertexHandle ovh = path[i][j];
                path[i][j] = mesh.vertex_handle( alt_param_id[ovh] );
              }
          }
      }

    // delete right_faces and add new right faces
    mesh.request_face_status();
    // mesh.request_edge_status();
    // mesh.request_vertex_status();

    for ( int i = 0; i < right_faces.size(); ++i )
      {
        // cout << "delete face " << right_faces[i].idx() << endl;
        mesh.delete_face( right_faces[i], false );
        MyMesh::FaceHandle new_face = mesh.add_face( faces_vhandles[i] );
        // cout << "new face " << new_face.idx() << endl;
      }

    // cut done.
    mesh.garbage_collection();
 
    return 0;
  };

};

#endif // _PARAM_HXX
