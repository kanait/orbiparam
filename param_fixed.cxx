////////////////////////////////////////////////////////////////////
//
//
// parameterization with fixed boundary
//
// $Id: $
//
// Copyright (c) 2013-2015 by Takashi Kanai. All rights reserved.
//
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <string>
// --------------------
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;
typedef OpenMesh::Vec3d Vector3d;

#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>
typedef Eigen::Triplet<double> T;

//#include "spd.hxx"

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

// 0: vertex i, 1-N: its neighbor vertices j
class VVWeights {

public:

  VVWeights(){};
  ~VVWeights(){ w_.clear(); };

  double w( int i ) { return w_[i]; };
  int size() const { return w_.size(); };
  std::vector<double>& ws() { return w_; };
  void addWeight( double d ) { w_.push_back(d); };
  bool setWeight( int i, double d ) {
    if ( i >= (int) w_.size() ) return false;
    w_[i] = d;
    return true;
  };

private:

  std::vector<double> w_;

};

double cotw( MyMesh::Point& p0, MyMesh::Point& p1, MyMesh::Point& p2 )
{
  Vector3d c1( p1 - p0 );
  Vector3d c2( p2 - p0 );
  double sin = OpenMesh::cross( c1, c2 ).length();
  double angle = std::fabs( std::atan2(sin, OpenMesh::dot(c1, c2)) );

  // std::cout << "\tcot " << 1.0 / std::tan( angle ) << " tan " << std::tan( M_PI_2 - angle ) << std::endl;
  // if ( 1.0 / std::tan( angle ) - std::tan( M_PI_2 - angle ) > 1.0e-05 )
  //   std::cout << "diff!" << std::endl;
  return std::tan( M_PI_2 - angle );
}

double tan2w( MyMesh::Point& p0, MyMesh::Point& p1, MyMesh::Point& p2 )
{
  Vector3d c1( p1 - p0 );
  Vector3d c2( p2 - p0 );
  double sin = OpenMesh::cross( c1, c2 ).length();
  double angle = std::fabs( std::atan2(sin, OpenMesh::dot(c1, c2)) );

  // std::cout << "\tcot " << 1.0 / std::tan( angle ) << " tan " << std::tan( M_PI_2 - angle ) << std::endl;
  // if ( 1.0 / std::tan( angle ) - std::tan( M_PI_2 - angle ) > 1.0e-05 )
  //   std::cout << "diff!" << std::endl;
  return std::tan( angle / 2.0 );
}

double computeCotangentWeight( MyMesh& mesh, MyMesh::VertexIHalfedgeIter& vih_it )
{
  MyMesh::HalfedgeHandle iheh = *vih_it; // vih_it.handle();
  // if ( mesh.from_vertex_handle( iheh ) == vh ) std::cout << "i from" << std::endl;
  // else if ( mesh.to_vertex_handle( iheh ) == vh ) std::cout << "i to" << std::endl;
  MyMesh::HalfedgeHandle next_iheh = mesh.next_halfedge_handle( *vih_it ); 
  // if ( mesh.from_vertex_handle( next_iheh ) == vh ) std::cout << "next i from" << std::endl;
  // else if ( mesh.to_vertex_handle( next_iheh ) == vh ) std::cout << "next i to" << std::endl

  MyMesh::Point p0, p1, p2;
  p0 = mesh.point( mesh.to_vertex_handle( next_iheh ) );
  p1 = mesh.point( mesh.from_vertex_handle( iheh ) );
  p2 = mesh.point( mesh.to_vertex_handle( iheh ) );
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

  return (al + be) / 2.0;
}

double computeMeanValueWeight( MyMesh& mesh, MyMesh::VertexIHalfedgeIter& vih_it )
{
  MyMesh::Point p0, p1, p2;

  // gamma
  MyMesh::HalfedgeHandle iheh = *vih_it; // vih_it.handle();
  MyMesh::HalfedgeHandle next_iheh = mesh.next_halfedge_handle( iheh );
  p0 = mesh.point( mesh.from_vertex_handle( next_iheh ) );
  p1 = mesh.point( mesh.to_vertex_handle( next_iheh ) );
  p2 = mesh.point( mesh.from_vertex_handle( iheh ) );
  double ga = tan2w( p0, p1, p2 );

  // delta
  MyMesh::HalfedgeHandle mheh = mesh.opposite_halfedge_handle( iheh );
  MyMesh::HalfedgeHandle next_mheh = mesh.next_halfedge_handle( mheh );
  p0 = mesh.point( mesh.from_vertex_handle( mheh ) ); // i
  p1 = mesh.point( mesh.from_vertex_handle( next_mheh ) ); //j
  p2 = mesh.point( mesh.to_vertex_handle( next_mheh ) );
  double de = tan2w( p0, p1, p2 );

  return (ga + de) / (p0 - p1).length();
}

void constructBoundary( MyMesh& mesh, std::vector<MyMesh::VertexHandle>& bverts )
{
  // find a boundary vertex
  MyMesh::VertexHandle vh0;
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

  // construct boundary vertices
  MyMesh::VertexHandle vh = vh0, prev_vh;
  do {
    bverts.push_back( vh );
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
  std::cout << "#boundary vertices: " << bverts.size() << std::endl;
}

void computeBoundaryMapping( MyMesh& mesh, std::vector<MyMesh::VertexHandle>& bverts,
                             std::vector<double>& paramx,  std::vector<double>& paramy )
{
  // corner coordinates
  std::vector<double> cornerx(RECTANGLE);
  std::vector<double> cornery(RECTANGLE);
  std::vector<MyMesh::VertexHandle> cornerv;
  cornerx[0] = 0.0; cornerx[1] = 1.0; cornerx[2] = 1.0; cornerx[3] = 0.0;
  cornery[0] = 0.0; cornery[1] = 0.0; cornery[2] = 1.0; cornery[3] = 1.0;

  // compute boundary length
  double blength = 0.0;
  int bvn = bverts.size();
  for ( int i = 0; i < bvn-1; ++i )
    {
      blength += (mesh.point( bverts[i] ) - mesh.point( bverts[i+1] )).length();
    }
  blength += (mesh.point( bverts[bvn-1] ) - mesh.point( bverts[0] )).length();
  std::cout << "boundary length " << blength << std::endl;

  // determine corner vertices
  double l4 = blength / (double) RECTANGLE;
  double tl = 0.0;
  cornerv.push_back( bverts[0] );
  paramx[ bverts[0].idx() ] = cornerx[0];
  paramy[ bverts[0].idx() ] = cornery[0];

  int count = 1;
  for ( int i = 0; i < bvn-1; ++i )
    {
      tl += (mesh.point( bverts[i] ) - mesh.point( bverts[i+1] )).length();
      if ( (tl > l4) && (count < RECTANGLE) )
        {
          cornerv.push_back( bverts[i] );
          paramx[ bverts[i].idx() ] = cornerx[count];
          paramy[ bverts[i].idx() ] = cornery[count];
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
              pll += (mesh.point( innerv[j] ) - mesh.point( bverts[j+1] )).length();
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
              cll += (mesh.point( innerv[j] ) - mesh.point( bverts[j+1] )).length();
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
      pll += (mesh.point( innerv[j] ) - mesh.point( bverts[j+1] )).length();
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
      cll += (mesh.point( innerv[j] ) - mesh.point( bverts[j+1] )).length();
    }
}

int main(int argc, char **argv)
{
  MyMesh  mesh;
  // check command line options
  if (argc != 6)
    {
      std::cerr << "Usage:  " << argv[0] << " infile outfile solver weight save" << std::endl;
      std::cerr << "\tsolver ... 1: SparseLU, 2: BICGSTAB" << std::endl;
      std::cerr << "\tweight ... 1: Cotangent Weight, 2: Mean Value Weight" << std::endl;
      std::cerr << "\tsave   ... 1: as 2D Vertex Coordinates, 2: as 2D Texture Coordinates" << std::endl;
      return 1;
    }

  sol = atoi(argv[3]);
  if ( sol == SPARSELU ) std::cout << "solver: SparseLU" << std::endl;
  else if ( sol == BICGSTAB) std::cout << "solver: BiCGSTAB" << std::endl;
  wei = atoi(argv[4]);
  if ( wei == COTW ) std::cout << "weights: Cotangent" << std::endl;
  else if ( wei == MVW ) std::cout << "weights: Mean Value" << std::endl;
  sav = atoi(argv[5]);
  if ( sav == SAVE_VERTEX ) std::cout << "save: as 2D vertex coords" << std::endl;
  else if ( sav == SAVE_TEXCOORD ) std::cout << "save: as 2D texture coords" << std::endl;

  // read mesh from stdin
  if ( ! OpenMesh::IO::read_mesh(mesh, argv[1]) )
    {
      std::cerr << "Error: Cannot read mesh from " << argv[1] << std::endl;
      return 1;
    }

  // initialize 2D parameters
  std::cout << "#vertices: " << mesh.n_vertices() << std::endl;
  std::cout << "#faces: " << mesh.n_faces() << std::endl;
  int n_vt = mesh.n_vertices();
  std::vector<double> paramx( n_vt );
  std::vector<double> paramy( n_vt );
  for ( int i = 0; i < n_vt; ++i )
    {
      paramx[i] = paramy[i] = 0.0;
    }

  //
  // boundary mapping
  //
  // construct boundary vertices
  std::vector<MyMesh::VertexHandle> bverts;
  constructBoundary( mesh, bverts );

  // compute 2D parameters for boundary
  computeBoundaryMapping( mesh, bverts, paramx, paramy );

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

      MyMesh::VertexIHalfedgeIter vih_it;
      for (vih_it=mesh.vih_iter( *v_it ); vih_it.is_valid(); ++vih_it)
        {
          double wdc;
          if ( wei == COTW )
            wdc = computeCotangentWeight( mesh, vih_it );
          else if ( wei == MVW )
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
  if ( sol == BICGSTAB )
    {
      Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver(spmat);
      std::cout << "solve xcoords ..." << std::endl;
      xx = solver.solve(bx);
      std::cout << "solve ycoords ..." << std::endl;
      xy = solver.solve(by);
    }
  else if ( sol == SPARSELU )
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
  std::cout << "done." << std::endl;

  if ( sav == SAVE_VERTEX )
    {
      // set point to 2D parameters
      for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
        {
          MyMesh::VertexHandle vh = *v_it;
          int id = vh.idx();
          MyMesh::Point p( paramx[id], paramy[id], 0.0 );
          mesh.set_point( vh, p );
        }
    }
  else if ( sav == SAVE_TEXCOORD )
    {
      mesh.request_vertex_texcoords2D();
      for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
        {
          MyMesh::VertexHandle vh = *v_it;
          int id = vh.idx();
          MyMesh::TexCoord2D t( paramx[id], paramy[id] );
          mesh.set_texcoord2D( vh, t );
        }
    }

  // write mesh to stdout
  OpenMesh::IO::Options opt;
  if ( mesh.has_vertex_texcoords2D() )
    {
      opt += OpenMesh::IO::Options::VertexTexCoord;
    }
  if ( ! OpenMesh::IO::write_mesh(mesh, argv[2], opt) )
    {
      std::cerr << "Error: cannot write mesh to " << argv[2] << std::endl;
      return 1;
    }

  return 0;
}
