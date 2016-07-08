#ifndef _MYMESH_HXX
#define _MYMESH_HXX

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;
typedef OpenMesh::Vec3d OMVector3d;

class Elen {
public:
  void setLen( double l ) { len_ = l; };
  double len() const { return len_; };
private:
  double len_;
};

OpenMesh::EPropHandleT<Elen> elen;

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

OpenMesh::VPropHandleT<Fixed> fffs;

#endif // _MYMESH_HXX
