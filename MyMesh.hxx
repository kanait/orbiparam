#ifndef _MYMESH_HXX
#define _MYMESH_HXX

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>

struct MyTraits : public OpenMesh::DefaultTraits
{
  typedef OpenMesh::Vec3d Point; // use double-values points
  typedef OpenMesh::Vec3d Normal; // use double-values points
  typedef OpenMesh::Vec2d TexCoord; // use double-values points
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  MyMesh;
//typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

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
