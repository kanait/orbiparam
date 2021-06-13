////////////////////////////////////////////////////////////////////
//
// $Id: VertexPQ.hxx 2021/06/13 22:14:01 kanai Exp $
//
// Copyright (c) 2021 Takashi Kanai
// Released under the MIT license
//
////////////////////////////////////////////////////////////////////

#ifndef _VERTEXPQ_HXX
#define _VERTEXPQ_HXX 1

#include <queue>
using namespace std;

#include "MyMesh.hxx"

struct compare {
  bool operator()(const std::pair<MyMesh::VertexHandle, double>& x,
                  const std::pair<MyMesh::VertexHandle, double>& y) const {
    return x.second > y.second;
  };
};

typedef std::pair<MyMesh::VertexHandle, double> value_type;
typedef std::vector<value_type> container_type;

class VertexPQ {

public:

  bool empty() { return pq_.empty(); };
  int size() const { return pq_.size(); };

  void push( double& err, MyMesh::VertexHandle& vh ) {
    pq_.push( std::make_pair( vh, err ) );
  };

  void pop() { pq_.pop(); };

  double top_err() const { return pq_.top().second; };
  MyMesh::VertexHandle top_vh() const { return pq_.top().first; };

private:

  std::priority_queue<value_type, container_type, compare> pq_;

};

#endif // _VERTEXPQ_HXX
