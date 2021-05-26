#ifndef CGAL_DELAUNAY_VERTEX_BASE_WITH_INFO_2_H
#define CGAL_DELAUNAY_VERTEX_BASE_WITH_INFO_2_H

#include <CGAL/license/Mesh_2.h>


#include <CGAL/Triangulation_vertex_base_2.h>

namespace CGAL {

template <typename Info_, class Gt,
          class Vb = Triangulation_vertex_base_2<Gt> >

class Delaunay_mesh_vertex_base_with_info_2 : public Vb
{
  Info_ _info;
	  
public:
  typedef typename Gt::FT          FT;
  typedef typename Vb::Point       Point;
  typedef typename Vb::Face_handle Face_handle;
  typedef Info_                                      Info;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
    typedef Delaunay_mesh_vertex_base_with_info_2<Info, Gt,Vb2> Other;
  };

protected:
  FT sizing_info_;

public:
  Delaunay_mesh_vertex_base_with_info_2()
    : Vb()
    , sizing_info_(0.)
  {}

  Delaunay_mesh_vertex_base_with_info_2(Point p)
    : Vb(p)
    , sizing_info_(0.)
  {}

  Delaunay_mesh_vertex_base_with_info_2(Point p,
                              Face_handle f)
    : Vb(p, f)
    , sizing_info_(0.)
  {}

  void set_sizing_info(const FT& s)
  { 
    sizing_info_ = s;
  }
  const FT& sizing_info() const { return sizing_info_; }
  
  const Info& info() const { return _info; }
  Info&       info()       { return _info; }
};

} // namespace CGAL

#endif
