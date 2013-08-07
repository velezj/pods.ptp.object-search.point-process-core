
#include "marked_grid.hpp"


namespace point_process_core {



  
  // Description:
  // Output stream operator
  std::ostream& operator<< (std::ostream& os,
			    const marked_grid_cell_t& cell)
  {
    os << "Cell{";
    for( size_t i = 0; i < cell.n; ++i ) {
      os << cell.coordinate[i];
      if( i < cell.n - 1 ) {
	os << ",";
      }
    }
    os << "}";
    return os;
  }


  // Description:
  // The hash function for cells
  size_t hash_value( const marked_grid_cell_t& cell ) {
    // std::cout << "hash(" << cell.n << ", ";
    // for( size_t i = 0; i < cell.n; ++i ) {
    //   std::cout << cell.coordinate[i] << " ";
    // }
    // std::cout << ")" << std::endl;
    boost::hash<std::vector<long> > hasher;
    return hasher( cell.coordinate );
  }

  // Description:
  // Eqaulity for cells
  bool operator== (const marked_grid_cell_t& a,
		   const marked_grid_cell_t& b )
  {
    return a.coordinate == b.coordinate;
  }
  bool operator!= (const marked_grid_cell_t& a,
		   const marked_grid_cell_t& b )
  {
    return !(a == b);
  }

  // Description:
  // Less than for cells
  bool operator< (const marked_grid_cell_t& a,
		  const marked_grid_cell_t& b )
  {
    return a.coordinate < b.coordinate;
  }



}
