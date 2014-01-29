
#include "marked_grid.hpp"
#include <cimg/CImg.h>
#include <stdexcept>

using namespace cimg_library;


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


  //====================================================================

  void save_png( const std::string& filename, 
		 const marked_grid_t<double>& grid )
  {

    // make sure this is 2D
    assert( grid.window().n == 2 );
    if( grid.window().n != 2 ) {
      throw std::runtime_error("cannot save PNG for marked_grid of dimmension other that 2");
    }
    
    // ok, calculate the width,height of image
    int width, height;
    width = (grid.window().end.coordinate[0] - grid.window().start.coordinate[0]) / grid.cell_sizes()[0];
    height = (grid.window().end.coordinate[1] - grid.window().start.coordinate[1]) / grid.cell_sizes()[1];
    
    // create the CImg
    CImg<double> image( width, height, 1, 1, 0 );

    // get hte max element
    double max_value = -std::numeric_limits<double>::infinity();
    for( auto cell : grid.all_cells() ) {
      if( grid(cell) &&
	  max_value < *grid(cell) ) {
	max_value = *grid(cell);
      }
    }
    
    // ok, now go thorugh every cell and push it's mark onto hte image
    for( auto cell : grid.all_cells() ) {
      long x = cell.coordinate[0];
      long y = cell.coordinate[1];
      if( grid( cell ) ) {
	image( x, y ) = ( *grid(cell) / max_value ) * 255;
      }
    }
    
    // save as PNG
    image.save_png( filename.c_str() );
  }

}
