
#include "marked_grid.hpp"
#include <stdexcept>

#define cimg_use_magick
#include <cimg/CImg.h>

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

  void save_bmp( const std::string& filename, 
		 const marked_grid_t<double>& grid )
  {

    // make sure this is 2D
    assert( grid.window().n == 2 );
    if( grid.window().n != 2 ) {
      throw std::runtime_error("cannot save BMP for marked_grid of dimmension other that 2");
    }
    
    // ok, calculate the width,height of image
    int width, height;
    width = (grid.window().end.coordinate[0] - grid.window().start.coordinate[0]) / grid.cell_sizes()[0];
    height = (grid.window().end.coordinate[1] - grid.window().start.coordinate[1]) / grid.cell_sizes()[1];
    
    // create the CImg
    CImg<unsigned char> image( width, height, 1, 3, 0 );

    // get hte max element
    double max_value = -std::numeric_limits<double>::infinity();
    std::vector<marked_grid_cell_t> cells = grid.all_cells();
    for( marked_grid_cell_t cell : cells ) {
      if( grid(cell) &&
	  max_value < *grid(cell) ) {
	max_value = *grid(cell);
      }
      if( grid( cell ) ) {
	//std::cout << "  grid val: " << *grid( cell ) << std::endl;
      }
    }
    std::cout << "saving BMP: (" << width << "x" << height << ") max_value = " << max_value << std::endl;

    
    // ok, now go thorugh every cell and push it's mark onto hte image
    cells = grid.all_cells();
    for( marked_grid_cell_t cell : cells ) {
      assert( cell.n == 2 );
      assert( cell.coordinate.size() == cell.n );
      long x = cell.coordinate[0];
      long y = cell.coordinate[1];
      if( grid( cell ) ) {
	int v = 55.0 + 200.0 * ( (*grid(cell)) / max_value );
	image( x, y, 0, 0 ) = v;
	image( x, y, 0, 1 ) = v;
	image( x, y, 0, 2 ) = v;
	//std::cout << "  (" << x << "," << y << "): " << (int)image(x,y) << "  [" << v << "]" << std::endl;
      }
    }
    
    // save as BMP
    image.save_bmp( filename.c_str() );
  }

}
