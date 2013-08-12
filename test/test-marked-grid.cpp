
#include <point-process-core/marked_grid.hpp>
#include <math-core/geom.hpp>
#include <math-core/io.hpp>
#include <iostream>


using namespace math_core;
using namespace point_process_core;

int main( int argc, char** argv )
{

  // create a 2d grid
  nd_aabox_t window;
  window.n = 2;
  window.start = point( 0.0,0.0 );
  window.end = point( 10.0,10.0 );
  marked_grid_t<int> grid( window, 1.0 );
  marked_grid_t<int> grid1 ( window, 1.0 );

  // set some marks
  grid.set( point( 3.5, 3.5), 1 );
  grid.set( point( 5.5, 3.5), -3 );
  grid1.set( point( 5.5, 3.5), -3 );
  grid1.set( point( 3.5, 3.5), 1 );
  
  
  // get marked cells
  std::vector<marked_grid_cell_t> marked_cells = grid.all_marked_cells();
  
  for( size_t i = 0; i < marked_cells.size(); ++i ) {
    std::cout << "cell " << marked_cells[i] << std::endl;
  }

  std::cout << "Grids equal: " << (grid == grid1) << std::endl;

  return 0;
}
