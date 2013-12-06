
#if !defined( __P2L_POINT_PROCESS_CORE_griding_HPP__ )
#define __P2L_POINT_PROCESS_CORE_griding_HPP__


#include "marked_grid.hpp"
#include <boost/function.hpp>


namespace point_process_core {


  //-------------------------------------------------------------------------
  

  // Description:
  // Given a function and a grid, evaluates the function at
  // *every* grid cell center and stores the result as the mark
  //  on the grid
  template< class T_Mark, class T_Result = T_Mark>
  void evaluate_function_at_grid_center
  ( const boost::function< T_Result( const math_core::nd_point_t& ) >& f,
    marked_grid_t<T_Mark>& grid )
  {
    for( auto cell : grid.all_cells() ) {
      math_core::nd_point_t x = math_core::centroid( grid.region( cell ) );
      T_Mark m = f( x );
      grid.set( cell, m );
    }
  }
  
  //-------------------------------------------------------------------------

  // Description:
  // Plot a marked grid suing linear grid-cells axis
  // and the height (y) as teh mark
  template<class T>
  void plot_linear_grid( const std::string& filename,
			 const marked_grid_t<T>& grid )
  {
    
  }

  //-------------------------------------------------------------------------
  //-------------------------------------------------------------------------
  
}


#endif

