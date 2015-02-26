
#if !defined( __POINT_PROCESS_CORE_MARKED_GRID_HPP__ )
#define __POINT_PROCESS_CORE_MARKED_GRID_HPP__

#include <boost/unordered_map.hpp>
#include <map>
#include <boost/optional.hpp>
#include <math-core/types.hpp>
#include <iostream>
#include <math-core/geom.hpp>


namespace point_process_core {

  using namespace math_core;

  // Description:
  // A cell for a marked grid
  struct marked_grid_cell_t
  {
    size_t n;
    std::vector<long> coordinate;
    marked_grid_cell_t() 
      : n(0),
	coordinate()
    {}
  };

  
  
  // Description:
  // A simple gridded n-dimensional space
  // with a potential mark at each grid cell
  template<typename T_Mark>
  class marked_grid_t
  {
  public:

    //typedef std::map<marked_grid_cell_t,T_Mark> map_t;
    typedef boost::unordered_map<marked_grid_cell_t, T_Mark> map_t;

    // Description:
    // Default constructor which creates an empty invalid grid
    marked_grid_t()
    {}

    // Description:
    // Creates a new marked grid with given resolution and
    // with given window
    marked_grid_t( const math_core::nd_aabox_t& window,
		   const double resolution )
    {
      _init( window, window.start,
	     std::vector<double>( window.start.n, resolution ) );
    }
    
    // Description:
    // Creates a new marked grid with given resolution and window,
    // Where the origin cell is at the given point
    marked_grid_t( const math_core::nd_aabox_t& window,
		   const math_core::nd_point_t& origin,
		   const double resolution )
    {
      _init( window,
	     origin,
	     std::vector<double>( origin.n, resolution ) );
    }
    

    // Description:
    // Creates a new marked grid with given resolution and window,
    // Where the origin cell is at the given point. The resolutions 
    // are variable between doimensions
    marked_grid_t( const math_core::nd_aabox_t& window,
		   const math_core::nd_point_t& origin,
		   const std::vector<double>& resolutions )
      : _bounds( window ),
	_origin( origin ),
	_cell_sizes( resolutions ),
	_map()
    {
      _init( window, origin, resolutions );
    }

    // Description:
    // Destroy a marked grid
    virtual ~marked_grid_t()
    {}


    // Description:
    // Returns a new grid of the given type (whichmight be different!)
    // that has the same structure (window, cellsie, origin) as this
    // grid
    template<typename T_New_Mark>
    marked_grid_t<T_New_Mark> copy_structure() const
    {
      return marked_grid_t<T_New_Mark>( _bounds,
					_origin,
					_cell_sizes );
    }


    // Description:
    // Access a given grid cell
    boost::optional<T_Mark> operator() ( const marked_grid_cell_t& cell ) const
    {
      typename map_t::const_iterator fiter
	= _map.find( cell );
      
      // debug
      //std::cout << "    mk() called" << std::endl;

      if( fiter != _map.end() ) {
	return boost::optional<T_Mark>( fiter->second );
      }
      return boost::optional<T_Mark>();
    }
    boost::optional<T_Mark> operator() ( const math_core::nd_point_t& point ) const
    {
      return this->operator()( this->cell( point) );
    }

    // // Description:
    // // Accessor with default if no mark
    // T_Mark operator()( const marked_grid_cell_t& cell,
    // 		       const T_Mark& default_mark ) 
    // {
    //   std::cout << "    mk(def) called" << std::endl;
    //   boost::optional<T_Mark> mark = this->operator()( cell );
    //   if( mark )
    // 	return *mark;
    //   std::cout << "    mk(def) using def" << std::endl;
    //   this->set( cell, default_mark );
    //   return default_mark;
    // }
    // T_Mark operator() (const math_core::nd_point_t& point,
    // 		       const T_Mark& default_mark )
    // {
    //   return this->operator()( point, default_mark );
    // }

    // Description:
    // Set a particular mark
    void set( const marked_grid_cell_t& cell,
	      const T_Mark& mark )
    {
      _map[ cell ] = mark;
    }
    void set( const math_core::nd_point_t& point,
	      const T_Mark& mark )
    {
      this->set( this->cell( point) , mark );
    }

    // Description:
    // Remove a mark
    void clear_mark( const marked_grid_cell_t& cell ) 
    {
      _map.erase( cell );
    }
    void clear_mark( const math_core::nd_point_t& point )
    {
      this->clear_mark( this->cell( point ) );
    }
    

    // Description:
    // Returns a vector with all possible cells in this grid
    // Note: Be Carefull, this could be HUGE!
    std::vector<marked_grid_cell_t> all_cells() const
    {
      std::vector<marked_grid_cell_t> cells;
      marked_grid_cell_t min_cell = this->cell( _bounds.start );
      marked_grid_cell_t max_cell = this->cell( _bounds.end );
      marked_grid_cell_t c = min_cell;
      while( is_cell_in_lower_iteration( c, max_cell ) ) {
	cells.push_back( c );
	c = cell_next_higher_iteration( c, min_cell, max_cell );
      }
      cells.push_back( max_cell );
      return cells;
    }

    // Description:
    // Returns a vector with all possible cells in a canonical ordering.
    // Note: be careful, this could be HUGE!
    std::vector<marked_grid_cell_t> all_cells_ordered() const
    {
      return all_cells();
    }
 
    // Description:
    // Returns a vector with all marked cells
    std::vector<marked_grid_cell_t> all_marked_cells() const
    {
      std::vector<marked_grid_cell_t> cells;
      typename map_t::const_iterator iter;
      for( iter = _map.begin(); iter != _map.end(); ++iter ) {
	cells.push_back( iter->first );
      }
      return cells;
    }

    // Description:
    // Clear all marks from this grid
    void clear()
    {
      _map.clear();
    }

    // Description:
    // Returns the box which encompases a particular cell in this grid
    math_core::nd_aabox_t region( const marked_grid_cell_t& cell ) const
    {
      assert( cell.n == _origin.n );
      math_core::nd_point_t start, end;
      start.n = cell.n;
      start.coordinate.resize( start.n );
      end.n = cell.n;
      end.coordinate.resize( end.n );
      for( size_t i = 0; i < cell.n; ++i ) {
	start.coordinate[i] = cell.coordinate[i] * _cell_sizes[i] + _origin.coordinate[i];
	end.coordinate[i] = start.coordinate[i] + _cell_sizes[i];
      }
      math_core::nd_aabox_t reg;
      reg.n = cell.n;
      reg.start = start;
      reg.end = end;
      return reg;
    }
    math_core::nd_aabox_t region( const math_core::nd_point_t& point ) const
    {
      return this->region( this->cell( point ) );
    }

    // Description:
    // Returns the cell for a given point
    marked_grid_cell_t cell( const math_core::nd_point_t& point ) const
    {
      assert( _origin.n == point.n );
      marked_grid_cell_t c;
      c.n = _origin.n;
      c.coordinate.resize( c.n );
      
      for( size_t i = 0; i < c.n; ++i ) {
	double x = point.coordinate[i];
	c.coordinate[i] = floor( (x - _origin.coordinate[i]) / _cell_sizes[i] );
      }
      return c;
    }


    // Description:
    // Equality for grids
    bool operator== (const marked_grid_t<T_Mark>& b ) const
    {
      if( ! (
	     _origin.coordinate == b._origin.coordinate &&
	     _cell_sizes == b._cell_sizes &&
	     _bounds.start.coordinate == b._bounds.start.coordinate &&
	     _bounds.end.coordinate == b._bounds.end.coordinate ) )
	return false;
      

      return _map == b._map;
    }


    // Description:
    // Returns the cell sizes
    std::vector<double> cell_sizes() const
    {
      return _cell_sizes;
    }

    // Description:
    // Returns the bounds for this grid
    nd_aabox_t window() const
    {
      return _bounds;
    }
    
  protected:

    
    // Description:
    // The origin for grid cells
    math_core::nd_point_t _origin;
    
    // Description:
    // The size of cell for each dimension
    std::vector<double> _cell_sizes;

    // Description:
    // The bounds of the grid
    math_core::nd_aabox_t _bounds;
    
    // Description:
    // The mapping between cells and marks
    map_t _map;

    // Description:
    // init this objest
    void _init( const math_core::nd_aabox_t& window,
		const math_core::nd_point_t& origin,
		const std::vector<double>& resolutions )
    {
      _bounds = window;
      _origin = origin;
      _cell_sizes = resolutions;
      _map = map_t();
      
      // debug
      //std::cout << "    mk created!" << std::endl;

    }

    // Description:
    // Returns whetehr the cell is is a lower iteration 
    // that another (in terms of iterating over a window of cells)
    bool is_cell_in_lower_iteration( const marked_grid_cell_t& a,
				     const marked_grid_cell_t& b ) const
    {
      assert( a.n == b.n );
      for( size_t i = 0; i < a.n; ++i ) {
	if( a.coordinate[i] < b.coordinate[i] )
	  return true;
      }
      return false;
    }

    // Descriptiopn:
    // Returns the next cell in the iteration 
    // given a point in the iteration and the maximum cell
    marked_grid_cell_t
    cell_next_higher_iteration( const marked_grid_cell_t& c,
				const marked_grid_cell_t& min,
				const marked_grid_cell_t& max ) const
    {
      return cell_next_higher_iteration( c, min, max, c.n - 1 );
    }

    // Descriptiopn:
    // Returns the next cell in the iteration 
    // given a point in the iteration and the maximum cell
    marked_grid_cell_t
    cell_next_higher_iteration( const marked_grid_cell_t& c,
				const marked_grid_cell_t& min,
				const marked_grid_cell_t& max,
				const size_t index ) const
    {
      assert( c.n == max.n );
      marked_grid_cell_t next = c;
      if( next.coordinate[index] + 1 >
	  max.coordinate[index] ) {
	
	// we are done with iteration if index == 0 and we cannot add 1
	if( index == 0 ) {
	  return max;
	}

	// ok, we need to carry over into hte next index
	// so we clamp down on our index variuable and 
	// try to increase the previous index
	next.coordinate[index] = min.coordinate[index];
	return cell_next_higher_iteration( next,
					   min,
					   max,
					   index - 1 );
      } else {
	next.coordinate[index] += 1;
	return next;
      }
    }
				     

  };

  //====================================================================
  

  // Description:
  // Output stream operator
  std::ostream& operator<< (std::ostream& os,
			    const marked_grid_cell_t& cell);

  //====================================================================
  
  // Description:
  // The hash function for cells
  size_t hash_value( const marked_grid_cell_t& cell );

  //====================================================================
  
  // Description:
  // Eqaulity for cells
  bool operator== (const marked_grid_cell_t& a,
		   const marked_grid_cell_t& b );
  bool operator!= (const marked_grid_cell_t& a,
		   const marked_grid_cell_t& b );

  // Description:
  // Less than for cells
  bool operator< (const marked_grid_cell_t& a,
		  const marked_grid_cell_t& b );

  //====================================================================
  

  // Description:
  // Returns the cells which are inside of the given window within
  // the marked grid.
  template<class T>
  std::vector<marked_grid_cell_t> 
  cells_fully_inside_window( const marked_grid_t<T>& grid,
			     const math_core::nd_aabox_t& window )
  {
    std::vector< marked_grid_cell_t > res;
    std::vector< marked_grid_cell_t > all_cells = grid.all_cells();
    std::vector< marked_grid_cell_t >::const_iterator iter;
    for( iter = all_cells.begin();
	 iter != all_cells.end();
	 ++iter ) {
      if( is_fully_inside( grid.region( *iter ), window ) ) {
	res.push_back( *iter );
      }
    }
    return res;
  }

  //====================================================================

  // Description:
  // Returns the smallest window which includes all the given cells
  template<class T>
  math_core::nd_aabox_t
  enclosing_window_for_cells( const marked_grid_t<T>& grid,
			      const std::vector<marked_grid_cell_t>& cells )
  {
    if( cells.empty() )
      return nd_aabox_t();

    // simply find max start and min end points for the cell regions
    nd_point_t max_point = grid.region( cells[0] ).end;
    nd_point_t min_point = grid.region( cells[0] ).start;
    for( size_t i = 1; i < cells.size(); ++i ) {
      nd_point_t s = grid.region( cells[i] ).start;
      nd_point_t e = grid.region( cells[i] ).end;

      // swap points if for some reason end < start
      if( point_lexicographical_compare( e, s ) ) {
	s = e;
	e = grid.region( cells[i] ).start;
      }
      
      // update max/min points
      if( point_lexicographical_compare( max_point, e ) ) {
        max_point = e;
      }
      if( point_lexicographical_compare( s, min_point ) ) {
	min_point = s;
      }
    }

    return aabox( min_point, max_point );
  }


  //====================================================================
  
  // Description:
  // Calculate the 0-1 loss distance between two marked grids
  // *with the same cells*
  template< class T >
  double marked_grid_distance( const marked_grid_t<T>& a,
			       const marked_grid_t<T>& b )
  {
    if( a.cell_sizes() != b.cell_sizes() ) {
      throw std::runtime_error( "cannot calculate distance between different sized markd grids" );
    }
    
    // ok, jsut calculate the 0-1 sum loss ( L1 )
    double dist = 0.0;
    for( auto cell : a.all_cells() ) {
      if( a( cell ) != b( cell ) ) {
	dist += 1.0;
      }
    }
    
    return dist;
  }
  
  //====================================================================
  
  // Description:
  // Computes a boolean marked grid from a point set, where the
  // cell is true if there was at least one point in the point set
  // inisde the cell region
  inline
  marked_grid_t<bool>
  point_set_as_grid( const std::vector<nd_point_t>& points,
		     const nd_aabox_t& window,
		     const double& epsilon )
  {
    marked_grid_t<bool> grid( window, epsilon );
    for( auto p : points ) {
      grid.set( p, true );
    }
    return grid;
  }
  
  //====================================================================

  // Description:
  // Create image from 2D marked grid
  void save_bmp( const std::string& filename, 
		 const marked_grid_t<double>& grid );
		 

  //====================================================================

  
}

#endif

