
#if !defined( __POINT_PROCESS_CORE_CONTEXT_HPP__ )
#define __POINT_PROCESS_CORE_CONTEXT_HPP__


#include <string>
#include <fstream>
#include <boost/optional.hpp>

namespace point_process_core {


  // Description:
  // An identifier for a context
  struct context_t
  {
    std::string id;  
    context_t( const std::string& id )
      : id(id)
    {}
  };

  // Description:
  // "Chain" the current context by adding thes given one at end and
  // return this new chained context
  context_t chain_context( const context_t& new_context );

  // Description:
  // Switches the current context to the given one.
  // Returns the previuous context_t object, if any
  boost::optional<context_t> set_context( const boost::optional<context_t>& c );

  // Description:
  // Returns the current context, if any
  boost::optional<context_t> get_current_context();
  
  // Description:
  // Returns a new filename prepended with the current context
  // and with the given name
  std::string context_filename( const std::string& name );

  
  // Description:
  // A scoped context switch.
  // When created, switches to the given context,
  // then switches to the previous context that was current when
  // the object was created at scope exit of the object
  class scoped_context_switch
  {
  public:
    boost::optional<context_t> _original_context;
    scoped_context_switch( const boost::optional<context_t>& new_context )
    {
      _original_context = set_context( new_context );
    }
    virtual ~scoped_context_switch()
    {
      set_context( _original_context );
    }
  };

}

#endif

