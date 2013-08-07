
#include <point-process-core/context.hpp>
#include <iostream>

using namespace point_process_core;

int main()
{

  std::cout << context_filename( "file.txt" ) << std::endl;

  // create a new scope
  {
    scoped_context_switch context( chain_context( context_t( "test-context" ) ) );
    
    std::cout << context_filename( "file.txt" ) << std::endl;
    
    {
      scoped_context_switch context( chain_context( context_t( "context-b" ) ) );

      std::cout << context_filename( "file.txt" ) << std::endl;
    }

    std::cout << context_filename( "file.txt" ) << std::endl;
			   
  }

  std::cout << context_filename( "file.txt" ) << std::endl;

  return 0;
}
