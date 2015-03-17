
add_executable( object-search.point-process-core-test-marked-grid
  test-marked-grid.cpp )
pods_use_pkg_config_packages( object-search.point-process-core-test-marked-grid
  object-search.math-core 
  object-search.point-process-core
  object-search.probability-core
  )
pods_install_executables( object-search.point-process-core-test-marked-grid )



# add_executable( point-process-core-test-context
#   test-context.cpp )
# pods_use_pkg_config_packages( point-process-core-test-context
#   math-core lcmtypes_math-core
#   point-process-core )
# pods_install_executables( point-process-core-test-context )


add_executable( object-search.point-process-core-test-histogram
  test-histogram.cpp )
pods_use_pkg_config_packages( object-search.point-process-core-test-histogram
  boost-1.54.0
  object-search.math-core 
  object-search.point-process-core
  object-search.probability-core
  )
pods_install_executables( object-search.point-process-core-test-histogram )
