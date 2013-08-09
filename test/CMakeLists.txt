
add_executable( p2l-point-process-core-test-marked-grid
  test-marked-grid.cpp )
pods_use_pkg_config_packages( p2l-point-process-core-test-marked-grid
  p2l-math-core lcmtypes_p2l-math-core
  p2l-point-process-core
  p2l-probability-core
  lcmtypes_p2l-probability-core )
pods_install_executables( p2l-point-process-core-test-marked-grid )



# add_executable( point-process-core-test-context
#   test-context.cpp )
# pods_use_pkg_config_packages( point-process-core-test-context
#   math-core lcmtypes_math-core
#   point-process-core )
# pods_install_executables( point-process-core-test-context )

