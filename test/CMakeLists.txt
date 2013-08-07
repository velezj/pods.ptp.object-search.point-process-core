
add_executable( test-marked-grid
  test-marked-grid.cpp )
pods_use_pkg_config_packages( test-marked-grid
  math-core lcmtypes_math-core
  point-process-core )
pods_install_executables( test-marked-grid )



# add_executable( point-process-core-test-context
#   test-context.cpp )
# pods_use_pkg_config_packages( point-process-core-test-context
#   math-core lcmtypes_math-core
#   point-process-core )
# pods_install_executables( point-process-core-test-context )

