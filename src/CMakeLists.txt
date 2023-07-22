set(SOURCES
  autogen.cpp
  autogen.hpp
  config.hpp
  ccd.cpp
  ccd.hpp
  cubic.cpp
  cubic.hpp
  geometry.cpp
  geometry.hpp
  math.cpp
  math.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(far_ccd PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################
