set(SOURCES
  autogen.cpp
  autogen.hpp
  config.hpp
  ccd.hpp
  ccd.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(far_ccd PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################
