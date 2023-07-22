if(TARGET sample_queries::sample_queries)
  return()
endif()

include(ExternalProject)
include(FetchContent)

set(SAMPLE_QUERIES_ROOT "${PROJECT_SOURCE_DIR}/tests/data/" CACHE PATH "Where should we download sample CCD queries?")

ExternalProject_Add(
  sample_queries_download
  PREFIX "${FETCHCONTENT_BASE_DIR}/sample_queries"
  SOURCE_DIR ${SAMPLE_QUERIES_ROOT}

  GIT_REPOSITORY https://github.com/Continuous-Collision-Detection/Sample-Queries.git
  GIT_TAG 4d6cce33477d8d5c666c31c8ea23e1aea97be371

  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  LOG_DOWNLOAD ON
)

# Create a dummy target for convenience
add_library(sample_queries INTERFACE)
add_library(sample_queries::sample_queries ALIAS sample_queries)

add_dependencies(sample_queries sample_queries_download)

target_compile_definitions(sample_queries INTERFACE SAMPLE_QUERIES_DIR=\"${SAMPLE_QUERIES_ROOT}\")