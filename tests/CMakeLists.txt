################################################################################
# Tests
################################################################################

add_executable(far_ccd_tests
  test_edge_edge_ccd.cpp
  test_point_triangle_ccd.cpp
)

target_include_directories(far_ccd_tests PUBLIC ".")

################################################################################
# Required Libraries
################################################################################

target_link_libraries(far_ccd_tests PUBLIC far_ccd::far_ccd)

include(catch2)
target_link_libraries(far_ccd_tests PUBLIC Catch2::Catch2WithMain)

# include(finite_diff)
# target_link_libraries(far_ccd_tests PUBLIC finitediff::finitediff)

# include(json)
# target_link_libraries(far_ccd_tests PUBLIC nlohmann_json::nlohmann_json)

# Extra warnings (link last for highest priority)
include(far_ccd_warnings)
target_link_libraries(far_ccd_tests PRIVATE far_ccd::warnings)

################################################################################
# Compiler options
################################################################################

target_compile_definitions(far_ccd_tests PUBLIC CATCH_CONFIG_ENABLE_BENCHMARKING)

# Test data directory definition
# if(WIN32)
#   target_compile_definitions(far_ccd_tests PUBLIC TEST_DATA_DIR_CSTR="${CMAKE_CURRENT_SOURCE_DIR}\\\\data\\\\")
# else()
#   target_compile_definitions(far_ccd_tests PUBLIC TEST_DATA_DIR_CSTR="${CMAKE_CURRENT_SOURCE_DIR}/data/")
# endif()

################################################################################
# Register tests
################################################################################

FetchContent_GetProperties(catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
catch_discover_tests(far_ccd_tests)
