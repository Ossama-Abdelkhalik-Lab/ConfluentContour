cmake_minimum_required(VERSION 3.20)
project(ConfluentContour VERSION 1.0.0 LANGUAGES Fortran)
enable_testing()
option(CONFLUENTCONTOUR_OPENMP "Parallelize contour_batch chunks with OpenMP" ON)
if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR "ConfluentContour is built and tested with GNU Fortran")
endif()

# The solver selects its sine/cosine entry point with the preprocessor:
# __sincos_stret on Apple platforms, the C sincos elsewhere.
set_source_files_properties(fortran/confluentcontour.f90 PROPERTIES Fortran_PREPROCESS ON)
if(APPLE)
    set_source_files_properties(fortran/confluentcontour.f90 PROPERTIES COMPILE_DEFINITIONS CONTOUR_APPLE_SINCOS)
endif()

add_library(confluentcontour fortran/confluentcontour.f90 fortran/contour_simple.f90 fortran/lambert_api.f90)
set_target_properties(confluentcontour PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
target_include_directories(confluentcontour PUBLIC "${CMAKE_BINARY_DIR}/modules")
target_compile_options(confluentcontour PRIVATE -O3 -ffree-line-length-none -fno-fast-math)
if(CONFLUENTCONTOUR_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
    target_link_libraries(confluentcontour PUBLIC OpenMP::OpenMP_Fortran)
endif()

function(contour_test name)
    add_executable(${name} ${ARGN})
    set_target_properties(${name} PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/${name}_modules")
    target_link_libraries(${name} PRIVATE confluentcontour)
    target_compile_options(${name} PRIVATE -O3 -ffree-line-length-none -fno-fast-math)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

# Native port against analytic truth; deterministic serial/OpenMP output.
contour_test(test_contour_port tests/test_contour_port.f90)
# Upper-boundary coordinates, merged/infeasible roots, parabolic continuity.
contour_test(test_contour_coordinates tests/test_contour_coordinates.f90)
# Table-based multi-revolution feasibility against an independent minimum.
contour_test(test_contour_tables tests/test_contour_tables.f90)
# contour_batch against contour_one, bit for bit.
contour_test(test_contour_batch tests/test_contour_batch.f90)
# ContourSimple against contour_one across regimes.
contour_test(test_contour_simple tests/test_contour_simple.f90)
# Dimensional API: units, statuses, transfer plane.
contour_test(test_lambert_api tests/test_lambert_api.f90)
# Near-antipodal and near-coincident geometry by two-body propagation.
contour_test(test_near_antipodal_api tests/two_body_reference.f90 tests/test_near_antipodal_api.f90)
contour_test(test_near_coincident_api tests/two_body_reference.f90 tests/test_near_coincident_api.f90)
if(CONFLUENTCONTOUR_OPENMP)
    set_tests_properties(test_contour_port test_contour_batch PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=4")
endif()
