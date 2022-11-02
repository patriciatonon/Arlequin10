project(runArlequin)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")

find_package(MPI)
find_package(Boost COMPONENTS system filesystem REQUIRED)
find_package(PETSc REQUIRED)

include_directories(include ${MPI_INCLUDE_PATH} ${PETSC_INCLUDES})
add_definitions(${PETSC_DEFINITIONS})

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

include(cmake-modules/EnableEigen.cmake)
enable_eigen()

add_subdirectory(src)

link_directories(${PETSC_DIR}/${PETSC_ARCH}/lib)

add_executable(${PROJECT_NAME} main.cpp ${CXX_SOURCE_FILES})

target_link_libraries(runArlequin metis Eigen3::Eigen ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${PETSC_LIBRARIES})
