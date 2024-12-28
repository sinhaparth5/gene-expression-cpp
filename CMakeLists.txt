cmake_minimum_required(VERSION 3.10)
project(GeneExpressionAnalysis)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen package
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Define source files directly
add_executable(gene_analysis 
    src/main.cpp
    src/GeneExpressionAnalysis.cpp
    src/ClusteringMethods.cpp
)

# Include directories
target_include_directories(gene_analysis PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${EIGEN3_INCLUDE_DIR}
)

# Link libraries
target_link_libraries(gene_analysis Eigen3::Eigen)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Install targets
install(TARGETS gene_analysis DESTINATION bin)
install(FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/include/GeneExpressionAnalysis.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/DataStructures.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ClusteringMethods.hpp
    DESTINATION include/gene_analysis
)