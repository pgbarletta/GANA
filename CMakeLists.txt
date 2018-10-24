cmake_minimum_required(VERSION 3.9.0)
project(GANA VERSION 0.1.0 LANGUAGES CXX CUDA)

if(MSVC)
    add_compile_options(/W3 /WX)
else()
    add_compile_options(-W -Wakk)
endif()
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")

find_package(CGAL)
find_package(chemfiles)
find_package(NetCDF)

add_executable(GANA
    src/main.cu
    src/continuous.cpp
    src/grid.cpp
    src/utils.cpp)
target_include_directories(GANA
    PRIVATE
        include)

target_link_libraries(GANA
    PRIVATE
        chemfiles
        netcdf
        CGAL
        mpfr
        gmp
        lm
    INTERFACE
        fmt::fmt-header-only)

