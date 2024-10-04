cmake_minimum_required(VERSION 3.10.0)
project(geom-virt)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(third-party/pmp-library)

add_executable(${PROJECT_NAME} src/main.cpp)

configure_file(${CMAKE_SOURCE_DIR}/inc/cwd.h.in ${CMAKE_SOURCE_DIR}/inc/cwd.h)

target_link_libraries(${PROJECT_NAME} pmp pmp_vis)
target_include_directories(${PROJECT_NAME} PUBLIC inc/)

