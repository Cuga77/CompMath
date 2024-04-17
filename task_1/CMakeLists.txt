cmake_minimum_required(VERSION 3.10)
project(Brusselator)

set(CMAKE_CXX_STANDARD 14)

set(SFML_DIR "C:/vcpkg/installed/x64-windows/share/sfml/lib/cmake/SFML")
find_package(SFML 2.5 COMPONENTS system window graphics CONFIG REQUIRED)

add_executable(Brusselator main.cpp brusselator.cpp)

target_link_libraries(Brusselator)