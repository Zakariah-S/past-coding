cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(Practise)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(prac prac.cc ${sources} ${headers})

target_link_libraries(prac ${Geant4_LIBRARIES})

add_custom_target(Practise DEPENDS prac)