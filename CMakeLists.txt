cmake_minimum_required (VERSION 3.10)

# boltzi uses newer additions to the language like unordered_map::contains.
set(CMAKE_CXX_STANDARD 20)

set(SRC_PATH ${CMAKE_HOME_DIRECTORY}/src/)
set(CMAKE_CXX_FLAGS "-DTETLIBRARY -O3 -Wall -Wextra -Wpedantic -Wno-unused-variable")
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})

# uncomment for debugging
#set(CMAKE_CXX_FLAGS "-DTETLIBRARY -g -Og -Wall -Wextra -Wpedantic -Wno-unused-variable")
#set(CMAKE_BUILD_TYPE Debug)



project(boltzi)

# To print out timings in some time critical operations.
#add_definitions(-DBOLTZI_MEASURE_TIMINGS)
# To add a dashboard showing graphs of the phonon population and temperature.
#add_definitions(-DBOLTZI_DASHBOARD)
# Use incorrect way of calculating scatter probabilities in order to eliminate
# temperature drifts associated with doing it the correct way.
add_definitions(-DBOLTZI_ALTERNATIVE_SCATTER_PROBABILITY)
# Always randomize the directions of phonos after three-phonon-processes (even for umklapp processes).
#add_definitions(-DBOLTZI_ALWAYS_RANDOMIZE_THREE_PHONON_SCATTERING)
# Allow transitions between branches on scatter events.
#add_definitions(-DBOLTZI_ALLOW_TRANSITIONS)
# Save the trajectory of the first phonon in the linear algorithm in a file called trajectory.csv
#add_definitions(-DBOLTZI_FIRST_TRAJECTORY)
# Save phonon spawn times
#add_definitions(-DBOLTZI_PHONON_SPAWN_TIMES)
# Save Wavevectors of phonons transitioning boundaries
#add_definitions(-DBOLTZI_BOUNDARY_DIRS)
# Save collision points of phonons, inelastically scattered on sample surfaces
#add_definitions(-DBOLTZI_COLLISION_POINTS)



add_executable(boltzi ${SRC_PATH}/main.cpp
    ${SRC_PATH}/sample.cpp ${SRC_PATH}/boundary.cpp ${SRC_PATH}/simulation_controller.cpp ${SRC_PATH}/sample_geometry_factory.cpp
    ${SRC_PATH}/random_init.cpp ${SRC_PATH}/boltzi_config.cpp ${SRC_PATH}/pool.cpp ${SRC_PATH}/silencer.cpp)

target_include_directories(boltzi PRIVATE ${CMAKE_HOME_DIRECTORY}/src/headers/)
target_link_libraries(boltzi tet)
target_link_libraries(boltzi armadillo)
#target_link_libraries(boltzi tablengraph)
target_link_libraries(boltzi yaml-cpp)
