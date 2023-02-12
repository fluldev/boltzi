cmake_minimum_required (VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)
set(SRC_PATH ${CMAKE_HOME_DIRECTORY}/../src)
#set(CMAKE_CXX_FLAGS "-DTETLIBRARY -O3 -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS "-DTETLIBRARY -Wall -Wextra -Wpedantic -O3")

project(tests)

#add_library(boltzi STATIC ${SRC_PATH}/sample.cpp ${SRC_PATH}/boundary.cpp ${SRC_PATH}/simulation_controller.cpp ${SRC_PATH}/sample_geometry_factory.cpp)
#target_include_directories(boltzi PRIVATE ${CMAKE_HOME_DIRECTORY}/../src/headers/)
#target_link_libraries(boltzi armadillo)
#target_link_libraries(boltzi tet)

#add_executable(test_grid test_grid.cpp)
#target_include_directories(test_grid PRIVATE ${CMAKE_HOME_DIRECTORY}/../src/headers/)
#
#add_executable(test_random test_random.cpp ${SRC_PATH}/random_init.cpp ${SRC_PATH}/boltzi_config.cpp)
#target_include_directories(test_random PRIVATE ${CMAKE_HOME_DIRECTORY}/../src/headers/)
#target_link_libraries(test_random armadillo)
#
#add_executable(test_phonon test_phonon.cpp ${SRC_PATH}/random_init.cpp ${SRC_PATH}/boltzi_config.cpp)
#target_include_directories(test_phonon PRIVATE ${CMAKE_HOME_DIRECTORY}/../src/headers/)
#target_link_libraries(test_phonon boltzi)
#target_link_libraries(test_phonon armadillo)
#
#add_executable(test_boundary test_boundary.cpp ${SRC_PATH}/random_init.cpp ${SRC_PATH}/boltzi_config.cpp)
#target_include_directories(test_boundary PRIVATE ${CMAKE_HOME_DIRECTORY}/../src/headers/)
#target_link_libraries(test_boundary boltzi)
#target_link_libraries(test_boundary tet)
#target_link_libraries(test_boundary armadillo)
#
#add_executable(test_spawn test_spawn.cpp ${SRC_PATH}/random_init.cpp ${SRC_PATH}/boltzi_config.cpp)
#target_include_directories(test_spawn PRIVATE ${CMAKE_HOME_DIRECTORY}/../src/headers/)
#target_link_libraries(test_spawn boltzi)
#target_link_libraries(test_spawn tet)
#target_link_libraries(test_spawn armadillo)
#
#
#add_executable(test_lin_disp test_lin_disp.cpp ${SRC_PATH}/random_init.cpp ${SRC_PATH}/boltzi_config.cpp)
#target_include_directories(test_lin_disp PRIVATE ${CMAKE_HOME_DIRECTORY}/../src/headers/)
#target_link_libraries(test_lin_disp boltzi)
#target_link_libraries(test_lin_disp tet)
#target_link_libraries(test_lin_disp armadillo)
#
#add_executable(test_lin_material test_lin_material.cpp ${SRC_PATH}/random_init.cpp ${SRC_PATH}/boltzi_config.cpp ${SRC_PATH}/pool.cpp)
#target_include_directories(test_lin_material PRIVATE ${CMAKE_HOME_DIRECTORY}/../src/headers/)
#target_link_libraries(test_lin_material boltzi)
#target_link_libraries(test_lin_material tet)
#target_link_libraries(test_lin_material armadillo)
#
#add_executable(test_advection test_advection.cpp ${SRC_PATH}/random_init.cpp ${SRC_PATH}/boltzi_config.cpp ${SRC_PATH}/pool.cpp)
#target_include_directories(test_advection PRIVATE ${CMAKE_HOME_DIRECTORY}/../src/headers/)
#target_link_libraries(test_advection boltzi)
#target_link_libraries(test_advection tet)
#target_link_libraries(test_advection armadillo)
#
#add_executable(test_dsmc test_dsmc.cpp ${SRC_PATH}/random_init.cpp ${SRC_PATH}/boltzi_config.cpp ${SRC_PATH}/pool.cpp)
#target_include_directories(test_dsmc PRIVATE ${CMAKE_HOME_DIRECTORY}/../src/headers/)
#target_link_libraries(test_dsmc boltzi)
#target_link_libraries(test_dsmc tet)
#target_link_libraries(test_dsmc armadillo)
#
#add_executable(test_pool test_pool.cpp)
#target_include_directories(test_pool PRIVATE ${CMAKE_HOME_DIRECTORY}/../src/headers/)

#add_executable(test_dsmc_spatial test_dsmc_spatial.cpp ${SRC_PATH}/random_init.cpp ${SRC_PATH}/boltzi_config.cpp ${SRC_PATH}/pool.cpp)
#target_include_directories(test_dsmc_spatial PRIVATE ${CMAKE_HOME_DIRECTORY}/../src/headers/)
#target_link_libraries(test_dsmc_spatial boltzi)
#target_link_libraries(test_dsmc_spatial tet)
#target_link_libraries(test_dsmc_spatial armadillo)
#target_link_libraries(test_dsmc_spatial tablengraph)

#add_executable(test_loader test_loader.cpp ${SRC_PATH}/random_init.cpp ${SRC_PATH}/boltzi_config.cpp ${SRC_PATH}/pool.cpp)
#target_include_directories(test_loader PRIVATE ${CMAKE_HOME_DIRECTORY}/../src/headers/)
#target_link_libraries(test_loader boltzi)
#target_link_libraries(test_loader tet)
#target_link_libraries(test_loader armadillo)
add_executable(test_trace test_trace.cpp)
target_include_directories(test_trace PRIVATE ${CMAKE_HOME_DIRECTORY}/../src/headers/)
