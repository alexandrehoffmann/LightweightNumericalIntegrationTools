find_package(PkgConfig REQUIRED)

add_executable(demo_ClensawCurtis      demo_ClensawCurtis.cpp)
add_executable(demo_GaussLegendre      demo_GaussLegendre.cpp)
add_executable(demo_HermiteAndLaguerre demo_HermiteAndLaguerre.cpp)

target_include_directories(demo_ClensawCurtis      PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(demo_GaussLegendre      PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(demo_HermiteAndLaguerre PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(demo_ClensawCurtis      PRIVATE LNIT)
target_link_libraries(demo_GaussLegendre      PRIVATE LNIT)
target_link_libraries(demo_HermiteAndLaguerre PRIVATE LNIT)

