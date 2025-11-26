find_package(PkgConfig REQUIRED)

add_executable(demo_ClensawCurtis       demo_ClensawCurtis.cpp)
add_executable(demo_ClensawCurtisHybrid demo_ClensawCurtisHybrid.cpp)
add_executable(demo_GaussLegendre       demo_GaussLegendre.cpp)
add_executable(demo_GLCC                demo_GLCC.cpp)
add_executable(demo_HermiteAndLaguerre  demo_HermiteAndLaguerre.cpp)
add_executable(demo_LevermooreLikePDF   demo_LevermooreLikePDF.cpp)

target_include_directories(demo_ClensawCurtis       PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(demo_ClensawCurtisHybrid PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(demo_GaussLegendre       PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(demo_GLCC                PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(demo_HermiteAndLaguerre  PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_include_directories(demo_LevermooreLikePDF   PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(demo_ClensawCurtis       PRIVATE LNIT)
target_link_libraries(demo_ClensawCurtisHybrid PRIVATE LNIT)
target_link_libraries(demo_GaussLegendre       PRIVATE LNIT)
target_link_libraries(demo_GLCC                PRIVATE LNIT)
target_link_libraries(demo_HermiteAndLaguerre  PRIVATE LNIT)
target_link_libraries(demo_LevermooreLikePDF   PRIVATE LNIT)

