include(FetchContent)

FetchContent_Declare(
  benchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG v1.8.3 # or latest stable release
)

# Disable benchmark's own tests (requires gtest, which we don't need)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(benchmark)

find_package(PkgConfig REQUIRED)

add_executable(benchmark_addaptiveQuadratures  benchmark_addaptiveQuadratures.cpp)

target_include_directories(benchmark_addaptiveQuadratures PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(benchmark_addaptiveQuadratures PRIVATE LNIT benchmark::benchmark)
