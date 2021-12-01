cmake_minimum_required(VERSION 3.17)
project(KKpipiStrongPhase
	VERSION 1.0
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT 6.22 CONFIG REQUIRED)

set(KKPIPI_BINNED_FIT_PATH /data/lhcb/users/tat/KKpipi_Binned_Fit)

find_library(KKPIPI_BINNED_FIT_LIB KKpipiBinnedFit PATHS ${KKPIPI_BINNED_FIT_PATH}/build/src)

include_directories(${KKPIPI_BINNED_FIT_PATH}/include)
include_directories(SYSTEM /data/lhcb/users/tat/uncertainties-cpp/)
include_directories(/data/lhcb/users/tat/eigen-3.4.0/)

if(NOT KKPIPI_BINNED_FIT_LIB)
  message(FATAL_ERROR "KKpipiBinnedFit library not found")
endif()

add_compile_options(-Wall)
add_compile_options(-Wunused)
add_compile_options(-Wextra)
add_compile_options(-Werror)
add_compile_options(-O3)

add_subdirectory(${CMAKE_SOURCE_DIR}/apps)
include_directories(${CMAKE_SOURCE_DIR}/include)
add_subdirectory(${CMAKE_SOURCE_DIR}/src)

target_compile_definitions(KKpipiStrongPhase PUBLIC INITIAL_CUTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/InitialCuts/")
target_compile_definitions(KKpipiStrongPhase PUBLIC DELTAE_CUTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/DeltaECuts/")
target_compile_definitions(KKpipiStrongPhase PUBLIC TRUTH_MATCHING_CUTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/TruthMatchingCuts/")
