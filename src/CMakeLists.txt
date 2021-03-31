add_library(KKpipiStrongPhase ApplyCuts.cpp InitialCuts.cpp Utilities.cpp)

target_include_directories(KKpipiStrongPhase PUBLIC ../include)

target_link_libraries(KKpipiStrongPhase PUBLIC ROOT::Physics ROOT::Tree)

target_compile_definitions(KKpipiStrongPhase PUBLIC INITIAL_CUTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/InitialCuts/")
