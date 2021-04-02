add_executable(FitDeltaE FitDeltaE.cpp)
add_executable(MakeDeltaECuts MakeDeltaECuts.cpp)
add_executable(PrepareTagTree PrepareTagTree.cpp)

target_link_libraries(FitDeltaE PUBLIC KKpipiStrongPhase)
target_link_libraries(FitDeltaE PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(MakeDeltaECuts PUBLIC KKpipiStrongPhase)
target_link_libraries(MakeDeltaECuts PUBLIC ROOT::Physics)

target_link_libraries(PrepareTagTree PUBLIC KKpipiStrongPhase)
target_link_libraries(PrepareTagTree PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

install(TARGETS FitDeltaE MakeDeltaECuts PrepareTagTree DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../bin)
