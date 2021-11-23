set(CMAKE_INSTALL_RPATH ${KKPIPI_BINNED_FIT_PATH}/lib)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(AnalyzeTopoAna AnalyzeTopoAna.cpp)
add_executable(DetermineDoubleTagYields DetermineDoubleTagYields.cpp)
add_executable(FitDeltaE FitDeltaE.cpp)
add_executable(FitSingleTagPeaking FitSingleTagPeaking.cpp)
add_executable(FitSingleTagMBC FitSingleTagMBC.cpp)
add_executable(MakeDeltaECuts MakeDeltaECuts.cpp)
add_executable(PredictDoubleTaggedYields PredictDoubleTaggedYields.cpp)
add_executable(PrepareTagTree PrepareTagTree.cpp)

target_link_libraries(AnalyzeTopoAna PUBLIC KKpipiStrongPhase)
target_link_libraries(AnalyzeTopoAna PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(DetermineDoubleTagYields PUBLIC KKpipiStrongPhase)
target_link_libraries(DetermineDoubleTagYields PUBLIC ${KKPIPI_BINNED_FIT_LIB} -ldl)
target_link_libraries(DetermineDoubleTagYields PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(FitDeltaE PUBLIC KKpipiStrongPhase)
target_link_libraries(FitDeltaE PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(FitSingleTagPeaking PUBLIC KKpipiStrongPhase)
target_link_libraries(FitSingleTagPeaking PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(FitSingleTagMBC PUBLIC KKpipiStrongPhase)
target_link_libraries(FitSingleTagMBC PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(MakeDeltaECuts PUBLIC KKpipiStrongPhase)
target_link_libraries(MakeDeltaECuts PUBLIC ROOT::Physics)

target_link_libraries(PredictDoubleTaggedYields PUBLIC KKpipiStrongPhase ${KKPIPI_BINNED_FIT_LIB})

target_link_libraries(PrepareTagTree PUBLIC KKpipiStrongPhase)
target_link_libraries(PrepareTagTree PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

install(TARGETS AnalyzeTopoAna DetermineDoubleTagYields FitDeltaE FitSingleTagPeaking FitSingleTagMBC MakeDeltaECuts PredictDoubleTaggedYields PrepareTagTree DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../bin)
