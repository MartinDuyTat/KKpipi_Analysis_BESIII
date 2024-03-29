set(CMAKE_INSTALL_RPATH ${KKPIPI_BINNED_FIT_PATH}/lib)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(AnalyzeTopoAna AnalyzeTopoAna.cpp)
add_executable(BinDoubleTags BinDoubleTags.cpp)
add_executable(BinMigrationStudy BinMigrationStudy.cpp)
add_executable(CorrectFlavourTagYields CorrectFlavourTagYields.cpp)
add_executable(FitDeltaE FitDeltaE.cpp)
add_executable(FitDoubleTagMBC FitDoubleTagMBC.cpp)
add_executable(FitFPlus FitFPlus.cpp)
add_executable(FitPeakingShape FitPeakingShape.cpp)
add_executable(FitSingleTagMBC FitSingleTagMBC.cpp)
add_executable(GetDCSCorrections GetDCSCorrections.cpp)
add_executable(GetDoubleTagEfficiencies GetDoubleTagEfficiencies.cpp)
add_executable(GetSingleTagEfficiencies GetSingleTagEfficiencies.cpp)
add_executable(MakeDeltaECuts MakeDeltaECuts.cpp)
add_executable(MakeResolutionHistograms MakeResolutionHistograms.cpp)
add_executable(PredictDoubleTaggedYields PredictDoubleTaggedYields.cpp)
add_executable(PrepareTagTree PrepareTagTree.cpp)

target_link_libraries(AnalyzeTopoAna PUBLIC KKpipiStrongPhase)
target_link_libraries(AnalyzeTopoAna PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(BinDoubleTags PUBLIC KKpipiStrongPhase)
target_link_libraries(BinDoubleTags PUBLIC ${KKPIPI_BINNED_FIT_LIB} -ldl)
target_link_libraries(BinDoubleTags PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(BinMigrationStudy PUBLIC KKpipiStrongPhase)
target_link_libraries(BinMigrationStudy PUBLIC ${KKPIPI_BINNED_FIT_LIB} -ldl)
target_link_libraries(BinMigrationStudy PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(CorrectFlavourTagYields PUBLIC KKpipiStrongPhase)
target_link_libraries(CorrectFlavourTagYields PUBLIC ${KKPIPI_BINNED_FIT_LIB} -ldl)
target_link_libraries(CorrectFlavourTagYields PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(FitDeltaE PUBLIC KKpipiStrongPhase)
target_link_libraries(FitDeltaE PUBLIC ${KKPIPI_BINNED_FIT_LIB} -ldl)
target_link_libraries(FitDeltaE PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(FitDoubleTagMBC PUBLIC KKpipiStrongPhase)
target_link_libraries(FitDoubleTagMBC PUBLIC ${KKPIPI_BINNED_FIT_LIB} -ldl)
target_link_libraries(FitDoubleTagMBC PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree ROOT::RooStats)

target_link_libraries(FitFPlus PUBLIC KKpipiStrongPhase)
target_link_libraries(FitFPlus PUBLIC ${KKPIPI_BINNED_FIT_LIB} -ldl)
target_link_libraries(FitFPlus PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(FitPeakingShape PUBLIC KKpipiStrongPhase)
target_link_libraries(FitPeakingShape PUBLIC ${KKPIPI_BINNED_FIT_LIB} -ldl)
target_link_libraries(FitPeakingShape PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

target_link_libraries(FitSingleTagMBC PUBLIC KKpipiStrongPhase)
target_link_libraries(FitSingleTagMBC PUBLIC ${KKPIPI_BINNED_FIT_LIB} -ldl)
target_link_libraries(FitSingleTagMBC PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree ROOT::RooStats)

target_link_libraries(GetDCSCorrections PUBLIC KKpipiStrongPhase)
target_link_libraries(GetDCSCorrections PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)
target_link_libraries(GetDCSCorrections PUBLIC ${KKPIPI_BINNED_FIT_LIB} -ldl)

target_link_libraries(GetDoubleTagEfficiencies PUBLIC KKpipiStrongPhase)
target_link_libraries(GetDoubleTagEfficiencies PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)
target_link_libraries(GetDoubleTagEfficiencies PUBLIC ${KKPIPI_BINNED_FIT_LIB} -ldl)

target_link_libraries(GetSingleTagEfficiencies PUBLIC KKpipiStrongPhase)
target_link_libraries(GetSingleTagEfficiencies PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)
target_link_libraries(GetSingleTagEfficiencies PUBLIC ${KKPIPI_BINNED_FIT_LIB} -ldl)

target_link_libraries(MakeDeltaECuts PUBLIC KKpipiStrongPhase)
target_link_libraries(MakeDeltaECuts PUBLIC ROOT::Physics)

target_link_libraries(MakeResolutionHistograms PUBLIC KKpipiStrongPhase)
target_link_libraries(MakeResolutionHistograms PUBLIC ROOT::Physics)
target_link_libraries(MakeResolutionHistograms PUBLIC ${KKPIPI_BINNED_FIT_LIB} -ldl)

target_link_libraries(PredictDoubleTaggedYields PUBLIC KKpipiStrongPhase ${KKPIPI_BINNED_FIT_LIB})

target_link_libraries(PrepareTagTree PUBLIC KKpipiStrongPhase)
target_link_libraries(PrepareTagTree PUBLIC ${KKPIPI_BINNED_FIT_LIB} -ldl)
target_link_libraries(PrepareTagTree PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

install(TARGETS AnalyzeTopoAna
		BinDoubleTags
		BinMigrationStudy
		CorrectFlavourTagYields
		FitDeltaE
		FitDoubleTagMBC
		FitFPlus
		FitPeakingShape
		FitSingleTagMBC
		GetDCSCorrections
		GetSingleTagEfficiencies
		GetDoubleTagEfficiencies
		MakeDeltaECuts
		MakeResolutionHistograms
		PredictDoubleTaggedYields
		PrepareTagTree DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../bin)
