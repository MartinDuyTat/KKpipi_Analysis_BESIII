add_library(KKpipiStrongPhase ApplyCuts.cpp DeltaECut.cpp DeltaEFit.cpp InitialCuts.cpp Utilities.cpp)

target_include_directories(KKpipiStrongPhase PUBLIC ../include)

target_link_libraries(KKpipiStrongPhase PUBLIC ROOT::Physics ROOT::Tree ROOT::RooFit ROOT::Gpad)
