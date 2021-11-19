add_library(KKpipiStrongPhase
	    ApplyCuts.cpp
	    CutsFromFile.cpp
	    DeltaECut.cpp
	    DeltaEFit.cpp
	    DeltaEFitModel.cpp
	    DoubleTagYield.cpp
	    InitialCuts.cpp
	    PredictNumberEvents.cpp
	    Settings.cpp
	    SingleTagYield.cpp
	    TopoAnaReader.cpp
	    TruthMatchingCuts.cpp
	    Utilities.cpp
	    RooShapes/Chebychev_Shape.cpp
	    RooShapes/DoubleGaussian_Shape.cpp 
	    RooShapes/DoublePolynomial_Shape.cpp
	    RooShapes/FitShape.cpp)

target_include_directories(KKpipiStrongPhase PUBLIC ../include)

target_link_libraries(KKpipiStrongPhase PUBLIC ROOT::Physics ROOT::Tree ROOT::RooFit ROOT::Gpad)
