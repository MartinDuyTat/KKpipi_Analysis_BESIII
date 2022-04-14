add_library(KKpipiStrongPhase
	    ApplyCuts.cpp
	    BinnedDataLoader.cpp
	    BinnedFitModel.cpp
	    Category.cpp
	    cisiK0pipi.cpp
	    CutsFromFile.cpp
	    DeltaECut.cpp
	    DeltaEFit.cpp
	    DeltaEFitModel.cpp
	    DoubleTagYield.cpp
	    FPlusFitter.cpp
	    InitialCuts.cpp
	    PredictNumberEvents.cpp
	    Settings.cpp
	    SingleTagYield.cpp
	    TopoAnaReader.cpp
	    TruthMatchingCuts.cpp
	    Utilities.cpp
	    HadronicParameters/DCS_Parameters.cpp
	    HadronicParameters/Ki.cpp
	    HadronicParameters/cisi.cpp
	    PhaseSpace/DalitzUtilities.cpp
	    PhaseSpace/KKpipi_PhaseSpace.cpp
	    PhaseSpace/KKpipi_vs_CP_PhaseSpace.cpp
	    PhaseSpace/KKpipi_vs_Flavour_PhaseSpace.cpp
	    PhaseSpace/KKpipi_vs_K0hh_PhaseSpace.cpp
	    RooShapes/Chebychev_Shape.cpp
	    RooShapes/CrystalBall_Shape.cpp
	    RooShapes/DoubleCrystalBall_Shape.cpp 
	    RooShapes/DoubleGaussian_Shape.cpp 
	    RooShapes/DoublePolynomial_Shape.cpp
	    RooShapes/FitShape.cpp)

target_include_directories(KKpipiStrongPhase PUBLIC ../include)

target_link_libraries(KKpipiStrongPhase PUBLIC ROOT::Physics ROOT::Tree ROOT::RooFit ROOT::Gpad)
