add_library(KKpipiStrongPhase
	    ApplyCuts.cpp
	    Bes3plotstyle.cpp
	    BinnedDataLoader.cpp
	    BinnedCPTagYieldPrediction.cpp
	    BinnedDTData.cpp
	    BinnedDTYieldPrediction.cpp
	    BinnedFitModel.cpp
	    Category.cpp
	    CholeskySmearing.cpp
	    cisiFitter.cpp
	    cisiK0pipi.cpp
	    cisiLikelihood.cpp
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
	    RawBinnedCPTagYields.cpp
	    RawBinnedDTYields.cpp
	    RooShapes/Chebychev_Shape.cpp
	    RooShapes/CrystalBall_Shape.cpp
	    RooShapes/DoubleCrystalBall_Shape.cpp 
	    RooShapes/DoubleGaussian_Shape.cpp 
	    RooShapes/DoubleGaussianRatio_Shape.cpp 
	    RooShapes/DoublePolynomial_Shape.cpp
	    RooShapes/FitShape.cpp)

target_include_directories(KKpipiStrongPhase PUBLIC ../include)

target_link_libraries(KKpipiStrongPhase PUBLIC -ldl)

target_link_libraries(KKpipiStrongPhase PUBLIC OpenMP::OpenMP_CXX)

target_link_libraries(KKpipiStrongPhase PUBLIC ROOT::Physics ROOT::Tree ROOT::RooFit ROOT::Gpad)
