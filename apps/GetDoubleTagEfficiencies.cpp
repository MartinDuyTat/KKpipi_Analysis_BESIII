// Martin Duy Tat 26th November 2021
/**
 * GetDoubleTagEfficiencies is an application that calculates the double tag efficiency matrices of double tags
 * First it counts the number of generated events in each bin, then it finds the number of reconstructed events in each bin
 */

#include<iostream>
#include<map>
#include<string>
#include<utility>
#include<iterator>
#include"TChain.h"
#include"TFile.h"
#include"TMatrixT.h"
#include"Utilities.h"
#include"Settings.h"
#include"Category.h"

int main(int argc, char *argv[]) {
  std::cout << "Calculating double tag efficiency matrix from signal MC\n";
  Settings settings = Utilities::parse_args(argc, argv);
  const bool K0pipiQCMC = settings.contains("K0pipi_QCMC_TagMode");
  std::string K0pipiMode;
  if(K0pipiQCMC) {
    K0pipiMode = settings.get("K0pipi_QCMC_TagMode");
    std::cout << "Using QCMC reweighting for " << K0pipiMode << "\n";
  }
  std::cout << "Setting up all bin combinations...\n";
  Category category(settings);
  auto BinCombinations = category.GetBinCombinations();
  auto NumberBins = BinCombinations.size();
  std::cout << "Binning ready\n";
  std::cout << "Getting the number of events generated in each bin...\n";
  TChain TruthChain("TruthTuple");
  TruthChain.Add(settings.get("TruthTupleFilename").c_str());
  std::vector<double> GeneratedEvents_CPEven,
                      GeneratedEvents_CPOdd,
                      GeneratedEvents_K0pipi,
                      GeneratedEvents;
  for(const auto &Bin : BinCombinations) {
    const std::string BinCutString = settings.get("TagBin_variable")
                                   + "_true == "
                                   + std::to_string(Bin.second);
    TCut BinCut(BinCutString.c_str());
    if(!(settings.contains("Inclusive_fit") && settings.getB("Inclusive_fit"))) {
      const std::string SignalBinCut = settings.get("SignalBin_variable")
                                     + "_true == " + std::to_string(Bin.first);
      BinCut = BinCut && TCut(SignalBinCut.c_str());
    }
    const double Events_CPEven =
      Utilities::SumWeights(&TruthChain,
			    "ModelWeight_CPEven",
			    std::string(BinCut.GetTitle()));
    const double Events_CPOdd =
      Utilities::SumWeights(&TruthChain,
			    "ModelWeight_CPOdd",
			    std::string(BinCut.GetTitle()));
    const std::string WeightNameK0pipi = K0pipiQCMC ? 
      "ModelWeight_" + K0pipiMode + "_TagBin" + std::to_string(Bin.second) : "";
    const double Events_K0pipi =
      Utilities::SumWeights(&TruthChain,
			    WeightNameK0pipi,
			    std::string(BinCut.GetTitle()));
    const double Events =
      Utilities::SumWeights(&TruthChain, "", std::string(BinCut.GetTitle()));
    GeneratedEvents_CPEven.push_back(Events_CPEven);
    GeneratedEvents_CPOdd.push_back(Events_CPOdd);
    GeneratedEvents_K0pipi.push_back(Events_K0pipi);
    GeneratedEvents.push_back(Events);
  }
  std::cout << "True bin yields counted\n";
  std::cout << "Counting reconstructed and true bin numbers...\n";
  TFile Outfile(settings.get("EfficiencyMatrixFilename").c_str(), "RECREATE");
  TMatrixT<double> EffMatrix_CPEven(NumberBins, NumberBins),
                   EffMatrix_CPEven_err(NumberBins, NumberBins),
                   EffMatrix_CPOdd(NumberBins, NumberBins),
                   EffMatrix_CPOdd_err(NumberBins, NumberBins),
                   EffMatrix_K0pipi(NumberBins, NumberBins),
                   EffMatrix_K0pipi_err(NumberBins, NumberBins),
                   EffMatrix(NumberBins, NumberBins),
		   EffMatrix_err(NumberBins, NumberBins);
  std::string TreeName = settings.get("TreeName");
  TChain Chain(TreeName.c_str());
  Chain.Add(settings.get("SignalMCFilename").c_str());
  double ModelWeight_CPEven = 1.0,
         ModelWeight_CPOdd = 1.0;
  Chain.SetBranchAddress("ModelWeight_CPEven", &ModelWeight_CPEven);
  Chain.SetBranchAddress("ModelWeight_CPOdd", &ModelWeight_CPOdd);
  std::map<int, double> K0pipiWeights;
  if(K0pipiQCMC) {
    for(int TagBin = 1; TagBin <= 8; TagBin++) {
      K0pipiWeights.insert({TagBin, 0.0});
      const std::string WeightName = "ModelWeight_" + K0pipiMode
	                           + "_TagBin" + std::to_string(TagBin);
      Chain.SetBranchAddress(WeightName.c_str(), &K0pipiWeights[TagBin]);
    }
  }
  int SignalBin = 0,
      SignalBin_true = 0,
      TagBin = 0,
      TagBin_true = 0;
  if(!(settings.contains("Inclusive_fit") && settings.getB("Inclusive_fit"))) {
    Chain.SetBranchAddress(settings.get("SignalBin_variable").c_str(),
			   &SignalBin);
    Chain.SetBranchAddress((settings.get("SignalBin_variable") + "_true").c_str(),
			   &SignalBin_true);
  }
  Chain.SetBranchAddress(settings.get("TagBin_variable").c_str(), &TagBin);
  Chain.SetBranchAddress((settings.get("TagBin_variable") + "_true").c_str(),
			 &TagBin_true);
  for(int i = 0; i < Chain.GetEntries(); i++) {
    Chain.GetEntry(i);
    auto RecBin_index =
      std::distance(BinCombinations.begin(),
		    std::find(BinCombinations.begin(),
			      BinCombinations.end(),
			      std::make_pair(SignalBin, TagBin)));
    auto TrueBin_index =
      std::distance(BinCombinations.begin(),
		    std::find(BinCombinations.begin(),
			      BinCombinations.end(),
			      std::make_pair(SignalBin_true, TagBin_true)));
    EffMatrix_CPEven(RecBin_index, TrueBin_index) += ModelWeight_CPEven;
    EffMatrix_CPOdd(RecBin_index, TrueBin_index) += ModelWeight_CPOdd;
    if(K0pipiQCMC) {
      EffMatrix_K0pipi(RecBin_index, TrueBin_index) += K0pipiWeights[TagBin_true];
    } else {
      EffMatrix_K0pipi(RecBin_index, TrueBin_index) += 1.0;
    }
    EffMatrix(RecBin_index, TrueBin_index) += 1.0;
  }
  std::cout << "Efficiency matrix constructed!\n";
  std::cout << "Normalizing efficiency matrix...\n";
  for(std::size_t i = 0; i < GeneratedEvents_CPEven.size(); i++) {
    for(std::size_t j = 0; j < GeneratedEvents_CPEven.size(); j++) {
      double p_CPEven = EffMatrix_CPEven(i, j)/GeneratedEvents_CPEven[j];
      double p_CPOdd = EffMatrix_CPOdd(i, j)/GeneratedEvents_CPOdd[j];
      double p_K0pipi = EffMatrix_K0pipi(i, j)/GeneratedEvents_K0pipi[j];
      double p = EffMatrix(i, j)/GeneratedEvents[j];
      EffMatrix_CPEven(i, j) = p_CPEven;
      EffMatrix_CPOdd(i, j) = p_CPOdd;
      EffMatrix_K0pipi(i, j) = p_K0pipi;
      EffMatrix(i, j) = p;
      EffMatrix_CPEven_err(i, j) = TMath::Sqrt(p_CPEven*(1 - p_CPEven)/
					       GeneratedEvents_CPEven[j]);
      EffMatrix_CPOdd_err(i, j) = TMath::Sqrt(p_CPOdd*(1 - p_CPOdd)/
					      GeneratedEvents_CPOdd[j]);
      EffMatrix_K0pipi_err(i, j) = TMath::Sqrt(p_K0pipi*(1 - p_K0pipi)/
                                               GeneratedEvents_K0pipi[j]);
      EffMatrix_err(i, j) = TMath::Sqrt(p*(1.0 - p)/GeneratedEvents[j]);
    }
  }
  EffMatrix_CPEven.Write("EffMatrix_CPEven");
  EffMatrix_CPEven_err.Write("EffMatrix_CPEven_err");
  EffMatrix_CPOdd.Write("EffMatrix_CPOdd");
  EffMatrix_CPOdd_err.Write("EffMatrix_CPOdd_err");
  EffMatrix_K0pipi.Write("EffMatrix_K0pipi");
  EffMatrix_K0pipi_err.Write("EffMatrix_K0pipi_err");
  EffMatrix.Write("EffMatrix");
  EffMatrix_err.Write("EffMatrix_err");
  std::cout << "Efficiency matrix ready\n";
  std::cout << "Double tag efficiency studies done!\n";
  return 0;
}
