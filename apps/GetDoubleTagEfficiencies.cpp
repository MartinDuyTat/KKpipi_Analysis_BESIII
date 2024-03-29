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
  const bool ReweightMC = settings.getB("ReweightMC");
  const bool QCMCReweighting = settings.getB("QCMCReweighting");
  std::map<int, double> CPEvenFractions;
  if(QCMCReweighting) {
    const std::string Tag = settings.get("BinningSchemeTagMode");
    for(int i = 1; i <= 8; i++) {
      const std::string Name = Tag + "_c" + std::to_string(i);
      const double ci = settings[Tag + "_BinningScheme"]["cisi"].getD(Name);
      const double FPlus = 0.5*(1 + ci*(Tag == "KLpipi" ? -1.0 : 1.0));
      CPEvenFractions.insert({i, FPlus});
    }
  }
  std::cout << "Setting up all bin combinations...\n";
  Category category(settings);
  auto BinCombinations = category.GetBinCombinations();
  auto NumberBins = BinCombinations.size();
  std::cout << "Binning ready\n";
  std::cout << "Getting the number of events generated in each bin...\n";
  TChain TruthChain("TruthTuple");
  TruthChain.Add(settings.get("TruthTupleFilename").c_str());
  std::vector<double> GeneratedEvents;
  for(const auto &Bin : BinCombinations) {
    TCut BinCut((settings.get("TagBin_variable") + "_true == " + std::to_string(Bin.second)).c_str());
    if(!(settings.contains("Inclusive_fit") && settings.getB("Inclusive_fit"))) {
      BinCut = BinCut && TCut((settings.get("SignalBin_variable") + "_true == " + std::to_string(Bin.first)).c_str());
    }
    if(ReweightMC && QCMCReweighting) {
      const double Events_CPEven = Utilities::SumWeights(&TruthChain, "ModelWeight_CPEven", std::string(BinCut.GetTitle()));
      const double Events_CPOdd = Utilities::SumWeights(&TruthChain, "ModelWeight_CPOdd", std::string(BinCut.GetTitle()));
      GeneratedEvents.push_back(Events_CPEven*(1.0 - CPEvenFractions[Bin.second]) + Events_CPOdd*CPEvenFractions[Bin.second]);
    } else if(ReweightMC) {
      const double Events = Utilities::SumWeights(&TruthChain, "ModelWeight", std::string(BinCut.GetTitle()));
      GeneratedEvents.push_back(Events);
    } else {
      const double Events = TruthChain.GetEntries(BinCut.GetTitle());
      GeneratedEvents.push_back(Events);
    }
  }
  std::cout << "True bin yields counted\n";
  std::cout << "Counting reconstructed and true bin numbers...\n";
  TFile Outfile(settings.get("EfficiencyMatrixFilename").c_str(), "RECREATE");
  TMatrixT<double> EffMatrix(NumberBins, NumberBins), EffMatrix_err(NumberBins, NumberBins);
  std::string TreeName = settings.get("TreeName");
  TChain Chain(TreeName.c_str());
  Chain.Add(settings.get("SignalMCFilename").c_str());
  double ModelWeight, ModelWeight_CPEven, ModelWeight_CPOdd, DataMCWeight;
  Chain.SetBranchAddress("ModelWeight", &ModelWeight);
  Chain.SetBranchAddress("ModelWeight_CPEven", &ModelWeight_CPEven);
  Chain.SetBranchAddress("ModelWeight_CPOdd", &ModelWeight_CPOdd);
  bool DataMCMismatchWeight = settings.getB("DataMCMismatchWeight");
  if(DataMCMismatchWeight) {
    Chain.SetBranchAddress("DataMCMismatchWeight", &DataMCWeight);
  }
  int SignalBin, SignalBin_true, TagBin, TagBin_true;
  if(settings.contains("Inclusive_fit") && settings.getB("Inclusive_fit")) {
    SignalBin = 0;
    SignalBin_true = 0;
  } else {
    Chain.SetBranchAddress(settings.get("SignalBin_variable").c_str(), &SignalBin);
    Chain.SetBranchAddress((settings.get("SignalBin_variable") + "_true").c_str(), &SignalBin_true);
  }
  Chain.SetBranchAddress(settings.get("TagBin_variable").c_str(), &TagBin);
  Chain.SetBranchAddress((settings.get("TagBin_variable") + "_true").c_str(), &TagBin_true);
  for(int i = 0; i < Chain.GetEntries(); i++) {
    Chain.GetEntry(i);
    if(DataMCMismatchWeight) {
      ModelWeight *= DataMCWeight;
      ModelWeight_CPEven *= DataMCWeight;
      ModelWeight_CPOdd *= DataMCWeight;
    }
    auto RecBin_index = std::distance(BinCombinations.begin(), std::find(BinCombinations.begin(), BinCombinations.end(), std::make_pair(SignalBin, TagBin)));
    auto TrueBin_index = std::distance(BinCombinations.begin(), std::find(BinCombinations.begin(), BinCombinations.end(), std::make_pair(SignalBin_true, TagBin_true)));
    if(ReweightMC && QCMCReweighting) {
      EffMatrix(RecBin_index, TrueBin_index) += ModelWeight_CPEven*(1.0 - CPEvenFractions[TagBin_true]) + ModelWeight_CPOdd*CPEvenFractions[TagBin_true];
    } else if(ReweightMC) {
      EffMatrix(RecBin_index, TrueBin_index) += ModelWeight;
    } else {
      EffMatrix(RecBin_index, TrueBin_index) += 1.0;
    }
  }
  std::cout << "Efficiency matrix constructed!\n";
  std::cout << "Normalizing efficiency matrix...\n";
  for(std::size_t i = 0; i < GeneratedEvents.size(); i++) {
    for(std::size_t j = 0; j < GeneratedEvents.size(); j++) {
      double p = EffMatrix(i, j)/GeneratedEvents[j];
      EffMatrix(i, j) = p;
      EffMatrix_err(i, j) = TMath::Sqrt(p*(1 - p)/GeneratedEvents[j]);
    }
  }
  EffMatrix.Write("EffMatrix");
  EffMatrix_err.Write("EffMatrix_err");
  std::cout << "Efficiency matrix ready\n";
  std::cout << "Double tag efficiency studies done!\n";
  return 0;
}
