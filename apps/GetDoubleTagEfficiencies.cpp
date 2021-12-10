// Martin Duy Tat 26th November 2021
/**
 * GetDoubleTagEfficiencies is an application that calculates the double tag efficiency matrices of double tags
 * First it counts the number of generated events in each bin, then it finds the number of reconstructed events in each bin
 */

#include<iostream>
#include<map>
#include<string>
#include"TChain.h"
#include"TFile.h"
#include"TMatrixT.h"
#include"Utilities.h"
#include"Settings.h"

/**
 * This functor class converts bin -8, -7, ..., -1, +1, ..., +8 to indices 0, 1, 2, ..., 15
 * @param Bin Bin number
 */
class BinIndex {
  public:
    BinIndex(int NumberBins): m_NumberBins(NumberBins) {}
    int operator ()(int Bin) { return Bin < 0 ? Bin + m_NumberBins : Bin + m_NumberBins - 1; }
  private:
    const int m_NumberBins;
};

int main(int argc, char *argv[]) {
  std::cout << "Calculating double tag efficiency matrix from signal MC\n";
  Settings settings = Utilities::parse_args(argc, argv);
  int NumberBins = settings["BinningScheme"].getI("NumberBins");
  std::cout << "First get the number of events generated in each bin...\n";
  TChain TruthChain("TruthTuple");
  TruthChain.Add(settings.get("TruthTupleFilename").c_str());
  std::map<int, int> GeneratedEvents;
  for(int i = -NumberBins; i <= NumberBins; i++) {
    if(i == 0) {
      continue;
    }
    int Events = TruthChain.GetEntries((settings.get("SignalBin_variable") + "_true == " + std::to_string(i)).c_str());
    GeneratedEvents.insert({i, Events});
  }
  std::cout << "True bin yields counted\n";
  std::cout << "Second get the number of reconstructed events in each bin and construct an efficiency matrix\n";
  TFile Outfile(settings.get("EfficiencyMatrixFilename").c_str(), "RECREATE");
  // For now I assume flavour tags, so we need both positive and negative bins, but this will be generalized in the future
  TMatrixT<double> EffMatrix(2*NumberBins, 2*NumberBins);
  std::string TreeName = settings.get("TreeName");
  TChain Chain(TreeName.c_str());
  Chain.Add(settings.get("SignalMCFilename").c_str());
  int SignalBin, SignalBin_true;
  Chain.SetBranchAddress(settings.get("SignalBin_variable").c_str(), &SignalBin);
  Chain.SetBranchAddress((settings.get("SignalBin_variable") + "_true").c_str(), &SignalBin_true);
  BinIndex binIndex(NumberBins);
  std::cout << "Counting reconstructed and true bin numbers...\n";
  for(int i = 0; i < Chain.GetEntries(); i++) {
    Chain.GetEntry(i);
    EffMatrix(binIndex(SignalBin), binIndex(SignalBin_true))++;
  }
  // Finally normalize all elements by the number of generated events
  for(int i = -NumberBins; i <= NumberBins; i++) {
    if(i == 0) {
      continue;
    }
    for(int j = -NumberBins; j <= NumberBins; j++) {
      if(j == 0) {
	continue;
      }
      EffMatrix(binIndex(i), binIndex(j)) /= GeneratedEvents[j];
    }
  }
  std::cout << "Efficiency matrix constructed!\n";
  EffMatrix.Write("EffMatrix");
  std::cout << "Double tag yields (copy to some file for now):\n";
  std::string Mode = settings.get("Mode");
  double SignalMC_SampleSize = static_cast<double>(settings["SignalMCSampleSize"].getI("DoubleTag"));
  double Entries = static_cast<double>(Chain.GetEntries());
  double Efficiency = Entries/SignalMC_SampleSize;
  std::cout << Mode << "_DoubleTagEfficiency     " << Efficiency << "\n";
  std::cout << Mode << "_DoubleTagEfficiency_err " << TMath::Sqrt(Efficiency*(1 - Efficiency)/SignalMC_SampleSize) << "\n";
  std::cout << "Double tag efficiency studies done!\n";
  return 0;
}
