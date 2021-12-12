// Martin Duy Tat 4th April 2021
/**
 * FitPeakingShape is an application that determines the shape and yield (relative to signal) of any peaking backgrounds for single and double tag fits
 */

#include<iostream>
#include<string>
#include<stdexcept>
#include<memory>
#include"TChain.h"
#include"TCanvas.h"
#include"RooRealVar.h"
#include"RooDataSet.h"
#include"RooFitResult.h"
#include"RooPlot.h"
#include"SingleTagYield.h"
#include"Utilities.h"
#include"Settings.h"
#include"RooShapes/FitShape.h"
#include"RooShapes/DoubleGaussian_Shape.h"
#include"RooShapes/DoubleCrystalBall_Shape.h"

int main(int argc, char *argv[]) {
  using namespace RooFit;
  Settings settings = Utilities::parse_args(argc, argv);
  std::cout << "Peaking background shape fit\n";
  std::cout << "Loading ROOT files...\n";
  std::string TreeName = settings.get("TreeName");
  std::string Mode = settings.get("Mode");
  std::string TagType = settings.get("TagType");
  int PeakingBackgrounds = settings["MBC_Shape"].getI(Mode + "_PeakingBackgrounds");
  for(int i = 0; i < PeakingBackgrounds; i++) {
    std::cout << "Fitting peaking background " << i << "\n";
    std::string Name = Mode + "_PeakingBackground" + std::to_string(i);
    TChain Chain(TreeName.c_str());
    std::string SignalMCMode = settings["MBC_Shape"].get(Name + "_Mode");
    std::string RecSignalMCMode = TagType == "ST" ? Mode : settings["MBC_Shape"].get(Name + "_RecMode");
    std::string Filename = settings["Datasets_WithDeltaECuts"].get("SignalMC_Peaking_" + TagType);
    Filename = Utilities::ReplaceString(Filename, "BACKGROUND", SignalMCMode);
    Filename = Utilities::ReplaceString(Filename, "TAG", RecSignalMCMode);
    Filename = Utilities::ReplaceString(Filename, "MODE", Mode);
    Chain.Add(Filename.c_str());
    RooRealVar MBC(settings.get("FitVariable").c_str(), "", 1.83, 1.8865);
    RooDataSet Data("Data", "", &Chain, MBC);
    std::unique_ptr<FitShape> PDF;
    if(settings["MBC_Shape"].get(Name + "_Shape") == "DoubleGaussian") {
      PDF = std::unique_ptr<FitShape>{new DoubleGaussian_Shape(Name, settings["MBC_Shape"][Name + "_FitSettings"], &MBC)};
    } else if(settings["MBC_Shape"].get(Name + "_Shape") == "DoubleCrystalBall") {
      PDF = std::unique_ptr<FitShape>{new DoubleCrystalBall_Shape(Name, settings["MBC_Shape"][Name + "_FitSettings"], &MBC)};
    } else {
      throw std::invalid_argument("Unknown peaking background shape");
    }
    RooFitResult *Result = PDF->GetPDF()->fitTo(Data, Save());
    TCanvas c("c", "", 1600, 1200);
    RooPlot *Frame = MBC.frame();
    Data.plotOn(Frame, Binning(100));
    PDF->GetPDF()->plotOn(Frame, LineColor(kBlue));
    Frame->SetTitle((SignalMCMode + " peaking background in " + Mode + " single tag;m_{BC} (GeV);Events").c_str());
    Frame->Draw();
    c.SaveAs(settings["MBC_Shape"].get(Name + "_PlotFilename").c_str());
    Result->Print();
    std::cout << "Fit results for peaking background " << i << ":\n";
    std::cout << "Status " << Result->status() << "\n";
    std::cout << "covQual " << Result->covQual() << "\n";
    RooArgList floating_param = Result->floatParsFinal();
    for(int i = 0; i < floating_param.getSize(); i++) {
      RooRealVar *param = static_cast<RooRealVar*>(floating_param.at(i));
      std::cout << param->GetName() << " " << param->getVal() << "\n";
    }
    std::cout << "\n";
  }
  std::cout << "Peaking backgrounds are now accounted for" << "\n";
  return 0;
}
