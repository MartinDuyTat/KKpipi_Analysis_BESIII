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
#include"RooShapes/CrystalBall_Shape.h"

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
    std::string SignalMode("");
    if(TagType == "DT") {
      SignalMode = settings["MBC_Shape"].get(Name + "_SignalMode");
    }
    std::string TagMode = settings["MBC_Shape"].get(Name + "_TagMode");
    std::string RecSignalMode, RecTagMode;
    if(settings["MBC_Shape"].contains(Name + "_ReconstructedSignalMode")) {
      RecSignalMode = settings["MBC_Shape"].get(Name + "_ReconstructedSignalMode");
    } else {
      RecSignalMode = SignalMode;
    }
    if(settings["MBC_Shape"].contains(Name + "_ReconstructedTagMode")) {
      RecTagMode = settings["MBC_Shape"].get(Name + "_ReconstructedTagMode");
    } else {
      RecTagMode = TagMode;
    }
    TChain Chain(TreeName.c_str());
    std::string Filename = settings["Datasets_WithDeltaECuts"].get("SignalMC_Peaking_" + TagType);
    if(TagType == "ST") {
      Filename = Utilities::ReplaceString(Filename, "BACKGROUND", TagMode);
      Filename = Utilities::ReplaceString(Filename, "TAG", RecTagMode);
    } else if(TagType == "DT") {
      Filename = Utilities::ReplaceString(Filename, "SIGNAL1", SignalMode);
      Filename = Utilities::ReplaceString(Filename, "TAG1", TagMode);
      Filename = Utilities::ReplaceString(Filename, "SIGNAL2", RecSignalMode);
      Filename = Utilities::ReplaceString(Filename, "TAG2", RecTagMode);
      Filename = Utilities::ReplaceString(Filename, "MODE", Mode);
    }
    Chain.Add(Filename.c_str());
    RooRealVar MBC(settings.get("FitVariable").c_str(), "", settings.getD("FitRange_low"), settings.getD("FitRange_high"));
    RooDataSet Data("Data", "", &Chain, MBC);
    std::unique_ptr<FitShape> PDF;
    std::string PDFShape = settings["MBC_Shape"].get(Name + "_Shape");
    if(PDFShape == "DoubleGaussian") {
      PDF = std::unique_ptr<FitShape>{new DoubleGaussian_Shape(Name, settings["MBC_Shape"][Name + "_FitSettings"], &MBC)};
    } else if(PDFShape == "DoubleCrystalBall") {
      PDF = std::unique_ptr<FitShape>{new DoubleCrystalBall_Shape(Name, settings["MBC_Shape"][Name + "_FitSettings"], &MBC)};
    } else if(PDFShape == "CrystalBall") {
      PDF = std::unique_ptr<FitShape>{new CrystalBall_Shape(Name, settings["MBC_Shape"][Name + "_FitSettings"], &MBC)};
    } else {
      throw std::invalid_argument("Unknown peaking background shape: " + PDFShape);
    }
    RooFitResult *Result = PDF->GetPDF()->fitTo(Data, Save());
    TCanvas c("c", "", 1600, 1200);
    RooPlot *Frame = MBC.frame();
    Data.plotOn(Frame, Binning(100));
    PDF->GetPDF()->plotOn(Frame, LineColor(kBlue));
    if(TagType == "ST") {
      Frame->SetTitle((TagMode + " peaking background in " + RecTagMode + " single tag;m_{BC} (GeV);Events").c_str());
    } else if(TagType == "DT") {
      Frame->SetTitle((SignalMode  + " vs "  + TagMode + " peaking background in " + RecSignalMode + " vs " + RecTagMode + " double tag;m_{BC} (GeV);Events").c_str());
    }
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
