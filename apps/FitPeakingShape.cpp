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
#include"RooArgList.h"
#include"RooDataSet.h"
#include"RooFitResult.h"
#include"RooPlot.h"
#include"RooArgusBG.h"
#include"SingleTagYield.h"
#include"Utilities.h"
#include"Unique.h"
#include"Settings.h"
#include"RooShapes/FitShape.h"
#include"RooShapes/DoubleGaussian_Shape.h"
#include"RooShapes/DoubleCrystalBall_Shape.h"
#include"RooShapes/CrystalBall_Shape.h"
#include"RooShapes/Chebychev_Shape.h"

int main(int argc, char *argv[]) {
  using namespace RooFit;
  Settings settings = Utilities::parse_args(argc, argv);
  std::cout << "Peaking background shape fit\n";
  std::cout << "Loading ROOT files...\n";
  std::string TreeName = settings.get("TreeName");
  std::string Mode = settings.get("Mode");
  std::string TagType = settings.get("TagType");
  int PeakingBackgrounds = settings["MBC_Shape"].getI(Mode + "_PeakingBackgrounds");
  std::vector<RooFitResult*> Results;
  std::vector<double> Yields;
  for(int i = 0; i < PeakingBackgrounds; i++) {
    std::string Name = Mode + "_PeakingBackground" + std::to_string(i);
    if(settings["MBC_Shape"].contains(Name + "_FitShape") && !settings["MBC_Shape"].getB(Name + "_FitShape")) {
      Results.push_back(nullptr);
      continue;
    }
    std::cout << "Fitting peaking background " << i << "\n";
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
    std::string Filename;
    if(settings["MBC_Shape"].contains(Name + "_Filename")) {
      Filename = settings["MBC_Shape"].get(Name + "_Filename");
    } else {
      Filename = settings["Datasets_WithDeltaECuts"].get("SignalMC_Peaking_" + TagType);
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
    }
    TChain Chain(TreeName.c_str());
    Chain.Add(Filename.c_str());
    RooRealVar MBC(settings.get("FitVariable").c_str(), "", settings.getD("FitRange_low"), settings.getD("FitRange_high"));
    RooArgSet Variables;
    Variables.add(MBC);
    std::string WeightName("");
    RooRealVar WeightVar;
    if(settings["MBC_Shape"].contains(Name + "_Weight")) {
      WeightName = settings["MBC_Shape"].get(Name + "_Weight");
      WeightVar = RooRealVar(WeightName.c_str(), "", 0.0, 1.0);
      Variables.add(WeightVar);
    }
    std::string Cut("");
    RooRealVar iDcyTr("iDcyTr", "", 0, 100000);
    if(settings["MBC_Shape"].contains(Name + "_ComponentsToIgnore")) {
      Variables.add(iDcyTr);
      auto ComponentsToIgnore = Utilities::ConvertStringToVector(settings["MBC_Shape"].get(Name + "_ComponentsToIgnore"));
      for(const auto ComponentToIgnore : ComponentsToIgnore) {
	Cut += "iDcyTr != " + ComponentToIgnore;
	if(!(ComponentToIgnore == ComponentsToIgnore.back())) {
	  Cut += " && ";
	}
      }
    }
    RooDataSet Data("Data", "", &Chain, Variables, Cut.c_str(), WeightName.c_str());
    std::unique_ptr<FitShape> PDF;
    std::string PDFShape = settings["MBC_Shape"].get(Name + "_Shape");
    if(PDFShape == "DoubleGaussian") {
      PDF = std::unique_ptr<FitShape>{new DoubleGaussian_Shape(Name, settings["MBC_Shape"][Name + "_FitSettings"], &MBC)};
    } else if(PDFShape == "DoubleCrystalBall") {
      PDF = std::unique_ptr<FitShape>{new DoubleCrystalBall_Shape(Name, settings["MBC_Shape"][Name + "_FitSettings"], &MBC)};
    } else if(PDFShape == "CrystalBall") {
      PDF = std::unique_ptr<FitShape>{new CrystalBall_Shape(Name, settings["MBC_Shape"][Name + "_FitSettings"], &MBC)};
    } else if(PDFShape == "Chebychev") {
      PDF = std::unique_ptr<FitShape>{new Chebychev_Shape(Name, settings["MBC_Shape"][Name + "_FitSettings"], &MBC)};
    } else {
      throw std::invalid_argument("Unknown peaking background shape: " + PDFShape);
    }
    RooAbsPdf *Model = nullptr;
    if(settings["MBC_Shape"][Name + "_FitSettings"].contains(Name + "_ArgusBackground") &&
       settings["MBC_Shape"][Name + "_FitSettings"].getB(Name + "_ArgusBackground")) {
      auto Nsig = Utilities::load_param(settings["MBC_Shape"][Name + "_FitSettings"], Name + "_Nsig");
      auto Nbkg = Utilities::load_param(settings["MBC_Shape"][Name + "_FitSettings"], Name + "_Nbkg");
      auto c = Utilities::load_param(settings["MBC_Shape"][Name + "_FitSettings"], Name + "_c");
      auto End = Unique::create<RooRealVar*>("End", "", 1.8865);
      auto BackgroundModel = Unique::create<RooArgusBG*>("ArgusBackground", "", MBC, *End, *c);
      auto SignalModel = PDF->GetPDF();
      Model = Unique::create<RooAddPdf*>("Model", "", RooArgList(*SignalModel, *BackgroundModel), RooArgList(*Nsig, *Nbkg));
    } else {
      Model = PDF->GetPDF();
    }
    Results.push_back(Model->fitTo(Data, Save(), SumW2Error(false)));
    if(WeightName == "") {
      Yields.push_back(Chain.GetEntries());
    } else {
      double Total = 0.0;
      double ChainWeight;
      Chain.SetBranchAddress(WeightName.c_str(), &ChainWeight);
      for(int j = 0; j < Chain.GetEntries(); j++) {
	Chain.GetEntry(j);
	Total += ChainWeight;
      }
      Yields.push_back(Total);
    }
    TCanvas c("c", "", 1600, 1200);
    RooPlot *Frame = MBC.frame();
    Data.plotOn(Frame, Binning(100));
    Model->plotOn(Frame, LineColor(kBlue));
    if(settings["MBC_Shape"][Name + "_FitSettings"].contains(Name + "_ArgusBackground") &&
       settings["MBC_Shape"][Name + "_FitSettings"].getB(Name + "_ArgusBackground")) {
      Model->plotOn(Frame, LineColor(kBlue), LineStyle(kDashed), Components((Name + "_" + PDFShape).c_str()));
    }
    if(TagType == "ST") {
      Frame->SetTitle((TagMode + " peaking background in " + RecTagMode + " single tag;m_{BC} (GeV);Events").c_str());
    } else if(TagType == "DT") {
      Frame->SetTitle((SignalMode  + " vs "  + TagMode + " peaking background in " + RecSignalMode + " vs " + RecTagMode + " double tag;m_{BC} (GeV);Events").c_str());
    }
    Frame->Draw();
    c.SaveAs(settings["MBC_Shape"].get(Name + "_PlotFilename").c_str());
  }
  int n = 0;
  for(auto Result : Results) {
    if(!Result) {
      n++;
      continue;
    }
    Result->Print();
    std::cout << "Fit results for peaking background " << n << ":\n";
    std::cout << "Status " << Result->status() << "\n";
    std::cout << "covQual " << Result->covQual() << "\n";
    RooArgList floating_param = Result->floatParsFinal();
    for(int i = 0; i < floating_param.getSize(); i++) {
      RooRealVar *param = static_cast<RooRealVar*>(floating_param.at(i));
      std::string ParameterName(param->GetName());
      std::cout << ParameterName << " " << param->getVal() << "\n";
      if(ParameterName.substr(ParameterName.length() - 4, ParameterName.length()) == "Nsig") {
	std::cout << ParameterName + "_err " << param->getError() << "\n";
      }
    }
    std::cout << Mode << "_PeakingBackground" << n << "_Yield " << Yields[n] << "\n";
    std::cout << Mode << "_PeakingBackground" << n << "_Yield_err " << TMath::Sqrt(Yields[n]) << "\n";
    std::cout << "\n";
    n++;
  }
  std::cout << "Peaking backgrounds are now accounted for" << "\n";
  return 0;
}
