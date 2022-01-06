// Martin Duy Tat 25th November 2021

#include<string>
#include<memory>
#include"TTree.h"
#include"RooDataSet.h"
#include"RooArgSet.h"
#include"RooCategory.h"
#include"BinnedDataLoader.h"
#include"Settings.h"
#include"Category.h"

BinnedDataLoader::BinnedDataLoader(const Settings &settings,
				   TTree *Tree,
				   RooRealVar *SignalMBC,
				   RooRealVar *TagMBC): m_Settings(settings),
							m_Tree(Tree),
							m_SignalMBC(SignalMBC),
							m_TagMBC(TagMBC),
							m_Category(m_Settings) {
  MakeDataSet();
}

BinnedDataLoader::~BinnedDataLoader() {
  if(m_DataSet) {
    delete m_DataSet;
  }
}

void BinnedDataLoader::MakeDataSet() {
  std::string SignalBin_Name = m_Settings.get("SignalBin_variable");
  std::string TagBin_Name = m_Settings.get("TagBin_variable");
  RooRealVar SignalBin(SignalBin_Name.c_str(), "", -8, 8);
  RooRealVar TagBin(TagBin_Name.c_str(), "", -8, 8);
  RooArgSet Variables(*m_SignalMBC, *m_TagMBC, SignalBin, TagBin);
  std::unique_ptr<RooRealVar> InvMassVar;
  std::string MassCut("");
  if(m_Settings.contains("InvariantMassVariable")) {
    std::string MassVarName = m_Settings.get("InvariantMassVariable");
    m_Tree->SetBranchStatus(MassVarName.c_str(), 1);
    double LowMassCut = m_Settings.getD("InvariantMassVariable_low");
    double HighMassCut = m_Settings.getD("InvariantMassVariable_high");
    InvMassVar = std::unique_ptr<RooRealVar>(new RooRealVar(MassVarName.c_str(), "", LowMassCut, HighMassCut));
    MassCut = "(" + MassVarName + " > " + std::to_string(LowMassCut) + " && " + MassVarName + " < " + std::to_string(HighMassCut) + ")";
    Variables.add(*InvMassVar);
  }
  m_DataSet = new RooDataSet("InputData", "", m_Tree, Variables, MassCut.c_str());
  RooCategory *CategoryVariable = m_Category.GetCategoryVariable();
  RooDataSet CategorySet("InputData_Category", "", RooArgSet(*CategoryVariable));
  for(int i = 0; i < m_DataSet->numEntries(); i++) {
    const RooArgSet *Row = m_DataSet->get(i);
    std::string CategoryString;
    int SignalBinNumber;
    if(m_Settings.contains("Inclusive_fit") && m_Settings.getB("Inclusive_fit")) {
      SignalBinNumber = 0;
    } else {
      SignalBinNumber = static_cast<RooRealVar*>(Row->find(SignalBin_Name.c_str()))->getVal();
    }
    int TagBinNumber = static_cast<RooRealVar*>(Row->find(TagBin_Name.c_str()))->getVal();
    CategoryString = m_Category(SignalBinNumber, TagBinNumber);
    CategoryVariable->setLabel(CategoryString.c_str());
    CategorySet.add(RooArgSet(*CategoryVariable));
  }
  m_DataSet->merge(&CategorySet);
}

RooDataSet* BinnedDataLoader::GetDataSet() {
  return m_DataSet;
}

Category* BinnedDataLoader::GetCategoryObject() {
  return &m_Category;
}
