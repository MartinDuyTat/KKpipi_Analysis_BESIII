// Martin Duy Tat 25th November 2021

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
  m_DataSet = new RooDataSet("InputData", "", m_Tree, RooArgSet(*m_SignalMBC, *m_TagMBC, SignalBin, TagBin), "TagMBC > 1.86 && TagMBC < 1.87");
  RooCategory *CategoryVariable = m_Category.GetCategoryVariable();
  RooDataSet CategorySet("InputData_Category", "", RooArgSet(*CategoryVariable));
  for(int i = 0; i < m_DataSet->numEntries(); i++) {
    const RooArgSet *Row = m_DataSet->get(i);
    int SignalBinNumber = static_cast<RooRealVar*>(Row->find(SignalBin_Name.c_str()))->getVal();
    int TagBinNumber = static_cast<RooRealVar*>(Row->find(TagBin_Name.c_str()))->getVal();
    std::string CategoryString = m_Category(SignalBinNumber, TagBinNumber);
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
