// Martin Duy Tat 25th November 2021
/**
 * BinnedDataLoader takes in a TTree and separates the events into bins, before creating a RooDataSet with the correct category for each event
 */

#ifndef BINNEDDATALOADER
#define BINNEDDATALOADER

#include"TTree.h"
#include"RooRealVar.h"
#include"Settings.h"
#include"Category.h"

class BinnedDataLoader {
  public:
    /**
     * Constructor that saves a pointer to the TTree and the settings file
     * The TTree is then binned and saved to a RooDataSet
     * @param settings Contains all the fit settings
     * @param Tree The TTree with all the double tag events
     * @param SignalMBC The fit variable
     */
    BinnedDataLoader(const Settings &settings, TTree *Tree, RooRealVar *SignalMBC, RooRealVar *TagMBC);
    /**
     * Destructor that deletes dataset
     */
    ~BinnedDataLoader();
    /**
     * Get the binned dataset
     */
    RooDataSet* GetDataSet();
    /**
     * Get the category object
     */
    Category* GetCategoryObject();
  private:
    /**
     * Create the RooDataSet with the correct category variable
     */
    void MakeDataSet();
    /**
     * Settings for the fit
     */
    Settings m_Settings;
    /**
     * Original TTree with double tag events
     */
    TTree *m_Tree;
    /**
     * The fit variable
     */
    RooRealVar *m_SignalMBC;
    /**
     * Need this to apply a cut on the tag side
     */
    RooRealVar *m_TagMBC;
    /**
     * Dataset that we want to fit
     */
    RooDataSet *m_DataSet;
    /**
     * The object keeping track of all the categories
     */
    Category m_Category;
};

#endif
