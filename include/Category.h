// Martin Duy Tat 25th November 2021
/**
 * Category is a class that stores the tag mode and the type of tag (Flavour, CP, SCMB), and from this it can create a unique string that labels the tag mode and the bin number that each event belongs to
 */

#ifndef CATEGORY
#define CATEGORY

#include<string>
#include<vector>
#include"RooCategory.h"
#include"Settings.h"

class Category {
  public:
    /**
     * Constructor that reads and stores the tag mode and tag type from the settings file
     * @param settings Settings object
     */
    Category(const Settings &settings);
    /**
     * Get the unique string that identifies the phase space bin of a tag mode
     * @param Bin Bin number
     */
    std::string GetCategory(int SignalBin, int TagBin = 0) const;
    /**
     * Same as GetCateory()
     */
    std::string operator ()(int SignalBin, int TagBin = 0) const;
    /**
     * Get a vector of all the categories
     */
    std::vector<std::string> GetCategories() const;
    /**
     * Get the category variable
     */
    RooCategory* GetCategoryVariable();
    /**
     * Get the bin number from the unique string that describes the category
     */
    int GetSignalBinNumber(const std::string &category) const;
  private:
    /**
     * Tag mode
     * It is assumed that the signal mode is KKpipi, this string only identifies the tag itself
     */
    std::string m_TagMode;
    /**
     * Type of tag mode
     * The number of bins depends on what type of tag it is (these numbers are for 2x8 KKpipi bins)
     * CP tags: 8 bins on the signal side
     * Flavour tags: 16 bins on the signal side
     * SCMB with KKpipi: 64 bins, but only 36 are unique
     * SCMB with KSpipi: 16x8 bins
     * SCMB with KSKK: 16x2 bins
     */
    std::string m_Type;
    /**
     * Number of bins for KKpipi
     */
    int m_SignalBins;
    /**
     * Number of bins on the tag side
     */
    int m_TagBins;
    /**
     * Category variable used in the fit
     */
    RooCategory m_CategoryVar;
    /**
     * Helper function that checks whether or not the bins are valid and throw an appropriate exception if not
     */
    void CheckValidBins(int SignalBin, int TagBin) const;
};

#endif
