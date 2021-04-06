// Martin Duy Tat 6th April 2021
/**
 * TopoAnaReader is a class for interpreting single tag component analysis of the output from TopoAna
 * The TopoAnaReader constructor will take in a file with a list of decay final states and look for these in the output txt file from TopoAna
 * The result is saved into .cut files so that these can be loaded directly into ROOT to pick out the components necessary for background study
 */

#ifndef TOPOANAREADER
#define TOTOANAREADER

#include<string>
#include<vector>

class TopoAnaReader {
  public:
    /**
     * Constructor that takes in a filename which contains the different decay components and a label
     * The format should be "Signal(pi+ pi- K+ K-)", where "Signal" is the label and "pi+ pi- K+ K-" is the decay, in the same order as shown in the TopoAna output
     * The last line should have the format "Other" without a decay, which is the label for all other decays, if this is not given all other decays are labelled as "Other"
     * @param Filename Filename with list of components to study
     */
    TopoAnaReader(const std::string &Filename);
    /**
     * Function that takes in the TopoAna txt output and sorts out the different components
     * @param Filename Filename of txt file from TopoAna
     */
    void AnalyzeComponents(const std::string &Filename);
    /**
     * Function that creates the .cut files that separates the different components
     * The filenames will be the label and a .cut extension
     */
    void SaveComponentCuts() const;
  private:
    /**
     * A vector with the different labels, not including the "Other" label
     */
    std::vector<std::string> ComponentLabels;
    /**
     * Label for all other components not included in the list
     */
    std::string OtherLabel;
    /**
     * Map connecting the label to a vector containing all the decay topology numbers from TopoAna for that particular component
     */
    std::map<std::string, std::vector<int>> DecayTopologyNumbers;
};

#endif
