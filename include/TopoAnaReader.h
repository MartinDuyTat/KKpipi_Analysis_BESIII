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
#include<tuple>
#include<utility>

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
    void SaveAllComponentCuts() const;
    /**
     * Helper function that creates a .cut file with the appropriate cuts from a vector
     * @param Label Component label
     * @param DecayTopologies Vector of decay topology numbers from TopoAna
     */
    void SaveComponentCuts(const std::string &Label, const std::vector<int> &DecayTopologies) const;
  private:
    /**
     * A vector that a tuple that connects the label, for example "Signal" or "Other", to the decay descriptor and a vector of all the topology numbers from TopoAna for that particular component
     */
    std::vector<std::tuple<std::string, std::string, std::vector<int>>> m_DecayComponents;
    /**
     * Tuple with label for all other components not included in the list, by default it's "Other", and the topology numbers from TopoAna
     */
    std::pair<std::string, std::vector<int>> m_OtherComponent;
};

#endif
