// Martin Duy Tat 4th October 2022
/**
 * BinnedDTYieldPrediction is an abstract class that stores the single tag yield and efficiency matrix of a tag, in addition to the flavour single tag yield, and it is used to predict the double tag yield in each bin
 */

#ifndef BINNEDDTYIELDPREDICTION
#define BINNEDDTYIELDPREDICTION

#include<vector>
#include<string>
#include"TMatrixT.h"
#include"Settings.h"

class BinnedDTYieldPrediction {
 public:
  /**
   * Constructor that stores the single tag yield, efficiency matrix and Ki
   * @param Tag The tag mode
   * @param Ki The Ki parameters normalised by the ST yield
   * @param Kbari The Kbari parameters normalised by the ST yield
   * @param settings The settings file
   */
  BinnedDTYieldPrediction(const std::string &Tag,
			  const std::vector<double> &Ki,
			  const std::vector<double> &Kbari,
			  const Settings &settings);
  /**
   * Function that returns the predicted bin yield
   */
  virtual std::vector<double> GetPredictedBinYields(
    double BF_KKpipi,
    const std::vector<double> &ci,
    const std::vector<double> &si) const = 0;
 protected:
  /**
   * The single tag yield of this tag, after efficiency correction
   */
  const double m_SingleTagYield;
  /**
   * The efficiency matrix for the CP even KKpipi model
   */
  const TMatrixT<double> m_EfficiencyMatrix_CPEven;
  /**
   * The efficiency matrix for the CP odd KKpipi model
   */
  const TMatrixT<double> m_EfficiencyMatrix_CPOdd;
  /**
   * The Ki parameters normalised by the ST yield
   */
  const std::vector<double> &m_Ki;
  /**
   * The Kbari parameters normalised by the ST yield
   */
  const std::vector<double> &m_Kbari;
 private:
  /**
   * Helper function to get the single tag yield
   * @param Tag The tag mode
   * @param settings The settings file
   */
  double GetSTYield(const std::string &Tag, const Settings &settings) const;
  /**
   * Helper function to get efficiency matrix
   * @param Tag The tag mode
   * @param settings The settings file
   * @param CPEvenOdd +1 for CP even model reweighted, -1 for CP odd model reweighted
   */
  TMatrixT<double> GetEffMatrix(const std::string &Tag,
				const Settings &settings,
				int CPEvenOdd) const;
};

#endif
