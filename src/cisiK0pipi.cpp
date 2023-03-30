// Martin Duy Tat 13th February 2023

#include<stdexcept>
#include<vector>
#include<string>
#include"cisiK0pipi.h"
#include"Settings.h"

cisiK0pipi::cisiK0pipi(const Settings &settings):
  m_Ki_KSpipi(InitialiseKi(settings, "KSpipi")),
  m_Ki_KLpipi(InitialiseKi(settings, "KLpipi")),
  // Then initialise ci and si
  m_cisi(Initialisecisi(settings)) {
}

std::vector<std::pair<double, double>>
cisiK0pipi::InitialiseKi(const Settings &settings,
			 const std::string &TagMode) {
  if(TagMode != "KSpipi" && TagMode != "KLpipi") {
    throw std::invalid_argument(TagMode + " tag mode is not KSpipi or KLpipi");
  }
  std::vector<std::pair<double, double>> KiKbari(16);
  for(std::size_t i = 1; i <= 8; i++) {
    const std::string pName = TagMode + "_K_p" + std::to_string(i);
    const std::string mName = TagMode + "_K_m" + std::to_string(i);
    const std::string TagModeBinning = TagMode + "_BinningScheme";
    const double Ki = settings[TagModeBinning]["Ki"].getD(pName);
    const double Ki_err = settings[TagModeBinning]["Ki"].getD(pName + "_err");
    const double Kbari = settings[TagModeBinning]["Ki"].getD(mName);
    const double Kbari_err = settings[TagModeBinning]["Ki"].getD(mName + "_err");
    KiKbari[i - 1] = std::make_pair(Ki, Ki_err);
    KiKbari[i - 1 + 8] = std::make_pair(Kbari, Kbari_err);
  }
  return KiKbari;
}

std::vector<double> cisiK0pipi::Initialisecisi(const Settings &settings) {
  std::vector<double> cisi(32);
  const std::vector<std::string> TagModes{"KSpipi", "KLpipi"};
  const std::vector<std::string> c_or_s{"c", "s"};
  for(std::size_t i = 0; i < 2; i++) {
    for(std::size_t j = 0; j < 2; j++) {
      for(std::size_t Bin = 1; Bin <= 8; Bin++) {
	const std::string Name = TagModes[i] + "_" + c_or_s[j]
                               + std::to_string(Bin);
	const std::string TagBinning = TagModes[i] + "_BinningScheme";
	const double ci_or_si = settings[TagBinning]["cisi"].getD(Name);
	cisi[Bin - 1 + j*8 + i*16] = ci_or_si;
	//double cisi_err = settings[TagBinning]["cisi"].getD(Name + "_err");
      }
    }
  }
  /*auto CorrMatrix = ParseCorrelationMatrix(settings);
  TMatrixDSym CovMatrix(32);
  for(int i = 0; i < 32; i++) {
    for(int j = 0; j < 32; j++) {
      CovMatrix(i, j) = CorrMatrix(i, j)*cisi_Sigma[i]*cisi_Sigma[j];
    }
  }*/
  return cisi;
}

double cisiK0pipi::Get_Ki(std::size_t Bin, const std::string &TagMode) const {
  if(TagMode == "KLpipi") {
    return m_Ki_KLpipi[Bin - 1].first;
  } else {
    return m_Ki_KSpipi[Bin - 1].first;
  }
}

double cisiK0pipi::Get_Kbari(std::size_t Bin, const std::string &TagMode) const {
  if(TagMode == "KLpipi") {
    return m_Ki_KLpipi[Bin - 1 + 8].first;
  } else {
    return m_Ki_KSpipi[Bin - 1 + 8].first;
  }
}

double cisiK0pipi::Get_ci(std::size_t Bin, const std::string &TagMode) const {
  const std::size_t Index = Bin - 1 + (TagMode == "KSpipi" ? 0 : 16);
  return m_cisi[Index];
}

double cisiK0pipi::Get_si(std::size_t Bin, const std::string &TagMode) const {
  const std::size_t Index = Bin - 1 + 8 + (TagMode == "KSpipi" ? 0 : 16);
  return m_cisi[Index];
}

/*TMatrixT<double> cisiK0pipi::ParseCorrelationMatrix(const Settings &settings) const {
  TMatrixT<double> CorrMatrix(32, 32);
  std::ifstream CorrFile(settings.get("KSpipi_KLpipi_cisi_Correlation"));
  for(int i = 0; i < 32; i++) {
    std::string Line;
    std::getline(CorrFile, Line);
    std::stringstream ss(Line);
    for(int j = 0; j < 32; j++) {
      if(j >= i) {
	double Value;
	ss >> Value;
	Value /= 100.0;
	CorrMatrix(i, j) = Value;
      } else {
	CorrMatrix(i, j) = CorrMatrix(j, i);
      }
    }
  }
  CorrFile.close();
  return CorrMatrix;
}*/
