// Martin Duy Tat 14th February 2022

#include<stdexcept>
#include<vector>
#include<utility>
#include"TH2F.h"

namespace DalitzUtilities {
  int LookUpBinNumber(double M2Plus, double M2Minus, const TH2F *BinningScheme) {
    int x = BinningScheme->GetXaxis()->FindBin(M2Plus);
    int y = BinningScheme->GetYaxis()->FindBin(M2Minus);
    Float_t BinNumberFloat = BinningScheme->GetBinContent(x, y);
    return M2Minus > M2Plus ? static_cast<int>(BinNumberFloat) : -static_cast<int>(BinNumberFloat);
  }

  int GetMappedK0hhBin(double M2Plus, double M2Minus, const TH2F *BinningScheme) {
    int x = BinningScheme->GetXaxis()->FindBin(M2Plus);
    int y = BinningScheme->GetYaxis()->FindBin(M2Minus);
    int i = 1;
    // Loop in circles around this point until reaching the Dalitz boundary
    while(true) {
      if(i > 1000) {
	throw std::runtime_error("Dalitz point too far outside phase space");
      }
      for(int j = 0; j <= i; j++) {
	// All possible combinations of displacements at the same distance
	std::vector<std::pair<int, int>> BinList{{i, j}, {i, -j}, {-i, j}, {-i, -j}, {j, i}, {j, -i}, {-j, i}, {-j, -i}};
	for(auto iter = BinList.begin(); iter != BinList.end(); iter++) {
	  int NewBin = static_cast<int>(BinningScheme->GetBinContent(x + iter->first, y + iter->second));
	  // Once we reach the Dalitz boundary the bin number is non-zero
	  if(NewBin != 0) {
	    return M2Minus > M2Plus ? NewBin : -NewBin;
	  }
	}
      }
      i++;
    }
  }

}
