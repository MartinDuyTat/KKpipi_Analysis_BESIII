// Martin Duy Tat 7th April 2021
/**
 * AnalyzeTopoAna takes in a list of particular decay components in a TopoAna file and finds all the decay topology numbers that match those decay components, the rest are classified as "Other" decay components
 * @param 1 Filename of TopoAna txt output file
 * @param 2 Filename with decay components
 */

#include<iostream>
#include"TopoAnaReader.h"

int main(int argc, char *argv[]) {
  if(argc != 3) {
    std::cout << "Need 2 input arguments\n";
    return 0;
  }
  std::cout << "Loading the decay components...\n";
  TopoAnaReader Reader{std::string(argv[2])};
  std::cout << "Decay components ready\n";
  std::cout << "Going through TopoAna output...\n";
  Reader.AnalyzeComponents(std::string(argv[1]));
  std::cout << "TopoAna analysis complete\n";
  std::cout << "Saving cut files...\n";
  Reader.SaveAllComponentCuts();
  std::cout << "Cuts safely stored in files\n";
  return 0;
}
