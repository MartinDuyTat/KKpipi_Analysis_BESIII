// Martin Duy Tat 29th November 2021
/**
 * GeneratorKinematics is a simple struct containing all variables describing the generator kinematics
 */

#ifndef GENERATORKINEMATICS
#define GENERATORKINEMATICS

#include<vector>

struct GeneratorKinematics {
  GeneratorKinematics(int Particles = 100): ParticleIDs(Particles), MotherIndex(Particles), TruePx(Particles), TruePy(Particles), TruePz(Particles), TrueEnergy(Particles) {}
  /**
   * Number of particles generated in event
   */
  int NumberParticles;
  /**
   * List of particle IDs
   */
  std::vector<int> ParticleIDs;
  /**
   * List of mother indices
   */
  std::vector<int> MotherIndex;
  /**
   * List of generated momenta in the x-direction
   */
  std::vector<double> TruePx;
  /**
   * List of generated momenta in the y-direction
   */
  std::vector<double> TruePy;
  /**
   * List of generated momenta in the z-direction
   */
  std::vector<double> TruePz;
  /**
   * List of generated energies in the x-direction
   */
  std::vector<double> TrueEnergy;
};

#endif
