#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include "deterministicSolver.h"
#include "stochasticSolver.h"

struct Configuration {
  unsigned int vertices, edges, tests;
  long long minWeight, maxWeight;
  long edgesSeed, weightsSeed, solverSeed;
  bool metric;
  std::vector<DeterministicSolution> dSolutions;
  std::map<StochasticSolution, std::map<std::string, long long>> sSolutions;
  std::vector<long long> iterations;
  Configuration(){}
};

#endif
