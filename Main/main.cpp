#include <algorithm>
#include <chrono>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "dataManager.h"
#include "tester.h"

int main() {
  const bool VERBOSE = true;
  const bool OVERWRITE = false;
  Configuration conf = DataManager::readConfiguration();
  const std::string FILENAME = "sample";
  DeterministicSolver dSolver;
  StochasticSolver sSolver(conf.sSolutions);

  sSolver.setSeed(conf.solverSeed);

  GraphGenerator gg = GraphGenerator(conf.vertices, conf.edges, conf.minWeight, conf.maxWeight);
  gg.setEdgeSeed(conf.edgesSeed);
  gg.setWeightSeed(conf.weightsSeed);

  auto data = Tester::test(gg, conf.dSolutions, conf.sSolutions, dSolver, sSolver,
                           conf.tests, conf.metric, VERBOSE);
  DataManager::saveData(data, conf, FILENAME, OVERWRITE);

  getchar();

  return 0;
}
