#ifndef TESTER_H
#define TESTER_H

#include "graph.h"
#include "graphGenerator.h"
#include "deterministicSolver.h"
#include "stochasticSolver.h"

class Tester {
private:
  const long long INF = (1LL << 50);
  Tester(){}
public:
  static std::vector<std::vector<long long>> test(
        GraphGenerator &gg,
        std::vector<DeterministicSolution> &dSolutions,
        std::map<StochasticSolution, std::map<std::string, long long>> &sSolutions,
        DeterministicSolver &dSolver,
        StochasticSolver &sSolver,
        long long tests,
        bool metric=true,
        bool verbose=true)
        {
    std::vector<std::vector<long long>> results(tests);
    long double percentSize = tests / 100.0L;
    long long curPercent = 0;
    auto start = clock();
    if (verbose == true)
      std::cout << "Started testing" << std::endl;
    Graph graph;
    for(long long finished = 0; finished < tests; ++finished) {
      if (verbose == true && (long long)(finished / percentSize) > curPercent) {
        curPercent = (long long)(finished / percentSize);
        std::cout << "Finished " << curPercent << "%\n";
        std::cout.flush();
      }
      graph = (metric ? gg.generateMetricGraph() : gg.generateHamiltonianGraph());
      for(auto method : dSolutions)
        results[finished].push_back(dSolver.solve(graph, method));
      for(auto [method, params] : sSolutions)
        results[finished].push_back(sSolver.solve(graph, method));
    }
    if (verbose == true) {
      std::cout << "Finished testing" << std::endl;
      auto finish = clock();
      std::cout << "Execution time: " << (long double)(finish - start) / CLOCKS_PER_SEC << 's' << std::endl;
    }
    return results;
  }
};

#endif
