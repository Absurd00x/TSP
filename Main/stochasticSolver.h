#ifndef STOCHASTIC_H
#define STOCHASTIC_H

#include "graph.h"

enum StochasticSolution {MONTE_CARLO=0, ANNEALING=1, GENETIC=2, UNKNOWN_STOCHASTIC=3};

class StochasticSolver {
private:
  const long long INF = (1LL << 50);
  const long long DEFAULT_ITERATIONS = 1000LL;
  inline static const std::string ANNEALING_NAME = "ANNEALING";
  inline static const std::string MONTE_CARLO_NAME = "MONTE_CARLO";
  inline static const std::string GENETIC_NAME = "GENETIC";
  inline static const std::string UNKNOWN_NAME = "UNKNOWN_STOCHASTIC";
  inline static const std::string KEY_ITERATIONS = "iters";
  inline static const std::string KEY_MUTATIONS = "mutations";
  inline static const std::string KEY_POPULATION = "population";
  std::random_device rand_dev;
  std::mt19937 rng;
  std::map<StochasticSolution, std::map<std::string, long long>> params;
  const int ALGOS_IMPLEMENTED = 3;

  std::vector<unsigned int> solveMonteCarlo(Graph &graph) {
    const int N = graph.getNodes();
    const int iterations = (int)params[StochasticSolution::MONTE_CARLO][KEY_ITERATIONS];
    std::vector<unsigned int> result, perm(N);
    std::iota(perm.begin(), perm.end(), 0);
    long long best = INF;
    for(int passed = 0; passed < iterations; ++passed) {
      std::shuffle(perm.begin(), perm.end(), rng);
      bool valid = true;
      long long cost = graph.tryTraverse(perm, valid);
      if (valid && best > cost) {
        best = cost;
        result = perm;
      }
    }
    return result;
  }

  std::vector<unsigned int> solveAnnealing(Graph &graph) {
    const int N = graph.getNodes();
    const int iterations = (int)params[StochasticSolution::ANNEALING][KEY_ITERATIONS];
    std::vector<unsigned int> result(N);
    std::iota(result.begin(), result.end(), 0);
    std::shuffle(result.begin(), result.end(), rng);

    std::uniform_int_distribution<int> neighboursPool(0, N - 1);
    auto sampleNeighbour = [&neighboursPool, this]() {return neighboursPool(rng);};
    std::uniform_real_distribution<double> probabilityPool(0.0, 1.0);
    auto sampleProbability = [&probabilityPool, this]() {return (long double)probabilityPool(rng);};

    long double temperature = 1.0L;
    auto tryTransit = [&temperature, &sampleProbability]
                      (long long currentEnergy, long long neighbourEnergy) {
      if (currentEnergy < neighbourEnergy)
        return false;
      long double P = expl(-(neighbourEnergy - currentEnergy) / temperature);
      long double random = sampleProbability();
      return P > random;
    };

    auto getChanges = [&result, &graph, N](int first, int second, int sign=1) {
      std::vector<int> switches = {first, second};
      long long res = 0;
      for(int permIndex : switches) {
        if (permIndex > 0)
          res += graph.getEdge(result[permIndex], result[permIndex - 1]) * sign;
        if (permIndex < N - 1)
          res += graph.getEdge(result[permIndex], result[permIndex + 1]) * sign;
      }
      return res;
    };

    const long double coolingFactor = 1.0L / powl(10.0L, (6.0L / iterations));
    std::deque<int> neighbourQueue;
    long long best = graph.traverse(result);
    for(int passed = 0; passed < iterations; ++passed) {
      temperature *= coolingFactor;

      while (neighbourQueue.size() < 2)
        neighbourQueue.push_back(sampleNeighbour());
      while (neighbourQueue.front() == neighbourQueue.back())
        neighbourQueue.push_back(sampleNeighbour());
      int first = neighbourQueue.front(), second = neighbourQueue.back();

      neighbourQueue.pop_front(), neighbourQueue.pop_back();

      long long extractedCost = best + getChanges(first, second, -1);
      std::swap(result[first], result[second]);
      long long nextCost = extractedCost + getChanges(first, second, 1);

      if (!tryTransit(best, nextCost))
        std::swap(result[first], result[second]);
      else
        best = nextCost;
    }
    bool valid = true;
    graph.tryTraverse(result, valid);
    if (!valid)
      result.clear();
    return result;
  }

  std::vector<unsigned int> solveGenetic(Graph &graph) {
    const int N = graph.getNodes();
    const int iterations = (int)params[StochasticSolution::GENETIC][KEY_ITERATIONS];
    const bool mutations = (bool)params[StochasticSolution::GENETIC][KEY_MUTATIONS];
    const int populationSize = (int)params[StochasticSolution::GENETIC][KEY_POPULATION];

    std::vector<unsigned int> result;
    std::vector<std::pair<long long, std::vector<unsigned int>>> population(populationSize);
    for(int i = 0; i < populationSize; ++i) {
      std::vector<unsigned int> path(N);
      std::iota(path.begin(), path.end(), 0u);
      std::shuffle(path.begin(), path.end(), rng);
      long long cost = graph.traverse(path);
      population[i].first = cost;
      std::swap(population[i].second, path);
    }

    auto comp = [](const std::pair<long long, std::vector<unsigned int>> &a,
                   const std::pair<long long, std::vector<unsigned int>> &b) {
      return a.first < b.first;
    };
    std::sort(population.begin(), population.end(), comp);

    auto willSurvive = [populationSize, this](int rank) {
      std::uniform_int_distribution<int> getProb(rank, populationSize);
      int value = getProb(rng);
      return value >= rank;
    };

    auto breed = [N, this](std::vector<unsigned int> &one, std::vector<unsigned int> &other) {
      std::uniform_int_distribution<int> randIndex(0, N / 3 - 1);
      int ind1 = randIndex(rng);
      int ind2 = randIndex(rng);
      std::vector<unsigned int> child = one;
      for(int i = 0; i < N / 2; ++i)
        child[ind1 + i] = other[ind2 + i];
      std::vector<int> count(N);
      for(int i = 0; i < N; ++i)
        ++count[child[i]];
      int zeroIndex = 0;
      for(int i = 0; i < N; ++i) {
        if (count[child[i]] > 1) {
          while(count[zeroIndex] != 0)
            ++zeroIndex;
          --count[child[i]];
          child[i] = zeroIndex;
          ++count[child[i]];
        }
      }
      return child;
    };

    long long best = population.front().first;
    result = population.front().second;
    for(int passed = 0; passed < iterations; ++passed) {
      std::vector<std::pair<long long, std::vector<unsigned int>>> nextGeneration;
      // Естественный отбор
      for(int i = 0; i < populationSize; ++i) {
        if (willSurvive(i)) {
          auto &[value, genom] = population[i];
          nextGeneration.push_back({});
          nextGeneration.back().first = value;
          std::swap(nextGeneration.back().second, genom);
        }
      }
      // Скрещивание
      int parents = (int)nextGeneration.size();
      int need = populationSize - parents;
      std::uniform_int_distribution<int> randomIndex(0, parents - 1);
      for(int i = 0; i < need; ++i) {
        int parent1 = i;
        int parent2 = randomIndex(rng);
        std::pair<long long, std::vector<unsigned int>> child;
        child.second = breed(nextGeneration[parent1].second, nextGeneration[parent2].second);
      }
      // Мутации
      if (mutations) {
        randomIndex = std::uniform_int_distribution(0, N - 1);
        for(int i = 0; i < populationSize; ++i) {
          int index1 = randomIndex(rng);
          int index2 = randomIndex(rng);
          std::swap(nextGeneration[i].second[index1], nextGeneration[i].second[index2]);
          nextGeneration[i].first = graph.traverse(nextGeneration[i].second);
        }
      }
      std::swap(population, nextGeneration);
      std::sort(population.begin(), population.end(), comp);
      if (best < population.front().first) {
        best = population.front().first;
        result = population.front().second;
      }
    }

    return result;
  }

public:
  static const long long INVALID_RESULT = std::numeric_limits<long long>::max();

  StochasticSolver(const std::map<StochasticSolution, std::map<std::string, long long>> &params_) : rng(rand_dev()) {
    params = params_;
    for(auto [method, conf] : params) {
      if (params[method].count(KEY_ITERATIONS) == 0) {
        params[method][KEY_ITERATIONS] = DEFAULT_ITERATIONS;
      }
    }
    if (params[StochasticSolution::GENETIC].count(KEY_MUTATIONS) == 0)
      params[StochasticSolution::GENETIC][KEY_MUTATIONS] = 1;
    if (params[StochasticSolution::GENETIC].count(KEY_POPULATION) == 0)
      params[StochasticSolution::GENETIC][KEY_POPULATION] = 1;
  }

  std::vector<unsigned int> getCertificate(Graph &graph, StochasticSolution solution) {
    std::vector<unsigned int> result;
    switch (solution) {
      case ANNEALING:   result = solveAnnealing(graph); break;
      case MONTE_CARLO: result = solveMonteCarlo(graph); break;
      case GENETIC:     result = solveGenetic(graph); break;
      default:          break;
    }
    return result;
  }

  long long solve(Graph &graph, StochasticSolution solution) {
    std::vector<unsigned int> result = getCertificate(graph, solution);
    if (result.empty() == true)
      return INVALID_RESULT;
    return graph.traverse(result);
  }

  static std::string getName(StochasticSolution solution) {
    switch (solution) {
      case ANNEALING:   return ANNEALING_NAME;
      case MONTE_CARLO: return MONTE_CARLO_NAME;
      case GENETIC:     return GENETIC_NAME;
      default:          break;
    }
    return UNKNOWN_NAME;
  }

  static StochasticSolution getId(std::string name) {
    if (name == ANNEALING_NAME)
      return StochasticSolution::ANNEALING;
    else if (name == MONTE_CARLO_NAME)
      return StochasticSolution::MONTE_CARLO;
    else if (name == GENETIC_NAME)
      return StochasticSolution::GENETIC;
    return StochasticSolution::UNKNOWN_STOCHASTIC;
  }

  static std::vector<StochasticSolution> getAllSolutions() {
    std::vector<StochasticSolution> solutions = {
      StochasticSolution::ANNEALING,
      StochasticSolution::MONTE_CARLO,
      StochasticSolution::GENETIC
    };
    return solutions;
  }

  void setSeed(long seed) {
    rng = std::mt19937(seed);
  }

};

#endif
