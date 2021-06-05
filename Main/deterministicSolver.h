#ifndef SOLVER_H
#define SOLVER_H

#include "graph.h"

enum DeterministicSolution {BRUTE=0, DYNAMIC_MIN=1, DYNAMIC_MAX=2, NEAREST_NEIGHBOUR=3, MULTIFRAGMENT=4, CHRISTOFIDES=5, UNKNOWN_DETERMINISTIC=6};

class DeterministicSolver {
private:
  const long long INF = (1LL << 50);
  inline static const std::string BRUTE_NAME = "BRUTE_FORCE";
  inline static const std::string DYNAMIC_MIN_NAME = "DYNAMIC_PROGRAMMING_MIN";
  inline static const std::string DYNAMIC_MAX_NAME = "DYNAMIC_PROGRAMMING_MAX";
  inline static const std::string CHRISTOFIDES_NAME = "CHRISTOFIDES";
  inline static const std::string MULTIFRAGMENT_NAME = "MULTIFRAGMENT";
  inline static const std::string NEIGHBOUR_NAME = "NEAREST_NEIGHBOUR";
  inline static const std::string UNKNOWN_DETERMINISTIC = "UNKNOWN_DETERMINISTIC";

  std::vector<unsigned int> solveBrute(Graph &graph) {
    std::vector<unsigned int> permutation(graph.getNodes());
    std::iota(permutation.begin(), permutation.end(), 0u);
    long long best = std::numeric_limits<long long>::max();
    std::vector<unsigned int> result;
    do {
      bool valid = true;
      long long cost = graph.tryTraverse(permutation, valid);
      if (valid == true && cost < best) {
        best = cost;
        result = permutation;
      }
    } while(std::next_permutation(permutation.begin(), permutation.end()));
    return result;
  }

  std::vector<unsigned int> solveDynamicMin(Graph &graph) {
    const int N = graph.getNodes();
    const long long MASKS = (1LL << N);
    const int DUMMY_VERTEX = -1;
    const long long FINISH_MASK = MASKS - 1;
    std::vector<std::vector<long long>> dp(1LL << N, std::vector<long long>(N, INF));
    std::vector<std::vector<int>> prev(1LL << N, std::vector<int>(N, DUMMY_VERTEX));

    for(int i = 0; i < N; ++i)
      dp[1LL << i][i] = 0;

    for(long long mask = 1; mask < MASKS; ++mask) {
      for(int bitMain = 0; bitMain < N; ++bitMain) {
        if (((mask >> bitMain) & 1) == 0)
          continue;
        long long subMask = mask ^ (1LL << bitMain);
        for(int bitSub = 0; bitSub < N; ++bitSub) {
          if (((subMask >> bitSub) & 1) == 0)
            continue;
          if (graph.hasEdge(bitSub, bitMain) == false)
            continue;
          long long weight = graph.getEdge(bitSub, bitMain);
          if (dp[mask][bitMain] > dp[subMask][bitSub] + weight) {
            dp[mask][bitMain] = dp[subMask][bitSub] + weight;
            prev[mask][bitMain] = bitSub;
          }
        }
      }
    }
    auto min = std::min_element(dp[FINISH_MASK].begin(), dp[FINISH_MASK].end());
    int minIndex = (int)(min - dp[FINISH_MASK].begin());
    std::vector<unsigned int> result = {(unsigned int)minIndex};

    long long mask = FINISH_MASK;
    int cur = minIndex, next = prev[FINISH_MASK][minIndex];
    for(; next != DUMMY_VERTEX; cur = next, next = prev[mask][cur]) {
      result.push_back(next);
      mask ^= (1LL << cur);
    }
    return result;
  }

  std::vector<unsigned int> solveDynamicMax(Graph &graph) {
    const int N = graph.getNodes();
    const long long MASKS = (1LL << N);
    const int DUMMY_VERTEX = -1;
    const long long FINISH_MASK = MASKS - 1;
    std::vector<std::vector<long long>> dp(1LL << N, std::vector<long long>(N, -INF));
    std::vector<std::vector<int>> prev(1LL << N, std::vector<int>(N, DUMMY_VERTEX));

    for(int i = 0; i < N; ++i)
      dp[1LL << i][i] = 0;

    for(long long mask = 1; mask < MASKS; ++mask) {
      for(int bitMain = 0; bitMain < N; ++bitMain) {
        if (((mask >> bitMain) & 1) == 0)
          continue;
        long long subMask = mask ^ (1LL << bitMain);
        for(int bitSub = 0; bitSub < N; ++bitSub) {
          if (((subMask >> bitSub) & 1) == 0)
            continue;
          if (graph.hasEdge(bitSub, bitMain) == false)
            continue;
          long long weight = graph.getEdge(bitSub, bitMain);
          if (dp[mask][bitMain] < dp[subMask][bitSub] + weight) {
            dp[mask][bitMain] = dp[subMask][bitSub] + weight;
            prev[mask][bitMain] = bitSub;
          }
        }
      }
    }
    auto max = std::max_element(dp[FINISH_MASK].begin(), dp[FINISH_MASK].end());
    int maxIndex = (int)(max - dp[FINISH_MASK].begin());
    std::vector<unsigned int> result = {(unsigned int)maxIndex};

    long long mask = FINISH_MASK;
    int cur = maxIndex, next = prev[FINISH_MASK][maxIndex];
    for(; next != DUMMY_VERTEX; cur = next, next = prev[mask][cur]) {
      result.push_back(next);
      mask ^= (1LL << cur);
    }
    return result;
  }

  std::vector<unsigned int> solveNearestNeighbor(Graph &graph) {
    const unsigned int N = graph.getNodes();

    std::vector<Edge> edges = graph.getAllEdges();
    std::vector<std::vector<std::pair<long long, unsigned int>>> denseGraph(N);
    for(auto &[weight, from, to] : edges) {
      denseGraph[from].emplace_back(weight, to);
      denseGraph[to].emplace_back(weight, from);
    }
    std::function<void(int, std::vector<bool>&, std::vector<unsigned int>&)> dfs
            = [&dfs, &denseGraph, N](unsigned int cur, std::vector<bool> &visited, std::vector<unsigned int> &tour) {
      visited[cur] = true;
      tour.push_back(cur);
      if (N == tour.size())
        return;

      for(auto &[weight, to] : denseGraph[cur]) {
        if (visited[to] == true)
          continue;
        dfs(to, visited, tour);
        if (N == tour.size())
          return;
      }
      visited[cur] = false;
      tour.pop_back();
    };
    std::vector<bool> visited(N, false);
    std::vector<unsigned int> result;
    long long best = INF;
    for(unsigned int start = 0; start < N; ++start) {
      std::fill(visited.begin(), visited.end(), false);
      std::vector<unsigned int> tour;
      dfs(start, visited, tour);
      long long cost = graph.traverse(tour);
      if (best > cost) {
        best = cost;
        std::swap(result, tour);
      }
    }
    return result;
  }

  std::vector<unsigned int> solveMultiFragment(Graph &graph) {
    const unsigned int N = graph.getNodes();
    std::vector<int> leader(N);
    std::iota(leader.begin(), leader.end(), 0);
    std::vector<int> rank(N, 1);

    std::function<int(int)> findSet = [&leader, &findSet](int v) {
      if (v == leader[v])
        return v;
      return leader[v] = findSet(leader[v]);
    };

    auto unionSets = [&](unsigned int a, unsigned int b) {
      a = findSet(a);
      b = findSet(b);
      if (a != b) {
        if (rank[a] < rank[b])
          std::swap(a, b);
        leader[b] = a;
        if (rank[a] == rank[b])
          ++rank[a];
        return true;
      }
      return false;
    };

    std::vector<int> degrees(N, 0);
    std::vector<std::vector<int>> tourGraph(N);
    long long cost = 0;
    for(auto &[weight, from, to] : graph.getAllEdges()) {
      if (degrees[from] > 1 || degrees[to] > 1)
        continue;
      bool unite = unionSets(from, to);
      if (unite) {
        cost += graph.getEdge(from, to);
        ++degrees[from], ++degrees[to];
        tourGraph[from].push_back(to);
        tourGraph[to].push_back(from);
      }
    }
    int start = (int)(std::find(degrees.begin(), degrees.end(), 1) - degrees.begin());
    std::vector<unsigned int> result;
    std::function<void(int, int)> dfs = [&result, &tourGraph, &dfs](int cur, int prev) {
      result.push_back(cur);
      for(auto next : tourGraph[cur]) {
        if (next != prev) {
          dfs(next, cur);
        }
      }
    };
    dfs(start, -1);
    return result;
  }

  std::vector<unsigned int> solveCristofides(Graph &graph) {
    const unsigned int N = graph.getNodes();
    std::vector<int> leader(N);
    std::iota(leader.begin(), leader.end(), 0);
    std::vector<int> rank(N, 1);

    std::function<int(int)> findSet = [&leader, &findSet](int v) {
      if (v == leader[v])
        return v;
      return leader[v] = findSet(leader[v]);
    };

    auto unionSets = [&](unsigned int a, unsigned int b) {
      a = findSet(a);
      b = findSet(b);
      if (a != b) {
        if (rank[a] < rank[b])
          std::swap(a, b);
        leader[b] = a;
        if (rank[a] == rank[b])
          ++rank[a];
        return true;
      }
      return false;
    };

    std::vector<int> degrees(N, 0);
    std::vector<Edge> edgesMST;
    for(auto &[weight, from, to] : graph.getAllEdges()) {
      bool unite = unionSets(from, to);
      if (unite) {
        ++degrees[from], ++degrees[to];
        edgesMST.emplace_back(weight, from, to);
      }
    }

    std::vector<unsigned int> picked;
    unsigned int pickedVetices = 0;
    for (unsigned int vertex = 0; vertex < N; ++vertex) {
      if (degrees[vertex] % 2 == 1) {
        ++pickedVetices;
        picked.push_back(vertex);
      }
    }

    std::vector<Edge> edgesMatching;
    for(unsigned int i = 0; i < picked.size(); ++i) {
      for(unsigned int j = i + 1; j < picked.size(); ++j) {
        // 1-индексация
        edgesMatching.emplace_back(graph.getEdge(picked[i], picked[j]), i + 1, j + 1);
      }
    }

    Graph graphMatching(pickedVetices + 1, edgesMatching, INF, true);
    std::vector<long long> v(pickedVetices + 1), u(pickedVetices + 1);
    std::vector<unsigned int> p(pickedVetices + 1), way(pickedVetices + 1);
    for(unsigned int i = 1; i <= pickedVetices; ++i) {
      p[0] = i;
      unsigned int j0 = 0;
      std::vector<long long> minv(pickedVetices + 1, INF);
      std::vector<char> used(pickedVetices + 1, false);
      do {
        used[j0] = true;
        unsigned int i0 = p[j0], j1;
        long long delta = INF;
        for(unsigned int j = 1; j <= pickedVetices; ++j)
          if(!used[j]) {
            long long cur = graphMatching.getEdge(i0, j) - u[i0] - v[j];
            if (cur < minv[j])
              minv[j] = cur, way[j] = j0;
            if (minv[j] < delta)
              delta = minv[j], j1 = j;
          }
        for(unsigned int j = 0; j <= pickedVetices; ++j)
          if (used[j])
            u[p[j]] += delta, v[j] -= delta;
          else
            minv[j] -= delta;
        j0 = j1;
      } while(p[j0] != 0);
      do {
        unsigned int j1 = way[j0];
        p[j0] = p[j1];
        j0 = j1;
      } while(j0);
    }

    std::vector<Edge> edgesFound;
    for(unsigned int j = 1; j <= pickedVetices; ++j) {
      edgesFound.emplace_back(graphMatching.getEdge(j, p[j]), picked[j - 1], picked[p[j] - 1]);
    }

    edgesFound.insert(edgesFound.end(), edgesMST.begin(), edgesMST.end());
    std::sort(edgesFound.begin(), edgesFound.end());
    std::vector<std::vector<std::tuple<long long, unsigned int, unsigned int>>> multigraph(N);
    unsigned int counter = 0;
    std::vector<bool> used(edgesMST.size() + edgesFound.size(), false);
    for(auto &[weight, from, to] : edgesFound) {
      multigraph[from].emplace_back(weight, to, counter);
      multigraph[to].emplace_back(weight, from, counter);
      ++counter;
    }

    std::vector<unsigned int> eulerTour;
    std::function<void(int)>dfs = [&eulerTour, &multigraph, &used, &dfs](unsigned int cur) {
      for(auto &[weight, to, index] : multigraph[cur]) {
        if (!used[index]) {
          used[index] = true;
          dfs(to);
        }
      }
      eulerTour.push_back(cur);
    };
    dfs(0);
    std::vector<unsigned int> result;
    std::set<unsigned int> set;
    for(auto cur : eulerTour) {
      if (set.count(cur) == 0) {
        set.insert(cur);
        result.push_back(cur);
      }
    }
    return result;
  }

public:
  std::vector<unsigned int> getCertificate(Graph &graph, DeterministicSolution solution) {
    std::vector<unsigned int> result;
    switch (solution) {
      case BRUTE:             result = solveBrute(graph); break;
      case DYNAMIC_MIN:       result = solveDynamicMin(graph); break;
      case DYNAMIC_MAX:       result = solveDynamicMax(graph); break;
      case CHRISTOFIDES:      result = solveCristofides(graph); break;
      case MULTIFRAGMENT:     result = solveMultiFragment(graph); break;
      case NEAREST_NEIGHBOUR: result = solveNearestNeighbor(graph); break;
      default:                break;
    }
    return result;
  }

  long long solve(Graph &graph, DeterministicSolution solution) {
    std::vector<unsigned int> result = getCertificate(graph, solution);
    return graph.traverse(result);
  }

  static std::vector<DeterministicSolution> getAllSolutions() {
    std::vector<DeterministicSolution> solutions = {
      DeterministicSolution::DYNAMIC_MIN,
      DeterministicSolution::DYNAMIC_MAX,
      DeterministicSolution::MULTIFRAGMENT,
      DeterministicSolution::NEAREST_NEIGHBOUR,
      DeterministicSolution::CHRISTOFIDES
    };
    return solutions;
  }

  static std::string getName(DeterministicSolution solution) {
    switch (solution) {
      case BRUTE:             return BRUTE_NAME;
      case DYNAMIC_MIN:       return DYNAMIC_MIN_NAME;
      case DYNAMIC_MAX:       return DYNAMIC_MAX_NAME;
      case CHRISTOFIDES:      return CHRISTOFIDES_NAME;
      case MULTIFRAGMENT:     return MULTIFRAGMENT_NAME;
      case NEAREST_NEIGHBOUR: return NEIGHBOUR_NAME;
      default:                return UNKNOWN_DETERMINISTIC;
    }
  }

  static DeterministicSolution getId(std::string name) {
    if (name == BRUTE_NAME)
      return DeterministicSolution::BRUTE;
    else if (name == DYNAMIC_MIN_NAME)
      return DeterministicSolution::DYNAMIC_MIN;
    else if (name == DYNAMIC_MAX_NAME)
      return DeterministicSolution::DYNAMIC_MAX;
    else if (name == CHRISTOFIDES_NAME)
      return DeterministicSolution::CHRISTOFIDES;
    else if (name == MULTIFRAGMENT_NAME)
      return DeterministicSolution::MULTIFRAGMENT;
    else if (name == NEIGHBOUR_NAME)
      return DeterministicSolution::NEAREST_NEIGHBOUR;
    return DeterministicSolution::UNKNOWN_DETERMINISTIC;
  }
};

#endif
