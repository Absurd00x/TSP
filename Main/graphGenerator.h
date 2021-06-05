#ifndef GRAPH_GENERATOR
#define GRAPH_GENERATOR

#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

#include "graph.h"

class GraphGenerator {
private:
  const long long INF = (1LL << 50);
  static const long long DEFAULT_EDGE_VALUE = std::numeric_limits<int>::max();
  static const unsigned int DEFAULT_VERTICES = 10u;
  static const unsigned int DEFAULT_EDGES = 30u;
  static const long long DEFAULT_MIN_WEIGHT = 10L;
  static const long long DEFAULT_MAX_WEIGHT = 99L;
  std::random_device rand_dev;
  std::mt19937 edgeRng;
  std::mt19937_64 weightRng;
  unsigned int vertices, edges;
  long long minWeight, maxWeight;

public:
  GraphGenerator(unsigned int vertices_=DEFAULT_VERTICES,
                 unsigned int edges_=DEFAULT_EDGES,
                 long long minWeight_=DEFAULT_MIN_WEIGHT,
                 long long maxWeight_=DEFAULT_MAX_WEIGHT)
                 :
                 edgeRng(rand_dev()),
                 weightRng(rand_dev()),
                 vertices(vertices_),
                 edges(edges_),
                 minWeight(minWeight_),
                 maxWeight(maxWeight_) {}

  void setEdgeSeed(long seed) {
    edgeRng = std::mt19937(seed);
  }

  void setWeightSeed(long seed) {
    weightRng = std::mt19937_64(seed);
  }

  Graph generateHamiltonianGraph() {
    ASSERT(edges >= vertices, "Not enough edges for graph to contain hamiltonian cycle");
    ASSERT(vertices > 2, "Degenerate case");
    ASSERT(edges <= vertices * (vertices - 1) / 2, "Too many edges");
    std::uniform_int_distribution<long long> distribution(minWeight, maxWeight);
    auto sampleWeight = [&distribution, this]() { return distribution(edgeRng); };

    // Сначала создается граф-кольцо, т.к. он является связным
    // и гаранированно содержит хотя бы один гамильтонов цикл
    // При этом из выборки будут удаляться уже выбранные ребра
    std::vector<std::vector<bool>> picked(vertices, std::vector<bool>(vertices, false));
    std::vector<unsigned int> permutation(vertices);
    std::iota(permutation.begin(), permutation.end(), 0);
    std::shuffle(permutation.begin(), permutation.end(), edgeRng);
    permutation.push_back(permutation.front());

    std::vector<Edge> edgesList;

    for(unsigned int i = 1; i < permutation.size(); ++i) {
      unsigned int u = permutation[i], v = permutation[i - 1];
      edgesList.emplace_back(sampleWeight(), v, u);
      picked[v][u] = picked[u][v] = true;
    }

    std::vector<Edge> additionalEdges;
    for(unsigned int i = 0; i < vertices; ++i) {
      for(unsigned int j = i + 1; j < vertices; ++j) {
        if (picked[i][j] == false)
          additionalEdges.emplace_back(sampleWeight(), i, j);
      }
    }
    shuffle(additionalEdges.begin(), additionalEdges.end(), edgeRng);
    unsigned int edgesLeft = edges - vertices;
    additionalEdges.resize(edgesLeft);
    edgesList.insert(edgesList.end(), additionalEdges.begin(), additionalEdges.end());

    // Рёбра отсортированы по весу
    std::sort(edgesList.begin(), edgesList.end());

    return Graph(vertices, edgesList, DEFAULT_EDGE_VALUE, true);
  }

  Graph generateMetricGraph() {
    ASSERT(vertices > 2, "Degenerate case");
    // Необходимо удовлетворить неравенство треугольника для весов рёбер
    // В этом случае отрезок возможных значений может будет равен
    // [minWeight + max(0, maxWeight - 2 * minWeight); maxWeight]
    long long shift = std::max(0LL, maxWeight - 2 * minWeight);
    std::uniform_int_distribution<long long> distribution(minWeight, maxWeight - shift);
    auto sampleWeight = [&distribution, this]() { return distribution(edgeRng); };

    std::vector<Edge> edgesList;
    for(unsigned int i = 0; i < vertices; ++i) {
      for(unsigned int j = i + 1; j < vertices; ++j) {
        edgesList.emplace_back(sampleWeight() + shift, i, j);
      }
    }

    // Рёбра отсортированы по весу
    std::sort(edgesList.begin(), edgesList.end());

    return Graph(vertices, edgesList, DEFAULT_EDGE_VALUE, true);
  }

};

#endif
