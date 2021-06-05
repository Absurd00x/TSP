#ifndef GRAPH_H
#define GRAPH_H

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

typedef std::tuple<long long, unsigned int, unsigned int> Edge;

class Graph {
private:
  unsigned int nodes;
  std::vector<Edge> _edges;
  long long dummyValue;
  bool bidirectional;
  std::vector<std::vector<long long>> matrix;

public:
  Graph(unsigned int n,
        const std::vector<Edge> &edges,
        long long dummyValue_,
        bool bidirectional_=true)
        :
        nodes(n),
        _edges(edges.begin(), edges.end()),
        dummyValue(dummyValue_),
        bidirectional(bidirectional_)
        {
    matrix.resize(n);
    for(auto &row : matrix)
      row.resize(n, dummyValue);

    for(const auto& [weight, from, to] : edges) {
      ASSERT(from < n && to < n, "Edge out of range");
      matrix[from][to] = weight;
      if (bidirectional)
        matrix[to][from] = weight;
    }
  }

  Graph() :
      nodes(0),
      _edges(),
      dummyValue(std::numeric_limits<int>::max()),
      bidirectional(true){}

  Graph(Graph &other) :
        nodes(other.getNodes()),
        _edges(other.getAllEdges()),
        dummyValue(other.getDummyValue()),
        bidirectional(other.isBidirectional()){}

  Graph& operator=(Graph other) {
    nodes = other.nodes;
    matrix = other.matrix;
    _edges = other._edges;
    dummyValue = other.dummyValue;
    bidirectional = other.bidirectional;
    return (*this);
  }


  std::vector<Edge> getAllEdges() {
    return _edges;
  }

  bool isBidirectional() {
    return bidirectional;
  }

  long long traverse(const std::vector<unsigned int> &path) const {
    long long result = 0;
    for(unsigned int i = 1; i < path.size(); ++i)
      result += matrix[path[i - 1]][path[i]];
    return result;
  }

  long long tryTraverse(const std::vector<unsigned int> &path, bool &valid) const {
    long long result = 0;
    for(unsigned int i = 1; i < path.size(); ++i) {
      long long weight = matrix[path[i - 1]][path[i]];
      if (weight == dummyValue) {
        valid = false;
        return result;
      }
      result += weight;
    }
    return result;
  }

  long long getDummyValue() {
    return dummyValue;
  }

  long long getEdge(unsigned int from, unsigned int to) {
    return matrix[from][to];
  }

  bool hasEdge(unsigned int from, unsigned int to) {
    return matrix[from][to] != dummyValue;
  }

  int getNodes() {
    return nodes;
  }

  void print() {
    std::cout << "Граф из " << nodes << " вершин и " << _edges.size() << " ребер:\n";
    for(auto &[weight, from, to] : _edges)
      std::cout << from + 1 << ' ' << to + 1 << ' ' << weight << '\n';
    std::cout.flush();
  }
};

#endif
