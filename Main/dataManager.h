#ifndef DATA_MANAGER_H
#define DATA_MANAGER_H

#include "configuration.h"
#include "deterministicSolver.h"
#include "stochasticSolver.h"

class DataManager {
private:
  DataManager(){}
  inline static const std::string dataFolder = "../Data/";
  inline static const std::string dataFormat = ".csv";
  inline static const std::string confFormat = ".conf";
public:

  static void saveData(
            std::vector<std::vector<long long>>data,
            Configuration conf,
            std::string filename,
            bool overwrite=false) {
    std::filesystem::path dataPath{dataFolder + filename + dataFormat};
    std::filesystem::path confPath{dataFolder + filename + confFormat};
    if (std::filesystem::exists(dataPath) == true
     && overwrite == false) {
       bool finished = false;
       for(int label = 2; finished == false; ++label) {
         std::string suffix = "(" + std::to_string(label) + ")";
         std::filesystem::path newPath{dataFolder + filename + suffix + dataFormat};
         if (std::filesystem::exists(newPath) == false) {
           dataPath = newPath;
           confPath = std::filesystem::path{dataFolder + filename + suffix + confFormat};
           finished = true;
         }
      }
    }
    std::ofstream fout(dataPath);
    fout.tie(0);
    for(unsigned int i = 0; i < conf.dSolutions.size(); ++i) {
      fout << DeterministicSolver::getName(conf.dSolutions[i]);
      fout << (i + 1 == conf.dSolutions.size() ? "" : ",");
    }
    for (auto &[method, params] : conf.sSolutions) {
      fout << ",";
      fout << StochasticSolver::getName(method);
    }
    fout << '\n';
    for(auto &row : data) {
      for(unsigned int i = 0; i < row.size(); ++i) {
        fout << row[i];
        fout << (i + 1 == row.size() ? "\n" : ",");
      }
    }
    fout.flush();
    fout.close();
    fout = std::ofstream(confPath);
    fout << conf.vertices << "\n";
    fout << conf.edges << "\n";
    fout << conf.minWeight << "\n";
    fout << conf.maxWeight << "\n";
    fout << (int)(conf.metric) << "\n";
    fout << conf.tests << "\n";
    fout << conf.edgesSeed << "\n";
    fout << conf.weightsSeed << "\n";
    fout << conf.solverSeed << "\n";
    fout << conf.dSolutions.size() << "\n";
    for(auto dSolution : conf.dSolutions)
      fout << DeterministicSolver::getName(dSolution) << "\n";
    fout << conf.sSolutions.size() << "\n";
    for(auto &[sSolution, params] : conf.sSolutions) {
      fout << StochasticSolver::getName(sSolution) << "\n";
      fout << params.size() << "\n";
      for(auto &[key, value] : params)
        fout << key << ' ' << value << "\n";
    }
    fout << "\n";
    fout.flush();
    fout.close();
  }

  static Configuration readConfiguration(std::string filename="default") {
    std::ifstream fin(dataFolder + filename + confFormat);
    Configuration conf;
    fin >> conf.vertices;
    fin >> conf.edges;
    fin >> conf.minWeight;
    fin >> conf.maxWeight;
    fin >> conf.metric;
    fin >> conf.tests;
    fin >> conf.edgesSeed;
    fin >> conf.weightsSeed;
    fin >> conf.solverSeed;
    unsigned int sz;
    fin >> sz;
    for(unsigned int i = 0; i < sz; ++i) {
      std::string name; fin >> name;
      conf.dSolutions.push_back(DeterministicSolver::getId(name));
    }
    fin >> sz;
    for(unsigned int i = 0; i < sz; ++i) {
      std::string name; fin >> name;
      StochasticSolution id = StochasticSolver::getId(name);
      unsigned int paramsSize; fin >> paramsSize;
      for(unsigned int j = 0; j < paramsSize; ++j) {
        std::string key;
        long long value;
        fin >> key >> value;
        conf.sSolutions[id][key] = value;
      }
    }
    return conf;
  }
};

#endif
