#ifndef METRICS_H
#define METRICS_H

#include <log2plot/logger.h>
#include <map>

class Statistics
{

  enum class Descriptor{MEAN, MEDIAN,MAX,PERCENT};

  struct Metric
  {
    Descriptor descriptor;
    double thr = 0;
    std::string file_key;
    Metric(const std::string &key)
    {
      if(key.find("p_") != key.npos)
      {
        descriptor = Descriptor::PERCENT;
        thr = atof(key.substr(2).c_str());
        file_key="p" + key.substr(2);
        return;
      }
      else if(key == "mean")
        descriptor = Descriptor::MEAN;
      else if(key == "max")
        descriptor = Descriptor::MAX;
      else if(key == "median")
        descriptor = Descriptor::MEDIAN;
      else
        std::cout << "Unknown metric descriptor: " << key << std::endl;
      file_key=key;
    }
    double describe(std::vector<double> &e) const;
  };

  std::vector<std::pair<const Metric, std::array<double, 2>>> data;
  std::vector<std::vector<double>> last_terms;

  log2plot::Logger logger;


public:
  Statistics(std::string save_path) : logger(save_path) {}

  void output(std::vector<std::string> metrics, std::string xlabel, std::string symbol, std::string unit, std::string estim_method);

  void addTerms(const std::vector<double> &errors);
  void endTurns(double x_ref);

  void plot(bool display)
  {
    logger.plot(true, display);
  }




};

#endif // METRICS_H
