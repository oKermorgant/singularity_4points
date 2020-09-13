#include <statistics.h>
#include <algorithm>
#include <numeric>



std::vector<double> medianOfTurns(const std::vector<std::vector<double>> &errors)
{
  uint valid;
  const auto turns(errors.size());

  if(turns == 1)
    return errors[0];

  const auto measures(errors[0].size());
  for(valid = 0; valid < measures; ++valid)
  {
    if(errors[0][valid] == -1)
      break;
  }
  valid++;

  std::vector<double> e(valid);

  std::vector<double> e_point(turns);
  const auto n((turns+1)/2);
  for(uint j = 0; j < valid; ++j)
  {
    for(uint turn = 0; turn < turns; ++turn)
      e_point[turn] = errors[turn][j];

    std::nth_element(e_point.begin(), e_point.begin()+n, e_point.end());
    e[j] = e_point[n];
  }
  return e;
}

double Statistics::Metric::describe(std::vector<double> &e) const
{
  switch (descriptor)
  {
  case Descriptor::PERCENT:
    return std::count_if(e.begin(), e.end(),
                        [&](auto v){return v>thr;})*100./e.size();

  case Descriptor::MAX:
    return *std::max_element(e.begin(), e.end());

  case Descriptor::MEAN:
    return std::accumulate(e.begin(), e.end(), 0) / e.size();

  case Descriptor::MEDIAN:
    const auto n((e.size()+1)/2);
    std::nth_element(e.begin(), e.begin()+n, e.end());
    return e[n];
  }
  return -1;
}

void Statistics::output(std::vector<std::string> metrics,
                        std::string xlabel,
                        std::string symbol,
                        std::string unit,
                        std::string estim_method)
{
  for(auto key: metrics)
    data.push_back({key, {0,0}});

  // build legends
  const std::string full_symbol("|" + symbol + "_e|");

  for(auto &[metric, val]: data)
  {
    if(metric.descriptor == Descriptor::PERCENT)
    {
      std::stringstream ss;
      ss << "$" << full_symbol << " > " << metric.thr;
      ss << "$ " << unit << " [\\%]";
      logger.saveXY(val, metric.file_key, "[<" + estim_method +">]", xlabel, ss.str());
    }
    else
    {
      logger.saveXY(val, metric.file_key, "[<" + estim_method + ">]", xlabel,
                    metric.file_key + " $" + full_symbol + "$ [" + unit  +"]");
      logger.setPlotArgs("--logY");
    }
  }
}

void Statistics::addTerms(const std::vector<double> &errors)
{
  last_terms.push_back(errors);
}


void Statistics::endTurns(double x_ref)
{
  auto stats(medianOfTurns(last_terms));

  for(auto &[metric, val]: data)
  {
    val[0] = x_ref;
    val[1] = metric.describe(stats);
  }

  logger.update();
  last_terms.clear();
}

