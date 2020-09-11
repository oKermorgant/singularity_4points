#ifndef SINGCONFIG_H
#define SINGCONFIG_H

#include <visp/vpPoint.h>
#include <visp/vpIoTools.h>
#include <visp/vpTranslationVector.h>
#include <log2plot/logger.h>
#include <log2plot/config_manager.h>
#include <log2plot/dir_tools.h>

namespace
{
double dist(const vpPoint &p, const vpTranslationVector &t)
{
  vpTranslationVector t2(p.get_oX(), p.get_oY(), p.get_oZ());
  return (t-t2).frobeniusNorm();
}
}

class Scene
{
public:
  Scene(log2plot::ConfigManager &_config);

  void initBaseDir(bool use_estim = true)
  {
    auto dataPath = config.read<std::string>("dataPath");
    dataPath = vpIoTools::path(dataPath);
    if(!vpIoTools::checkDirectory(dataPath))
      dataPath = std::string(SRC_PATH) + "/results";

    dataPath +=  "/config" + scene_n + "/";
    if(use_estim)
      dataPath += config.read<std::string>("estim") + "/";

    if(config.read<double>("noise") != 0)
    {
      auto nz = config.read<std::string>("noise");
      dataPath += "noise_" + nz.substr(2) + "/";
    }
    config.setDirName(dataPath);
  }

  void buildBaseName(const vpTranslationVector &oTc, double d = 0);

  double distToSingularity(const vpTranslationVector &oTc) const
  {
    return (singular_point-oTc).frobeniusNorm();
  }
  double distToSingularity(const vpHomogeneousMatrix &cMo) const
  {
    return distToSingularity(cMo.inverse().getTranslationVector());
  }

  double distToClosestPoint(const vpTranslationVector &oTc) const
  {
    auto closest = *std::min_element(points.begin(),
                                     points.end(),
                                     [&](const vpPoint&p1,
                                     const vpPoint &p2)
    {
      return dist(p1, oTc) < dist(p2, oTc);});
    return dist(closest, oTc);
  }

  // full pose from single position
  vpHomogeneousMatrix cMoFrom(const vpTranslationVector &oTc,
                                     const vpHomogeneousMatrix &cMo_ref,
                                     bool use_ref = true);

  vpHomogeneousMatrix cMoFrom(const vpTranslationVector &oTc)
  {
    return cMoFrom(oTc, vpHomogeneousMatrix(), false);
  }

  void displayPoints(log2plot::Logger &logger,
                     std::string point_color, std::string singularity_color);

  // VS stuff
  std::pair<std::vector<vpHomogeneousMatrix>, vpHomogeneousMatrix>
  vsPoses();

  log2plot::ConfigManager &config;
  std::vector<vpPoint> points;
  vpTranslationVector singular_point;
  vpTranslationVector centroid;
  std::string scene_n;
  uint n_points;

  bool p0 = false;

  template <class T>
  T read(std::string key) const
  {
    return config.read<T>({scene_n, key});
  }
};



#endif
