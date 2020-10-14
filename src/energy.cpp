#include <log2plot/logger.h>
#include <task.h>
#include <scene.h>
#include <pose_estim.h>


int main()
{
  // init and read config file
  log2plot::ConfigManager config(std::string(SRC_PATH) + "/config.yaml");
  config.forceParameter("control", "none");
  config.forceParameter("estim", "none");

  config.forceParameter("vs_dist", "0");

  const bool singularScene(true);
  if(singularScene)
    config.forceParameter("scene", "1");
  else
     config.forceParameter("scene", "1o");

  Scene scene(config);
  PoseEstim estimator(scene, "energy");

  const auto [cMos, cdMo] = scene.vsPoses(); {}
  const auto cdRc(cdMo.getRotationMatrix());

  Task task(scene.points, cdMo);
  log2plot::Logger logger(config.fullName());


  constexpr double dist(1);
  constexpr double step(dist/50);
  constexpr int iter(2*(dist/step)+1);

  // energy is visual error at (x,y)
  std::array<double, iter> energy;
  double x;
  // a hack around timed logs - x used for the time...
  logger.setTime(x);
  logger.saveTimed(energy, "e" + config.read<std::string>("vs_dist"), "[]", "");

  // desired position wrt singular point
  vpTranslationVector sTcd(cdMo.inverse());
  sTcd = sTcd - scene.singular_point;
  logger.showFixedObject({{sTcd[0], sTcd[1]}},
                         "[[0,0]]", "gD");

  for(int xi = 0; xi < iter; ++xi)
  {
    x = -dist + xi*step;
    std::cout << "At X = " << x << std::endl;
    for(int yi = 0; yi < iter; ++yi)
    {
      double y = -dist + yi*step;
/*      vpHomogeneousMatrix cMo(-(cdRc.t()*(scene.singular_point + vpTranslationVector{x, y, 0})),
                              cdRc);
      task.computeFrom(cMo);*/
      task.computeFrom(scene.cMoFrom(scene.singular_point + vpTranslationVector{x, y, 0}));
      energy[yi] = log(task.e.frobeniusNorm());
    }
    logger.update();
  }

    std::cout << "Saving to " << config.fullName() +
                 "e" + config.read<std::string>("vs_dist") << ".yaml\n";
}
