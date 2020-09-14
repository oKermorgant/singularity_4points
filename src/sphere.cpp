
#include <stdlib.h>

#include <visp/vpImage.h>
#include <visp/vpImageIo.h>
//#include <visp/vpDisplayX.h>
#include <visp/vpDisplayOpenCV.h>
#include <visp/vpCameraParameters.h>

#include <visp/vpMath.h>
#include <visp/vpSubMatrix.h>
#include <visp/vpHomogeneousMatrix.h>
#include <visp/vpParseArgv.h>
#include <visp/vpIoTools.h>
#include <visp/vpWireFrameSimulator.h>
#include <visp/vpRobotCamera.h>

#include <visp/vpPose.h>

#include <log2plot/logger.h>

#include <scene.h>
#include <pose_estim.h>
#include <string>
#include <task.h>
#include "sphere_tools.h"

#include <statistics.h>

using std::string;
using std::cout;
using std::endl;


int main (int argc, char** argv)
{
  // init and read config file
  log2plot::ConfigManager config(std::string(SRC_PATH) + "/config.yaml");
  config.updateFrom(argc, argv);

  Scene scene(config);
  PoseEstim estimator(scene, "sphere");

  if(!estimator.use_estim())
  {
    std::cout << "Estimation method not valid\n";
    return 0;
  }

  Statistics transStats(config.fullName() + "t_"),
      rotStats(config.fullName()+"r_");

  transStats.output({"mean","median","p_1","p_10","p_50"}, "Radius [m]", "t", "mm", estimator.fullMethod());
  rotStats.output({"mean","median","p_1","p_10","p_90"}, "Radius [m]", "\\theta", "deg", estimator.fullMethod());

  vpPoseVector pose;
  const double rmax(config.read<double>("sphere_r"));
  const double dmin(config.read<double>("sphere_d"));

  const vpTranslationVector oTs(scene.singular_point);

  std::cout << "Sphere of " << rmax << " m around point " << oTs.t() << std::endl;

  vpHomogeneousMatrix cMo;

  pose.buildFrom(cMo);

  log2plot::Logger sphere3D(config.fullName() + "sphere");
  vpColVector dummy(6);
  sphere3D.save3Dpose(dummy, "3D", estimator.fullMethod());
  sphere3D.setLineType("[b]");

  std::vector<std::vector<double>> M(2, {0,0,0});

  const uint turns(config.read<double>("noise") == 0. ? 1 : config.read<uint>("sphere_turns"));

  double r = 0.001;
  const double dr(0.01);
  while(r < rmax+dr)
  {
    const bool do_3D(r < 0.1 && r+dr >= 0.1);

    const auto T = spherePoints(r, dmin);

    std::cout << "Using " << T.size() << " points for radius " << r << std::endl;
    std::vector<double> errT(T.size(),-1);
    std::vector<double> errR(T.size(),-1);

    for(uint turn = 0; turn < turns; ++turn)
    {
      uint err_idx = 0;
      for(const auto &sTc: T)
      {
        if(scene.distToClosestPoint(oTs + sTc) < 2)
          continue;

        cMo = scene.cMoFrom(oTs + sTc);

        vpHomogeneousMatrix cMo_e(estimator.computePose(cMo));

        if(do_3D)
          sphere3D.showFixedObject<vpTranslationVector>(
          {sTc, cMo_e.inverse().getTranslationVector()-oTs},
                "[[0,1]]", "b.-");

        const auto [te,re] = estimator.errorMetrics(cMo, cMo_e); {}
        errT[err_idx] = te*1000;
        errR[err_idx] = re;
        err_idx++;
      }
      transStats.addTerms(errT);
      rotStats.addTerms(errR);
    }

    transStats.endTurns(r);
    rotStats.endTurns(r);

    if(do_3D)
      sphere3D.update();

    if(r == 0.001)
      r = 0;
    r += dr;
  }

  const auto plot = config.read<std::string>("plot");
  if(plot != "none")
  {
    sphere3D.plot(true, plot == "display" || plot == "create");
    transStats.plot(plot=="display");
    rotStats.plot(plot=="display");
  }
}

