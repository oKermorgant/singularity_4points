
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

using std::string;
using std::cout;
using std::endl;

typedef std::tuple<std::string, double> Descriptor;

Descriptor readDes(std::string s)
{
  if(s.find("p_") != s.npos)
  {
    return {"percent", atof(s.substr(2).c_str())};
  }
  else
    return {s, 0.};
}

std::vector<double> medianOfTurns(const vpMatrix &E)
{
  uint valid;
  for(valid = 0; valid < E.getCols(); ++valid)
  {
    if(E[0][valid] == -1)
      break;
  }
  valid++;

  std::vector<double> e(valid);

  std::vector<double> e_point(E.getRows());
  const auto n((E.getRows()+1)/2);
  for(uint j = 0; j < valid; ++j)
  {
    for(uint turn = 0; turn < E.getRows(); ++turn)
      e_point[turn] = E[turn][j];

    std::nth_element(e_point.begin(), e_point.begin()+n, e_point.end());
    e[j] = e_point[n-1];
  }
  return e;
}

double descriptor(std::vector<double> &e, std::string s, double thr = 0)
{
  if(s == "percent")
  {
    double ret(0);
    for(auto v: e)
    {
      if(v > thr)
        ret += 1;
    }
    return ret/e.size()*100;
  }
  else if(s == "mean")
  {
    double ret(0);
    for(auto v: e)
      ret += v;
    return ret/e.size();
  }
  else if(s == "max")
  {
    return *std::max_element(e.begin(), e.end());
  }
  else if(s == "median")
  {
    const auto n((e.size()+1)/2);
    std::nth_element(e.begin(), e.begin()+n, e.end());
    return e[n-1];
  }
  return 0;
}

std::vector<vpTranslationVector> spherePoints(const int denom)
{
  const double dt(M_PI/denom);
  std::vector<vpTranslationVector> T;
  std::cout << "Building sphere points...\n";

  for(int iz = -denom; iz < denom+1; ++iz)
    //for(int iz = 0; iz < denom+1; ++iz)
  {
    const double z = static_cast<double>(iz)/denom;

    const double r = sqrt(1-z*z);
    const double dtr = .5*dt/r;
    const uint curr = T.size();
    for(double t = 0; t < 2*M_PI; t+= dtr)
      T.emplace_back(r*cos(t), r*sin(t), z);
    std::cout << "   z = " << z << " -> using " << T.size()-curr << " points" << std::endl;
  }
  std::cout << "Using " << T.size() << " sphere points\n";
  return T;
}

int main (int argc, char** argv)
{
  // init and read config file
  log2plot::ConfigManager config(std::string(SRC_PATH) + "/config.yaml");
  config.updateFrom(argc, argv);

  Scene scene(config);
  scene.initBaseDir(true);

  std::vector<Descriptor> sphere_des;
  std::vector<vpColVector> sphere_e;
  for(std::string s: {"mean","median","max","p_1","p_10","p_50","p_100"})
  {
    sphere_des.push_back(readDes(s));
    sphere_e.push_back(vpColVector(2));
  }
  const uint n_des(sphere_des.size());

  // file name comes from singular point and distance to it
  auto &point = scene.points;
  const vpTranslationVector oTs(scene.singular_point);

  config.addNameElement("sphere");
  PoseEstim estimator(scene);
  if(!estimator.use_estim())
  {
    std::cout << "Estimation method not valid\n";
    return 0;
  }

 // auto method = config.read<std::string>("estim");
 // if(config.read<bool>("vvs"))
 //   method += "+VVS";

  log2plot::Logger logger(config.fullName() + "_");
  vpPoseVector pose;
  const double rmax(config.read<double>("sphere_r"));
  double dr(std::max(0.02, rmax/20));

  std::cout << "Sphere of " << rmax << " m around point " << oTs.t() << std::endl;

  for(uint i = 0; i< n_des; ++i)
  {
    auto &[des, thr] = sphere_des[i]; {}
    if(thr)
    {
      std::stringstream ss, ss_l;
      ss << "$||t_e||>";
      ss << thr;
      ss << "$ mm [\\%]";
      ss_l << "p" << thr;
      logger.saveXY(sphere_e[i], ss_l.str()+"_err", "[<" + estimator.fullMethod() +">]", "sphere radius [m]", ss.str());
      thr /= 1000;
    }
    else
    {
      logger.saveXY(sphere_e[i], des + "_err", "[<" + estimator.fullMethod() + ">]", "sphere radius [m]", des + " $||t_e||$ [m]");
      logger.setPlotArgs("--logY");
    }
  }

  // logger.showMovingCamera();
  //available.displayPoints(logger, "C1-s", "rD");

  vpHomogeneousMatrix cMo;
  double err;

  const auto T(spherePoints(config.read<int>("sphere_n")));

  pose.buildFrom(cMo);

  log2plot::Logger sphere3D(config.fullName() + "_");
  vpColVector dummy(6);
  sphere3D.save3Dpose(dummy, "3D", estimator.fullMethod());
  sphere3D.setLineType("[b]");

  std::vector<std::vector<double>> M(2, {0,0,0});

  const uint turns(config.read<double>("noise") == 0. ? 1 : config.read<uint>("sphere_turns"));
  vpMatrix Errors(turns, T.size());
  double r = 0;
  while(r < rmax+dr)
  {
    std::cout << "Computing for radius = " << r << std::endl;

    const bool do_3D(r < 0.5 && r+dr >= 0.5);

    uint tc(0);
    Errors = -1;

    for(uint turn = 0; turn < turns; ++turn)
    {
      uint err_idx = 0;
      for(const auto &sTc: T)
      {
        if(scene.distToClosestPoint(oTs + sTc*r) < 2)
          continue;

        tc++;
        cMo = scene.cMoFrom(oTs + sTc*r);

        vpHomogeneousMatrix cMo_e(estimator.computePose(cMo));

        if(do_3D)
          sphere3D.showFixedObject<vpTranslationVector>(
          {sTc*r, cMo_e.inverse().getTranslationVector()-oTs},
                "[[0,1]]", "b.-");

        err = estimator.poseError(cMo, cMo_e).getTranslationVector().frobeniusNorm();
        Errors[turn][err_idx++] = err;
      }
    }

    auto errors = medianOfTurns(Errors);

    // save into given vectors
    for(uint i = 0; i < n_des; ++i)
    {
      const auto [des, thr] = sphere_des[i]; {}
      sphere_e[i][0] = r;
      sphere_e[i][1] = descriptor(errors, des, thr);
    }

    logger.update();
    if(do_3D)
      sphere3D.update();

    if(r < 0.2)
      r += 0.02;
    else
      r += dr;
  }

  const auto plot = config.read<std::string>("plot");
  if(plot != "none")
  {
    sphere3D.plot(true, plot == "display");
    logger.plot(true, plot == "display");
  }
}

