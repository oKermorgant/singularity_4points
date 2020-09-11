
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

#include <task.h>
#include <scene.h>
#include <pose_estim.h>

using std::string;
using std::cout;
using std::endl;

// rose params
const int k(7);
const int k2(k-k%2);
const double s(0.1);
const double s0(0.0);

const double dt(M_PI/(50*k));
const double tMax(1.1*(2-(k%2))*M_PI);

vpTranslationVector rose(double t)
{
  return {s0*cos(k2*t) + s*cos(k*t)*cos(t),
        s0*sin(k2*t) + s*cos(k*t)*sin(t),0};
}

int main (int argc, char **argv)
{
  // init and read config file
  log2plot::ConfigManager config(std::string(SRC_PATH) + "/config.yaml");
  config.updateFrom(argc, argv);
  config.forceParameter("control", "rose");

  Scene scene(config);
  scene.initBaseDir(true);

  const auto control = config.read<string>("control");

  auto &point = scene.points;
  const uint n = scene.n_points;

  config.addNameElement(control);

  PoseEstim estimator(scene);
  if(!estimator.use_estim())
  {
    std::cerr << "No estimator given\n";
    return 0;
  }

  log2plot::Logger logger(config.fullName() + "_");

  vpHomogeneousMatrix cMo(scene.cMoFrom(scene.singular_point));
  Task task(point, cMo);

  unsigned int iter = 0;
  vpColVector e(n*2), v(6), cond(2), xy(2*n), te(3);

  vpHomogeneousMatrix oMs;
  oMs.insert(scene.singular_point);
  vpPoseVector pose, pose_cMs, pose_ceMs, pose_inter, pose_cMce;

  logger.save(estimator.Zerr, "Z", "Z_", "<Z-depth error> \\frac{Z-\\tilde{Z}}{Z}[\\%]");

  // real vs estimated camera pose
  if(estimator.refinement == Refinement::None)
    logger.regroupNext(2);
  else
    logger.regroupNext(3);
  logger.save3Dpose(pose_cMs, "cMs", "Real pose", true);
  //logger.showMovingCamera();
  logger.setLineType("[C0,C0,C0,C0]");
  logger.showFixedObject({{0,0,0}}, log2plot::legendFullyDisconnected(1), "rD");

  if(estimator.refinement != Refinement::None)
  {
    logger.save3Dpose(pose_inter, "iMs", estimator.rawMethod(), true);
    logger.setLineType("[C2,C2,C2,C2]");
  }

  //logger.showMovingCamera();
  logger.save3Dpose(pose_ceMs, "ceMs", estimator.fullMethod(), true);
  logger.setLineType("[C1,C1,C1,C1]");


  logger.showFixedObject({{0,0,0}}, log2plot::legendFullyDisconnected(1), "rD");
  logger.setPlotArgs("--legendCol 2 --ae -60 42");
  //logger.showMovingCamera();

  logger.save(te, "te", "[t_x, t_y, t_z]","Position error [m]");
  //logger.setPlotArgs("--yLim -0.02 0.02");
  logger.save(pose_cMce, "cMce",
              "[t_x,t_y,t_z,\\theta u_x,\\theta u_y,\\theta u_z]","Pose error [m,rad]");
  logger.setPlotArgs("--legendCol 2 --ae -60 42");


  logger.save(cond, "cond", "['<distance>', '<cond>{}^{-1}\\mathbf{L}']", "distance to singularity [m]");
  logger.setPlotArgs("--twin 0 --logY2 --twinLabel cond${}^{-1}\\mathbf{L}$ --legendCol 1 --legendLoc out");
  //  logger.saveXY(d_det, "det_d_xy", "", "distance to singularity", "det L^tL");
  logger.save3Dpose(pose, "pose", "Camera pose", true);
  logger.setLineType("[C0,C2--,C3,C4]");

  scene.displayPoints(logger, "C1-s", "rD");

  vpMatrix L;
  bool Zok(true);
  std::vector<double> dist_hist;

  double t(0);
  while (t < tMax && Zok)
  {
    t = (iter++)*dt;
    // update pose
    cMo = scene.cMoFrom(scene.singular_point + rose(t));

    dist_hist.push_back(rose(t).frobeniusNorm());

    pose.buildFrom(cMo);

    auto [cMo_raw, cMo_final] = estimator.computePoseStep(cMo); {}

    if(estimator.refinement == Refinement::None)
    cMo_final = cMo_raw;

    estimator.computeDepthError(cMo, cMo_final);

    pose_ceMs.buildFrom(cMo_final * oMs);
    pose_inter.buildFrom(cMo_raw * oMs);
    pose_cMce.buildFrom(cMo * cMo_final.inverse());

    for(uint ti = 0; ti < 3; ti++)
      te[ti] = pose_cMce[ti];

    pose_cMs.buildFrom(cMo * oMs);


    // compute from cMo to get real L determinant
    task.computeFrom(cMo, cMo);
    // save real values

    // 2D positions
    for (uint i=0;i<n;++i )
    {
      xy[2*i] = point[i].get_x();
      xy[2*i+1] = -point[i].get_y();
    }

    // interaction matrix and determinant
    L = task.L;
    cond[0] = scene.distToSingularity(cMo);
    //d_det[1] = (L.transpose()*L).det();
    cond[1] = 1./L.cond(1e-13);
    //dist_hist.push_back(log10(cond[1]));

    logger.update();
  }

  // automatic threshold from det
  std::vector<double> steps;

  for(uint i = 1; i < dist_hist.size()-1; ++i)
  {
    if(dist_hist[i] < dist_hist[i-1] &&
       dist_hist[i] < dist_hist[i+1])
      steps.push_back(i);
  }
  logger.writeStepsAll(steps, {});

  config.saveConfig();
  logger.plot(true, config.read<std::string>("plot") == "display");
}
