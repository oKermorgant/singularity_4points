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

enum class Phase{
  BRING_TO_FRONT,
  CENTER_IMAGE_POINTS,
  ALIGN_REF
};


using AutoPose = std::pair<vpHomogeneousMatrix, vpColVector>;


std::vector<AutoPose> bestOrientationDetail(vpHomogeneousMatrix &oMc,
                                            std::vector<vpPoint>&P,
                                            const vpHomogeneousMatrix &cMo_ref,
                                            bool use_ref)
{
  const uint n = static_cast<uint>(P.size());

  vpMatrix L(n, 3), H(n,n);
  vpColVector s(n), sd(n);
  vpThetaUVector w;
  const double l = 0.1;
  double Zmin(0);
  vpRotationMatrix R;

  std::vector<AutoPose> res;
  vpColVector xy(2*n);

  std::vector<Phase> phases;
    phases.push_back(Phase::BRING_TO_FRONT);

  phases.push_back(Phase::CENTER_IMAGE_POINTS);
  if(use_ref)
    phases.push_back(Phase::ALIGN_REF);

  for(auto phase: phases)
  {
    uint iter(0);

    if(phase == Phase::CENTER_IMAGE_POINTS)
    {
      L.resize(2*n, 3);
      H.eye(2*n);
      s.resize(2*n);
      sd.resize(2*n);
    }
    else if(phase == Phase::ALIGN_REF)
    {
      L.resize(n, 3);
      H.eye(n);
      s.resize(n);
      sd.resize(n);
      for(uint i = 0; i < n; ++i)
      {
        P[i].track(cMo_ref);
        sd[i] = atan2(P[i].get_y(), P[i].get_x());
        L[i][2] = -1;
      }
    }
    std::pair<double, uint> maxCoord{-1, 0};
    while(iter++ < 1000)
    {
      // get interaction matrix and weights
      for(uint i = 0; i < n; ++i)
      {
        P[i].track(oMc.inverse());

        const double X = P[i].get_X();
        const double Y = P[i].get_Y();
        const double Z = P[i].get_Z();
        const double x = X/Z;
        const double y = Y/Z;

        if(Z < 0 || std::abs(x) > 1 || std::abs(y) > 1)
        {
          xy[2*i] = xy[2*i+1] = std::nan("");
        }
        else
        {
          xy[2*i] = x;
          xy[2*i+1] = y;
        }

        switch(phase)
        {
        case Phase::BRING_TO_FRONT:
          // bring Z in front of the camera
          L[i][0] = -Y;
          L[i][1] = X;
          s[i] = Z - 1;

          if(Z < 0)
            H[i][i] = 1;
          else if(Z > 0.5)
            H[i][i] = .1;
          else
            H[i][i] = 1 - Z/5;
          break;
        case Phase::CENTER_IMAGE_POINTS:
          // bring all points to the center
          s[2*i] = x;
          s[2*i+1] = y;
            L[2*i][0] = x*y;
            L[2*i][1] = -(1+x*x);
            L[2*i+1][0] = 1+y*y;
            L[2*i+1][1] = -x*y;
          L[2*i][2] = y;
          L[2*i+1][2] = -x;
          maxCoord = std::max(maxCoord, {std::abs(x), 2*i});
          maxCoord = std::max(maxCoord, {std::abs(y), 2*i+1});
          break;
        default:
          double rho = sqrt(x*x + y*y);
          double theta = atan2(y, x);
          s[i] = theta;
            L[i][0] = cos(theta)/rho;
            L[i][1] = sin(theta)/rho;
        }
      }

      Zmin = std::min_element(P.begin(), P.end(),
                              [](const vpPoint &p1, const vpPoint &p2)
      {
        return p1.get_Z() < p2.get_Z();
      })->get_Z();

      if(phase == Phase::BRING_TO_FRONT && Zmin > 0.3)
        break;

      if(phase == Phase::CENTER_IMAGE_POINTS)
      {
        //H = 0.;
        //H[maxCoord.second][maxCoord.second] = 1;
      }

      // angular velocity
      w = -l * (H*L).pseudoInverse() * H*(s-sd);

      if(std::abs(w.getTheta()) < 1e-3)
        break;

      R.buildFrom(w);

      oMc.insert(oMc.getRotationMatrix()*R);
      res.push_back({oMc.inverse(), xy});
    }
  }
  oMc = oMc.inverse();
  return res;
}

int main (int argc, char **argv)
{
  // init and read config file
  log2plot::ConfigManager config(std::string(SRC_PATH) + "/config.yaml");
  config.updateFrom(argc, argv);
  config.forceParameter("control", "none");
  config.forceParameter("estim", "none");

  Scene scene(config);
  PoseEstim estimator(scene, "align_demo");

  vpHomogeneousMatrix cdMo, cMo;

  cdMo.insert(scene.singular_point);
  cMo.insert(scene.singular_point + vpTranslationVector{0.5, -0.5, -0.5});

  auto &point = scene.points;

  const auto resDes = bestOrientationDetail(cdMo, point, cdMo, false);
  const auto resCur = bestOrientationDetail(cMo, point, resDes.back().first, true);

  log2plot::Logger logger(config.fullName());
  vpColVector XY(4*4);

  logger.saveXY(XY, "xy", "[$s$, '', '', '', $s^*$, '', '', '']", "", "");
  logger.setLineType("[C0, C0, C0, C0, C1, C1, C1, C1]");
  logger.showFixedObject({{-.4,.4},{.4,.4},{.4,-.4},{-.4,-.4}}, "[[0,1],[1,2],[2,3],[3,0]]", "k");
  logger.setPlotArgs("--xLim -1 0.5 --yLim -0.5 1");


  logger.regroupNext(2);
  vpPoseVector pose, poseD;
  logger.save3Dpose(pose, "cMo", "Initial pose", true);
  logger.showMovingCamera();
  logger.setLineType("[C0, C0, C0, C0]");
  logger.save3Dpose(poseD, "cdMo", "Desired pose", true);
  logger.showMovingCamera();
  logger.setLineType("[C1, C1, C1, C1]");
  scene.displayPoints(logger, "C2-s", "rD");

  log2plot::setNaN(XY);

  pose.buildFrom(resCur.front().first);

  // align desired pose
  for(const auto &[M, xy]: resDes)
  {
    poseD.buildFrom(M);
    for(uint i = 0; i < 4; ++i)
    {
      XY[8+2*i] = xy[2*i];
      XY[8+2*i+1] = xy[2*i+1];
    }
    logger.update();
  }

  // align current pose
  for(const auto &[M, xy]: resCur)
  {
    pose.buildFrom(M);
    for(uint i = 0; i < 4; ++i)
    {
      XY[2*i] = xy[2*i];
      XY[2*i+1] = xy[2*i+1];
    }
    logger.update();
  }

  logger.plot();




}
