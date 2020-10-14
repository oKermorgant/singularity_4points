#include <log2plot/logger.h>
#include <task.h>
#include <scene.h>
#include <pose_estim.h>
#include <visp/vpRobotCamera.h>

using std::string;
using std::cout;
using std::endl;

struct PoseLabel
{
  std::string camlab;
  std::string index;
  PoseLabel(uint i)
  {
    std::stringstream ss;
    ss << i+1;
    camlab = ss.str();
    ss.str("");
    ss << i;
    index = ss.str();
  }
  std::string color() const
  {
    return "C" + index;
  }
  std::string pref(std::string start) const
  {
    return lab(start, "");
  }
  std::string suf(std::string end) const
  {
    return lab("c", end);
  }
  std::string lab(std::string start, std::string end) const
  {
    return start + camlab + end;
  }
};

void showCam(log2plot::Logger &logger, vpHomogeneousMatrix cMo, std::string color, double scale = 0.07)
{
  cMo = cMo.inverse();
  const double x(1.5*scale), y(1*scale), z(4*scale);

  // init points
  std::vector< std::vector<double> > Mv(5, std::vector<double>(3,0));
  Mv[1][0] = x;    Mv[1][1] = -y;    Mv[1][2] = z;
  Mv[2][0] = -x;   Mv[2][1] = -y;    Mv[2][2] = z;
  Mv[3][0] = -x;   Mv[3][1] = y;     Mv[3][2] = z;
  Mv[4][0] = x;    Mv[4][1] = y;     Mv[4][2] = z;

  for(uint i = 0; i < 5; ++i)
  {
    vpPoint P(Mv[i][0], Mv[i][1], Mv[i][2]);
    P.track(cMo);
    Mv[i][0] = P.get_X();
    Mv[i][1] = P.get_Y();
    Mv[i][2] = P.get_Z();
  }
  logger.showFixedObject(Mv, "[[0,1],[0,2],[0,3],[0,4],[1,2],[2,3],[3,4],[4,1]]", color);
}

vpHomogeneousMatrix finalPose(vpRobotCamera &robot, Task &task, vpHomogeneousMatrix cMo,
                              double lambda = 1., double errMin = 1., uint maxIter = 1000)
{
  double err(2*errMin);
  uint iter(0);
  robot.setPosition(cMo);
  while( iter++ < maxIter && err > errMin)
  {
    robot.getPosition ( cMo );
    task.computeFrom(cMo, cMo);
    robot.setVelocity(vpRobot::CAMERA_FRAME, -lambda * task.L.pseudoInverse() * task.e);
    err = task.e.frobeniusNorm();
  }
  return cMo;
}

int main (int argc, char **argv)
{
  // init and read config file
  log2plot::ConfigManager config(std::string(SRC_PATH) + "/config.yaml");
  config.updateFrom(argc, argv);
  config.forceParameter("control", "pinv");

  Scene scene(config);
  PoseEstim estimator(scene, "visualservo");

  const auto errMin(config.read<double>("errMin"));
  const auto dt(config.read<double>( "dt"));
  const auto max_vel(config.read<double>( "max_vel"));
  const auto maxT(config.read<double>( "maxT"));

  auto [cMo, cdMo] = scene.vsPoses(); {}
  const auto m(cMo.size());

  auto &point = scene.points;

  const auto control = config.read<string>("control");
  const auto lambda = config.read<double>("lambda");

  const uint n = scene.n_points;



  log2plot::Logger logger(config.fullName() + "_");

  // robot
  std::vector<vpRobotCamera> robot(m);
  for(uint i = 0; i < m; ++i)
  {
    robot[i].setSamplingTime ( dt );
    robot[i].setPosition ( cMo[i]);
    robot[i].setMaxRotationVelocity(max_vel*robot[i].getMaxRotationVelocity());
    robot[i].setMaxTranslationVelocity(max_vel*robot[i].getMaxTranslationVelocity());
  }

  Task task(point, cdMo);

  std::vector<bool> local_min(m);
  for(uint i = 0; i < m; ++i)
  {
    const auto Merr = finalPose(robot[0], task, cMo[i], lambda, errMin, 2000)* cdMo.inverse();
    local_min[i] = Merr.getTranslationVector().frobeniusNorm() > 1e-1 ||
        std::abs(Merr.getThetaUVector().getTheta()) > 1e-1;
  }

  double err = 2*errMin;
  vpColVector e(n*2), v(6), cond(2), xy(2*m*n);

  vpHomogeneousMatrix oMs(scene.singular_point, vpRotationMatrix());

  std::vector<vpPoseVector> pose(m), cMs(m);

  // absolute poses including observed object
  logger.regroupNext(m);
  for(uint i = 0; i < m; ++i)
  {
    // camera poses from multi start
    const PoseLabel lab(i);
    logger.save3Dpose(pose[i], lab.pref("pose"), lab.lab("'Start ", "'"), true);
    pose[i].buildFrom(cMo[i]);
    logger.setLineType("[" + lab.color() + ",C2--,C3,C4]");
    showCam(logger, cMo[i], lab.color());

    if(local_min[i])
      showCam(logger, finalPose(robot[0], task, cMo[i], lambda, errMin, 2000),
          lab.color() + "--");

    if(i == m-1)
    {
      scene.displayPoints(logger, "C1-s", "rD");
      logger.setPlotArgs("--legendCol 2 --ae -134 49");
      showCam(logger, cdMo, "k-");
    }
  }

  // relative poses around singularity
  logger.regroupNext(m);
  std::vector<double> dist_singularity(m);

  std::stringstream xy_labels, xy_colors, ss_starts;
  xy_labels << "[";
  xy_colors << "[";
  ss_starts << "[";
  for(uint i = 0; i < m; ++i)
  {
    const PoseLabel lab(i);
    cMs[i].buildFrom(cMo[i] * oMs);
    dist_singularity[i] = vpHomogeneousMatrix(cMs[i]).inverse().getTranslationVector().frobeniusNorm();
    logger.save3Dpose(cMs[i], lab.lab("c", "Ms"), lab.lab("'Start ", "'"), true);
    logger.setLineType("[" + lab.color() + ",C2--,C3,C4]");
    showCam(logger, vpHomogeneousMatrix(cMs[i]), lab.color());

    vpPoseVector pose(cdMo * oMs);
    //logger.showMovingCamera({des[0], des[1], des[2], des[3], des[4], des[5]});

    if(local_min[i])
      showCam(logger,
              finalPose(robot[0], task, cMo[i], lambda, errMin, 2000) * oMs,
          lab.color() + "--");

    xy_labels << lab.lab("'<Start ", ">'");
    xy_colors << lab.color();

    ss_starts << lab.lab("'<Start ", ">'");

    if(i == m-1)
    {
      logger.showFixedObject({{0,0,0}}, log2plot::legendFullyDisconnected(1), "rD");
      logger.setPlotArgs("--legendCol 2 --ae 7 46");
      showCam(logger, cdMo*oMs, "k-");
    }
    else
    {
      ss_starts <<", ";

    }
    xy_colors <<", ";
    xy_labels << ", ";

    for(uint k = 1; k < 4; ++k)
    {
      std::string comma = (k == 3 && i == m-1)?"":", ";
      xy_labels << "'' " << comma;
      xy_colors << lab.color() <<  comma;
    }
  }
  xy_labels << "]";
  xy_colors << "]";
  ss_starts << "]";
  robot[0].setPosition(cMo[0]);

  // XY trajectories
  std::vector<std::vector<double>> xyd(n, {0,0});
  for (uint k=0;k<n;++k)
  {
    xyd[k][0] = task.pd[k].get_x();
    xyd[k][1] = -task.pd[k].get_y();
  }
  logger.saveXY(xy, "xy", xy_labels.str(), "", "");
  logger.setLineType(xy_colors.str());
  logger.showFixedObject(xyd, log2plot::legendFullyDisconnected(n), "kd");

  // get closest start from singularity
  const auto closest = std::distance(dist_singularity.begin(),
                                     std::min_element(
                                       dist_singularity.begin(),
                                       dist_singularity.end()));

  xy_labels.str("");
  xy_colors.str("");
  xy_labels << "[<Start 1>, <Start " << closest+1 << ">]";
  xy_colors << "[C0, C" << closest << "]";

  std::stringstream iter_zoom;
  iter_zoom << 10*static_cast<uint>(0.05/dt)+1;

  logger.save(cond, "det_d", xy_labels.str(), "<cond>{}^{-1}\\mathbf{L}");
  logger.setLineType(xy_colors.str());
  logger.setPlotArgs("--logY --index " + iter_zoom.str());

  vpMatrix L;
  bool Zok(true);
  std::vector<double> det_hist;
  vpColVector vi, v1, Ze(4);
  std::cout << "Plotting velocities from start #" << closest +1<< std::endl;
  logger.save(vi, "vel", "[v_x,v_y,v_z,\\omega_x,\\omega_y,\\omega_z]", "Velocity [m/s,rad/s]");
  logger.setLineType("[C0,C1,C2,C3,C4,C5]");
  logger.setPlotArgs("--legendCol 2 --index " + iter_zoom.str());
  logger.save(v1, "vel1", "[v_x,v_y,v_z,\\omega_x,\\omega_y,\\omega_z]", "Velocity [m/s,rad/s]");
  logger.setLineType("[C0,C1,C2,C3,C4,C5]");
  logger.setPlotArgs("--legendCol 2 --index " + iter_zoom.str());
  if(estimator.use_estim())
  {
    logger.save(Ze, "Z", "Z_", "<Z-depth error> \\frac{Z-\\tilde{Z}}{Z}[\\%]");
    logger.setPlotArgs("--index " + iter_zoom.str());
  }

  vpColVector e_all(m);
  logger.save(e_all, "err", ss_starts.str(), "Visual error norm");
  logger.setPlotArgs("--logY");
  vpColVector te(m);
  logger.save(te, "t", ss_starts.str(), "Position error norm [m]");
  logger.setPlotArgs("--legendCol 2");

  double t(0);
  while ( t < maxT && err > errMin && Zok)
  {
    err = 0.;
    t += dt;
    // update poses
    for(uint i = 0; i < m; ++i)
    {
      robot[i].getPosition ( cMo[i] );
      pose[i].buildFrom(cMo[i]);
      cMs[i].buildFrom(cMo[i] * oMs);
      te[i] = (cMo[i]*cdMo.inverse()).getTranslationVector().frobeniusNorm();

      // compute from cMo to get real L determinant
      vpHomogeneousMatrix cMo_e(cMo[i]);
      if(estimator.use_estim())
        cMo_e = estimator.computePose(cMo[i]);
      task.computeFrom(cMo[i], cMo_e);
      // 2D positions
      for (uint k=0;k<n;++k)
      {
        xy[2*(k+i*n)] = point[k].get_x();
        xy[2*(k+i*n)+1] = -point[k].get_y();
      }
      L = task.L;
      e = task.e;

      err = std::max(e.frobeniusNorm(), err);
      e_all[i] = e.frobeniusNorm();

      // gain
      v = -lambda * L.pseudoInverse() * e;
      robot[i].setVelocity(vpRobot::CAMERA_FRAME, v);
      if(i == closest)
      {
        cond[1] = 1./L.cond();
        vi = v;
        if(estimator.use_estim())
        {
          estimator.computeDepthError(cMo[i], cMo_e);
          Ze = estimator.Zerr;
        }
      }
      else if(i == 3)
      {
        cond[0] = 1./L.cond();
        v1 = v;
      }
    }
    logger.update();

  }

  for(uint i = 0; i < m; ++i)
  {
    task.computeFrom(cMo[i]);
    std::cout << "Final error / start " << i+1 << " -> " << task.e.frobeniusNorm() << std::endl;
  }

  config.saveConfig();
  logger.plot(true, config.read<std::string>( "plot") == "display");
}
