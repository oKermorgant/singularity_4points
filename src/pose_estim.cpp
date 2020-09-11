#include <pose_estim.h>
#include <opencv2/core.hpp>
#include <numeric>
#include <visp/vpExponentialMap.h>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>
#include <opengv/absolute_pose/methods.hpp>
#include <lambdatwist/pnp/p4p.h>

PoseEstim::PoseEstim(Scene &scene)
  : scene_points(&scene.points), noise(scene.config.read<double>("noise"),0)
{
  method = scene.config.read<std::string>("estim");
  auto ref(scene.config.read<std::string>("ref"));
  if(ref == "lm")
    refinement = Refinement::LM;
  else if(ref == "vvs")
    refinement = Refinement::VVS;

  std::string tag = "reorder";
  if(scene.config.has(scene.scene_n + ":" + tag))
    tag = scene.scene_n + ":" + tag;

  reorder = scene.config.read<std::vector<uint>>(tag);

  const std::vector<std::string> known_methods{
    "dem",
    "epnp",
    "p4p",
    "upnp"
  };

  if(std::find(known_methods.begin(), known_methods.end(), method)
     == known_methods.end())
  {
    method = "none";
    return;
  }

  n_points = static_cast<uint>(scene_points->size());
  Zerr.resize(n_points);

  std::string file_info(refineMethod());
  if(scene.config.read<std::string>("control") != "rose")
    file_info = method + file_info;
  scene.config.addNameElement(file_info);

  std::cout << "Using method: '" << fullMethod()
            << "' on scene " << scene.scene_n
            << " with " << scene.config.read<double>("noise")
            << " noise" << std::endl;
}

std::string PoseEstim::rawMethod() const
{
  if(method == "dem")
    return "DeMenthon";
  if(method == "upnp")
    return "UPnP";
  if(method == "epnp")
    return "EPnP";
  if(method == "p4p")
    return "P4P";
  return "unknown";
}

std::string PoseEstim::fullMethod() const
{
  const auto ext = refineMethod();
  return rawMethod() + (ext == "" ? "" : "+") + refineMethod();
}

void PoseEstim::computeDepthError(const vpHomogeneousMatrix &cMo, const vpHomogeneousMatrix &cMo_e)
{
  // compare real vs estimated Z
  for(uint i = 0; i < n_points; ++i)
  {
    auto &p = (*scene_points)[i];
    p.track(cMo);
    const double z(p.get_Z());
    p.track(cMo_e);
    const double ze(p.get_Z());
    Zerr[i] = 100*std::abs(z - ze)/z;
  }
}

vpHomogeneousMatrix PoseEstim::poseErrorAt(const vpHomogeneousMatrix &cMo)
{
  return poseError(cMo, computePose(cMo));
}

std::tuple<vpHomogeneousMatrix, vpHomogeneousMatrix>
PoseEstim::computePoseStep(const vpHomogeneousMatrix &cMo)
{
  const auto cMo_e = computePose(cMo, false);

  auto cMoRef(cMo_e);

  refine(cMoRef);
  return {cMo_e, cMoRef};
}

vpHomogeneousMatrix PoseEstim::computePose(vpHomogeneousMatrix cMo, bool with_refine)
{
  for(auto &p: *scene_points)
  {
    p.track(cMo);
    p.set_x(p.get_x() + noise());
    p.set_y(p.get_y() + noise());
  }

  if(method == "dem")
  {
    computePoseViSP(cMo, vpPose::DEMENTHON);
  }
  else if(method == "epnp")
  {
    computePoseOpenCV(cMo, cv::SOLVEPNP_EPNP);
  }
  else if(method == "p4p")
  {
    computePoseP4P(cMo);
  }
  else if(method == "upnp")
  {
    computePoseUPnP(cMo);
  }
  else
    std::cout << "Method " << method << " not implemented\n";

  if(with_refine)
    refine(cMo);

  return cMo;
}

void PoseEstim::computePoseViSP(vpHomogeneousMatrix &cMo, vpPose::vpPoseMethodType visp_method)
{
  vpPose pose_estim;
  for(auto p: *scene_points)
    pose_estim.addPoint(p);

  pose_estim.computePose(visp_method, cMo);
}

void PoseEstim::refine(vpHomogeneousMatrix &cMo)
{
  if(refinement == Refinement::None)
    return;

  if(refinement == Refinement::VVS)
  {
    const double v_min(1e-10);
    const double lambda (0.9);
    const double iter_max(500);


    const uint np(scene_points->size());
    vpMatrix L(2*np, 6);
    vpColVector sd(2 * np), s(2 * np);
    vpColVector v(6, 1);

    // create sd
    uint k(0);
    std::vector<vpPoint> vvs_points;
    for(auto P: *scene_points)
    {
      sd[2 * k] = P.get_x();
      sd[2 * k + 1] = P.get_y();
      vvs_points.push_back(P);
      k++;
    }

    vpHomogeneousMatrix cMoPrev = cMo;

    uint iter(0);
    while (iter++ < iter_max && v.frobeniusNorm() > v_min)
    {
      k = 0;
      for (auto &P: vvs_points)
      {
        cMo[2][3] = std::max(cMo[2][3], 1e-3);
        P.track(cMo);

        double x = s[2 * k] = P.get_x(); /* point projected from cMo */
        double y = s[2 * k + 1] = P.get_y();
        double Z = P.get_Z();
        L[2 * k][0] = -1 / Z;
        L[2 * k][1] = 0;
        L[2 * k][2] = x / Z;
        L[2 * k][3] = x * y;
        L[2 * k][4] = -(1 + x * x);
        L[2 * k][5] = y;

        L[2 * k + 1][0] = 0;
        L[2 * k + 1][1] = -1 / Z;
        L[2 * k + 1][2] = y / Z;
        L[2 * k + 1][3] = 1 + y * y;
        L[2 * k + 1][4] = -x * y;
        L[2 * k + 1][5] = -x;
        k += 1;
      }
      // compute the VVS control law
      // std::cout << "iter = " << iter << std::endl << L << std::endl << std::endl;
      v = -lambda * L.pseudoInverse() * (s - sd);

      cMoPrev = cMo;
      cMo = vpExponentialMap::direct(v).inverse() * cMo;
    }
  }
  else
  {
    computePoseViSP(cMo, vpPose::LOWE);
  }

}


void PoseEstim::computePoseOpenCV(vpHomogeneousMatrix &cMo, int cv_method)
{
  // build opencv matrices
  std::vector<cv::Point3d> P;
  std::vector<cv::Point2d> xy;

  for(auto i: reorder)
  {
    const auto p = scene_points->operator[](i);
    P.emplace_back(p.get_oX(), p.get_oY(), p.get_oZ());
    xy.emplace_back(p.get_x(), p.get_y());
  }

  cv::Mat K = cv::Mat::eye(3, 3, cv::DataType<double>::type);
  cv::Mat dist_coeffs = cv::Mat::zeros(4,1,cv::DataType<double>::type);
  cv::Mat r, t;

  // outputs oMc
  cv::solvePnP(P, xy, K ,dist_coeffs, r, t, false, cv_method);

  vpThetaUVector tu(r.at<double>(0), r.at<double>(1), r.at<double>(2));
  vpTranslationVector T(t.at<double>(0), t.at<double>(1), t.at<double>(2));

  cMo.buildFrom(T, tu);
}

void PoseEstim::computePoseUPnP(vpHomogeneousMatrix &cMo)
{
  opengv::bearingVectors_t bearings;
  opengv::points_t points;
  vpColVector sd(2*n_points);
  uint i(0);
  for(const auto &p: *scene_points)
  {
    bearings.push_back({p.get_x(), p.get_y(), 1});
    bearings.back() /= bearings.back().norm();
    points.push_back({p.get_oX(), p.get_oY(), p.get_oZ()});

    sd[2*i] = p.get_x();
    sd[2*i+1] = p.get_y();
    i++;
  }

  opengv::absolute_pose::CentralAbsoluteAdapter adapter(bearings, points);

  opengv::transformations_t solutions =
      opengv::absolute_pose::upnp(adapter);


  // keep best
  double best_score(1000);
  vpColVector s(2*n_points);
  cMo.eye();

  for(const auto &M_pnp: solutions)
  {
    // build this homogeneous matrix
    vpHomogeneousMatrix M;
    for(i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 4; ++j)
        M[i][j] = M_pnp(i,j);
    }
    M = M.inverse();
    i = 0;
    for(auto p: *scene_points)
    {
      p.track(M);

      if(p.get_Z() < 0)
        break;

      s[2*i] = p.get_x();
      s[2*i+1] = p.get_y();
      i++;
    }

    if(i == n_points)  // all points could project
    {
      const double score((s-sd).frobeniusNorm());
      if(score < best_score)
      {
        best_score = score;
        cMo = M;
      }
    }
  }
}



void PoseEstim::computePoseP4P(vpHomogeneousMatrix &cMo)
{
  std::vector<cvl::Vector3D> points;
  std::vector<cvl::Vector2D> bearings;
  for(const auto &p: *scene_points)
  {
    bearings.push_back({p.get_x(), p.get_y()});
    points.push_back({p.get_oX(), p.get_oY(), p.get_oZ()});
  }
  const auto p3p_M = cvl::p4p(points, bearings, {0,1,2,3}).get3x4();
  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 4; ++j)
      cMo[i][j] = p3p_M(i,j);
  }
}



