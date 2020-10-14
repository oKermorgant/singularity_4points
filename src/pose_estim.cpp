#include <pose_estim.h>
#include <numeric>
#include <visp/vpExponentialMap.h>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>
#include <opengv/absolute_pose/methods.hpp>
#include <lambdatwist/pnp/p4p.h>

PoseEstim::PoseEstim(Scene &scene, std::string sim_prefix)
  : scene_points(&scene.points), noise(scene.config.read<double>("noise"),0)
{
  auto ref(scene.config.read<std::string>("ref"));
  if(ref == "lm")
    refinement = Refinement::LM;
  else if(ref == "vvs")
    refinement = Refinement::VVS;

  const auto method_s = scene.config.read<std::string>("estim");
  if(method_s == "dem")
    method = Method::DEM;
  else if(method_s == "epnp")
    method = Method::EPnP;
  else if(method_s == "upnp")
    method = Method::UPnP;
  else if(method_s == "p4p")
    method = Method::P4P;

  // build base directory
  auto dataPath = scene.config.read<std::string>("dataPath");
  dataPath = vpIoTools::path(dataPath);
  if(!vpIoTools::checkDirectory(dataPath))
    dataPath = std::string(SRC_PATH) + "/results";

  // scene level
  dataPath +=  "/scene" + scene.scene_n + "/";

  // estimator level
  if(method != Method::None)
  {
    dataPath += method_s;
    if(refinement != Refinement::None)
      dataPath += "_" + ref;
  }
  else
    dataPath += "no_estim";
  dataPath += "/" + sim_prefix;

  if(method != Method::None)
  {
    if(scene.config.read<double>("noise") != 0)
    {
      auto nz = scene.config.read<std::string>("noise");
      dataPath +=  "_" + nz.substr(2) + "/";
    }
    else
      dataPath += "_noNoise";
  }
  scene.config.setDirName(dataPath);

  std::cout << "Saving to " << dataPath << std::endl;

  if(method == Method::None)
    return;

  n_points = static_cast<uint>(scene_points->size());
  Zerr.resize(n_points);

  std::cout << "Using method: '" << fullMethod()
            << "' on scene " << scene.scene_n
            << " with " << scene.config.read<double>("noise")
            << " noise" << std::endl;
}

std::string PoseEstim::rawMethod() const
{
  switch (method) {
  case Method::DEM:
    return "DeMenthon";
  case Method::P4P:
    return "P4P";
  case Method::UPnP:
    return "UPnP";
  case Method::EPnP:
    return "EPnP";
  default:
    return "unknown";
  }
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

  switch (method) {
  case Method::DEM:
    computePoseViSP(cMo, vpPose::DEMENTHON);
    break;
  case Method::P4P:
    computePoseP4P(cMo);
    break;
  case Method::EPnP:
  case Method::UPnP:
    computePoseOpenGV(cMo);
    break;

  default:
    std::cout << "Method " << rawMethod() << " not implemented\n";
  }

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
    computePoseViSP(cMo, vpPose::VIRTUAL_VS);
  else
    computePoseViSP(cMo, vpPose::LOWE);
}

void PoseEstim::computePoseOpenGV(vpHomogeneousMatrix &cMo)
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

  if(method == Method::EPnP)
  {
    auto M_epnp = opengv::absolute_pose::epnp(adapter);
    for(i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 4; ++j)
        cMo[i][j] = M_epnp(i,j);
    }
    cMo = cMo.inverse();
    return;
  }

  // UPnP case below
  auto solutions(opengv::absolute_pose::upnp(adapter));

  n_solutions[0] = solutions.size();

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
  vpColVector sd(2*n_points);
  uint i(0);

  for(const auto &p: *scene_points)
  {
    bearings.push_back({p.get_x(), p.get_y()});
    points.push_back({p.get_oX(), p.get_oY(), p.get_oZ()});
    sd[2*i] = p.get_x();
    sd[2*i+1] = p.get_y();
    i++;
  }

  // get all combinations of p3p
  cvl::Vector4<uint> indices{0,1,2,3};
  double err(std::numeric_limits<double>::max());
  vpHomogeneousMatrix Mc;

  for(int iter = 0; iter < 4; ++iter)
  {
    const auto p3p_M = cvl::p4p(points, bearings, indices).get3x4();
    // write this solution as homogeneous matrix
    for(int i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 4; ++j)
        Mc[i][j] = p3p_M(i,j);
    }
    // compute reproj error on unused point
    auto &p = scene_points->operator[](indices[3]);
    p.track(Mc);
    const double dx = p.get_x() - sd[2*indices[3]];
    const double dy = p.get_y() - sd[2*indices[3]+1];
    const double this_err(dx*dx + dy*dy);
    if(this_err < err)
    {
      cMo = Mc;
      err = this_err;
    }

    // rotate vector for next iteration
    for(size_t i = 0; i < 3; ++i)
      std::swap(indices[i], indices[i+1]);
  }
}

double PoseEstim::reprojectionError(const vpHomogeneousMatrix &cMo, const vpHomogeneousMatrix cdMo) const
{
  double s(0);

  for(auto &p: *scene_points)
  {
    p.track(cMo);
    double dx(p.get_x()), dy(p.get_y());
    p.track(cdMo);
    dx -= p.get_x();
    dy -= p.get_y();
    s += dx*dx + dy*dy;
  }
  return sqrt(s);

}
