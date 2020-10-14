#ifndef POSEESTIM_H
#define POSEESTIM_H

#include <visp/vpPose.h>
#include <visp/vpGaussRand.h>
#include <scene.h>

enum class Refinement{None, VVS, LM};
enum class Method{None, DEM, P4P, UPnP, EPnP};

class PoseEstim
{
public:

  PoseEstim(Scene &scene, std::string sim_prefix);

  std::string rawMethod() const;
  std::string fullMethod() const;
  std::tuple<vpHomogeneousMatrix, vpHomogeneousMatrix>
  computePoseStep(const vpHomogeneousMatrix &cMo);

  vpHomogeneousMatrix computePose(vpHomogeneousMatrix cMo, bool with_refine = true);
  void computePoseViSP(vpHomogeneousMatrix &cMo, vpPose::vpPoseMethodType visp_method);
  void computePoseOpenGV(vpHomogeneousMatrix &cMo);
  void computePoseEPnP(vpHomogeneousMatrix &cMo);
  void computePoseP4P(vpHomogeneousMatrix &cMo);

  void refine(vpHomogeneousMatrix &cMo);
  void computeDepthError(const vpHomogeneousMatrix &cMo, const vpHomogeneousMatrix &cMo_e);
  vpHomogeneousMatrix poseErrorAt(const vpHomogeneousMatrix &cMo);

  vpHomogeneousMatrix poseError(const vpHomogeneousMatrix &cMo,
                                const vpHomogeneousMatrix cMo_e)
  {
    return  cMo_e * cMo.inverse();
  }

  static std::pair<double, double> errorMetrics(const vpHomogeneousMatrix &cMo,
                                         const vpHomogeneousMatrix cMo_e)
  {
    const auto Merr(cMo_e*cMo.inverse());

    // norm of translation error [m] + abs value of 3D rotation error angle [deg]
    return {Merr.getTranslationVector().frobeniusNorm(),
            std::abs(180/M_PI*acos(0.5*(Merr[0][0] + Merr[1][1] + Merr[2][2]-1)))};
  }

  double reprojectionError(const vpHomogeneousMatrix &cMo,
                                  const vpHomogeneousMatrix cdMo) const;

  bool use_estim() const {return method != Method::None;}
  std::vector<vpPoint>* scene_points;
  vpColVector Zerr;
  std::array<int, 1> n_solutions;
  uint n_points;
  vpGaussRand noise;
  Method method = Method::None;
  Refinement refinement = Refinement::None;

  std::string refineMethod() const
  {
    if(refinement == Refinement::VVS)
      return "VVS";
    if(refinement == Refinement::LM)
      return "LM";
    return "";
  }

};

#endif // POSEESTIM_H
