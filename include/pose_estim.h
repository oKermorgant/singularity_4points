#ifndef POSEESTIM_H
#define POSEESTIM_H

#include <visp/vpPose.h>
#include <opencv2/calib3d.hpp>
#include <visp/vpGaussRand.h>
#include <scene.h>

enum class Refinement{None, VVS, LM};
enum class Method{None, DEM, P4P, UPnP, EPnP};

class PoseEstim
{
public:

  PoseEstim(Scene &scene);

  std::string rawMethod() const;
  std::string fullMethod() const;
  std::tuple<vpHomogeneousMatrix, vpHomogeneousMatrix>
  computePoseStep(const vpHomogeneousMatrix &cMo);

  vpHomogeneousMatrix computePose(vpHomogeneousMatrix cMo, bool with_refine = true);
  void computePoseViSP(vpHomogeneousMatrix &cMo, vpPose::vpPoseMethodType visp_method);
  void computePoseOpenCV(vpHomogeneousMatrix &cMo, int cv_method);
  void computePoseUPnP(vpHomogeneousMatrix &cMo);
  void computePoseP4P(vpHomogeneousMatrix &cMo);

  void refine(vpHomogeneousMatrix &cMo);
  void computeDepthError(const vpHomogeneousMatrix &cMo, const vpHomogeneousMatrix &cMo_e);
  vpHomogeneousMatrix poseErrorAt(const vpHomogeneousMatrix &cMo);

  vpHomogeneousMatrix poseError(const vpHomogeneousMatrix &cMo,
                                const vpHomogeneousMatrix cMo_e)
  {
    return  cMo_e * cMo.inverse();
  }

  bool use_estim() const {return method != "none";}
  std::vector<vpPoint>* scene_points;
  vpColVector Zerr;
  std::vector<uint> reorder;
  uint n_points;
  vpGaussRand noise;
  std::string method = "none";
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
