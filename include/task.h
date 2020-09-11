#ifndef TASK_H
#define TASK_H

#include <visp/vpServo.h>
#include <visp/vpFeaturePoint.h>
#include <visp/vpFeatureBuilder.h>

struct Task
{
  Task(std::vector<vpPoint> &point, const vpHomogeneousMatrix &cdMo)
    : P(&point)
  {
    task.setServo(vpServo::EYEINHAND_CAMERA);
    task.setInteractionMatrixType(vpServo::CURRENT);

    p.resize(P->size());
    pd.resize(P->size());

    for(uint i = 0; i < P->size(); ++i)
    {
      point[i].track (cdMo);
      vpFeatureBuilder::create ( pd[i],point[i] );
      task.addFeature(p[i], pd[i]);
    }
  }

  void computeFrom(const vpHomogeneousMatrix &M)
  {
    computeFrom(M, vpHomogeneousMatrix(), false);
  }

  void computeFrom(const vpHomogeneousMatrix &M, const vpHomogeneousMatrix &Me, bool use_estim = true)
  {
    for (uint i=0;i<P->size();++i )
    {
      auto &pt = P->at(i);
      if(use_estim)
      {
        pt.track(Me);
        const double z = pt.get_Z();
        pt.track(M);
        vpFeatureBuilder::create ( p[i],pt);
        p[i].set_Z(z);
      }
      else
      {
        pt.track(M);
        vpFeatureBuilder::create ( p[i],pt);
      }
    }
    task.computeControlLaw();
    e = task.computeError();
    L = task.computeInteractionMatrix();
  }


  std::vector<vpPoint> * P;
  std::vector<vpFeaturePoint> p, pd;
  vpServo task;
  vpMatrix L;
  vpColVector e;
};

#endif // TASK_H
