#ifndef SPHERE_TOOLS_H
#define SPHERE_TOOLS_H

#include <visp/vpTranslationVector.h>

std::vector<vpTranslationVector> goldenSphericalSpiral(uint n)
{
  static const double angle_incr(4*M_PI/(1 + sqrt(5)));

  std::vector<vpTranslationVector> points;
  points.reserve(n);

  const double n_inv(1./n);

  for(uint i = 0; i < n; ++i)
  {
    const double z(i*2*n_inv -1 + n_inv);
    double rho(sqrt(1-z*z));
    const double angle(i*angle_incr);
    points.emplace_back(rho*cos(angle), rho*sin(angle), z);
  }
  return points;
}


std::vector<vpTranslationVector> spherePoints(double r, double dmin)
{
  if(r == 0)
    return {{0,0,0}};

  // generate n points on the sphere
  // distance to 3 closest neighboors is less than dmin
  const uint n = std::max<uint>(100, pow(4.1*r/dmin, 2));
  std::vector<vpTranslationVector> points;
  points.reserve(n);

  // use golden spiral
  static const double angle_incr(4*M_PI/(1 + sqrt(5)));
  const double n_inv(1./n);

  for(uint i = 0; i < n; ++i)
  {
    const double z(i*2*n_inv + n_inv - 1);
    const double rho(sqrt(1-z*z));
    const double angle(i*angle_incr);
    points.emplace_back(r*rho*cos(angle), r*rho*sin(angle), r*z);
  }
  return points;
}


#endif // SPHERE_TOOLS_H
