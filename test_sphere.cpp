#include "src/sphere_tools.h"


int main()
{
  auto t = spherePoints(0.7, 0.1);
  std::cout << "n points: " << t.size() << std::endl;

  double m = 0.1;
  for(size_t i = 0; i < t.size()-1; ++i)
  {
    for(size_t j = i+1; j < t.size(); ++j)
    {
      m = std::min(m, (t[i]-t[j]).frobeniusNorm());
    }
  }
  std::cout << "min dist: " << m << std::endl;
}
