#include <scene.h>

using namespace std;

bool comparePoints(const vpPoint &p1, const vpPoint &p2)
{
  return p1.get_X() + 2*p1.get_Y() + 3*p1.get_Z() <
      p2.get_X() + 2*p2.get_Y() + 3*p2.get_Z();
}


Scene::Scene(log2plot::ConfigManager &_config):
  config(_config)
{
  if(config.read<bool>("kill"))
    system("killall python");

  scene_n = config.read<string>("scene");

  // read 3D points
  points.clear();
  const double d1 = read<double>("d1");
  const double d2 = read<double>("d2");
  const double d3 = read<double>("d3");
  const double d4 = read<double>("d4");
  const double d5 = read<double>("d5");
  points.push_back({0,0,0});
  points.push_back({1,0,0});
  points.push_back({d1,d2,0});
  points.push_back({d3,d4,d5});

  const double cx(1+d1+d3);
  const double cy(d2+d4);
  const double cz(d5);
  centroid.set(cx/4, cy/4, cz/4);

  // read singular point
  singular_point = read<vpTranslationVector>("singularity");

  n_points = static_cast<uint>(points.size());
}


std::pair<std::vector<vpHomogeneousMatrix>, vpHomogeneousMatrix>
Scene::vsPoses()
{
  std::vector<vpTranslationVector> oTc(1);
  vpTranslationVector oTcd;
  // read initial / final poses
  try {
    oTc[0] = read<vpTranslationVector>("oTc");
    oTcd = read<vpTranslationVector>("oTcd");
    std::cout << "Found initial / final translations\n";
  } catch (...) {

    // around singular point
    oTc[0] = oTcd = singular_point[0];
    const auto d = 0.7*config.read<double>("vs_dist");

    oTc[0][0] += d;
    oTc[0][1] += d;
    oTcd[0] -= d;
    oTcd[1] -= d;
  }

  // init VS desired pose
  vpHomogeneousMatrix cdMo(cMoFrom(oTcd));

  // read VS starting positions
  uint idx = 2;
  while(true)
  {
    std::stringstream ss;
    ss << "m" << idx;
    if(!config.has({scene_n, ss.str()}))
    {
      if(!config.has(ss.str()))
        break;
      oTc.push_back(singular_point + config.read<vpTranslationVector>(ss.str()));
    }
    else
    {
      oTc.push_back(singular_point + read<vpTranslationVector>(ss.str()));
    }
    idx++;
  }
  std::cout << "Found " << oTc.size() << " VS starts\n";

  std::vector<vpHomogeneousMatrix> cMo;
  for(auto &T: oTc)
    cMo.push_back(cMoFrom(T, cdMo));

  return {cMo, cdMo};
}

void Scene::displayPoints(log2plot::Logger &logger,
                          std::string point_color,
                          std::string singularity_color)
{
  // display points
  std::vector<std::vector<double>> Obj(n_points);
  for(uint i = 0; i < n_points; ++i)
    Obj[i] = {points[i].get_oX(), points[i].get_oY(), points[i].get_oZ()};

  logger.showFixedObject(Obj, log2plot::legendFullyConnected(n_points), point_color);

  // display singular points
  Obj = {{singular_point[0], singular_point[1], singular_point[2]}};

  logger.showFixedObject(Obj, log2plot::legendFullyDisconnected(1), singularity_color);
}


vpHomogeneousMatrix Scene::cMoFrom(const vpTranslationVector &oTc,
                               const vpHomogeneousMatrix &cMo_ref,
                               bool use_ref)
{
  enum class Phase{BRING_TO_FRONT,CENTER_IMAGE_POINTS,ALIGN_REF};

  vpHomogeneousMatrix oMc(oTc, vpRotationMatrix());


  vpMatrix L(n_points, 3), H(n_points,n_points);
  vpColVector s(n_points), sd(n_points);
  vpThetaUVector w;
  const double l = 0.1;
  double Zmin(0);
  vpRotationMatrix R;

  std::vector<Phase> phases{Phase::BRING_TO_FRONT,Phase::CENTER_IMAGE_POINTS};
  if(use_ref)
    phases.push_back(Phase::ALIGN_REF);

  for(auto phase: phases)
  {
    uint iter(0);

    if(phase == Phase::CENTER_IMAGE_POINTS)
    {
      L.resize(2*n_points, 3);
      H.eye(2*n_points);
      s.resize(2*n_points);
      sd.resize(2*n_points);
    }
    else if(phase == Phase::ALIGN_REF)
    {
      L.resize(n_points, 3);
      H.eye(n_points);
      s.resize(n_points);
      sd.resize(n_points);
      for(uint i = 0; i < n_points; ++i)
      {
        points[i].track(cMo_ref);
        sd[i] = atan2(points[i].get_y(), points[i].get_x());
        L[i][2] = -1;
      }
    }
    std::pair<double, uint> maxCoord{-1, 0};
    while(iter++ < 1000)
    {
      // get interaction matrix and weights
      for(uint i = 0; i < n_points; ++i)
      {
        points[i].track(oMc.inverse());

        const double X = points[i].get_X();
        const double Y = points[i].get_Y();
        const double Z = points[i].get_Z();
        const double x = X/Z;
        const double y = Y/Z;

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

      Zmin = std::min_element(points.begin(), points.end(),
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

      R.buildFrom(w);

      oMc.insert(oMc.getRotationMatrix()*R);
    }
  }

  return oMc.inverse();
}
