#include <log2plot/logger.h>
#include <task.h>
#include <scene.h>
#include <pose_estim.h>
#include <visp/vpRobotCamera.h>
#include <visp/vpQuadProg.h>
#include <numeric>


constexpr double side = 2;
constexpr int steps = 201;  // between -side and +side
constexpr double step = (2*side)/(steps-1);
constexpr double inf = std::numeric_limits<double>::max();

constexpr double Zoffset = 0.00;

struct Cell
{
  double err = inf;
  vpQuaternionVector q;
  static bool isEmpty(const Cell &cell)
  {
    return cell.err == inf;
  }
};

struct Grid
{
  std::array<Cell, steps*steps> grid;
  const vpTranslationVector &singular_point;
public:
  Grid(const vpTranslationVector &point) : singular_point(point){}

  void printPercentage() const
  {
    std::cout << std::fixed << std::setprecision(2)
              << 100-(100.*std::count_if(grid.begin(), grid.end(), Cell::isEmpty))/grid.size()
              << " % " << std::endl;
  }

  vpHomogeneousMatrix averageNeighboors(std::array<const Cell, steps*steps>::iterator cell) const
  {
    // get this cell coordinates
    int idx = std::distance(grid.begin(), cell);

    vpHomogeneousMatrix cMo;
    // find neighboors with defined rotations
    /*vpMatrix Q(4, 0);
    for(int dir: {-1, 1, -steps, steps})
    {
      if(idx + dir >= 0 && idx + dir < grid.size() && !Cell::isEmpty(grid[idx+dir]))
      {
        Q.resize(4, Q.getCols()+1,false);
        for(uint i = 0; i < 4; ++i)
          Q[i][Q.getCols()-1] = grid[idx+dir].q[i];
      }
    }

    Q = Q*Q.t();
    vpColVector eigenVals;
    vpMatrix eigenVectors;
    Q.eigenValues(eigenVals, eigenVectors);
    const double eigenMax(eigenVals.getMaxValue());


    for(uint i = 0; i < eigenVals.size(); ++i)
    {
      if(eigenVals[i] == eigenMax)
      {
        cMo.insert(vpQuaternionVector{Q[0][i], Q[1][i], Q[2][i], Q[3][i]}/Q.getCol(i).frobeniusNorm());
        break;
      }
    }*/
    for(int dir: {-1, 1, -steps, steps})
    {
      if(idx + dir >= 0 && idx + dir < grid.size() && !Cell::isEmpty(grid[idx+dir]))
      {
        cMo.insert(grid[idx+dir].q);
        break;
      }
    }

    cMo.insert(-(cMo.getRotationMatrix()*
                 (singular_point+vpTranslationVector{value(idx % steps), value(idx/steps), Zoffset})));
    return cMo;
  }

  static inline size_t index(double v)
  {
    return std::round((v+side)/step);
  }
  static inline double value(size_t i)
  {
    return -side + i*step;
  }
  Cell &at(size_t xi, size_t yi)
  {
    return grid[xi+yi*steps];
  }
  void update(double x, double y, double val, const vpHomogeneousMatrix &cMo)
  {
    auto &cur(at(index(x), index(y)));
    if(val < cur.err)
    {
      cur.err = val;
      cur.q.buildFrom(cMo.getRotationMatrix());
    }
  }
};

int main()
{
  // init and read config file
  log2plot::ConfigManager config(std::string(SRC_PATH) + "/config.yaml");
  config.forceParameter("control", "none");
  config.forceParameter("estim", "none");

  config.forceParameter("vs_dist", "0");

  const bool singularScene(true);
  const auto base_scene(config.read<std::string>("scene").substr(0,1));
  if(singularScene)
    config.forceParameter("scene", base_scene);
  else
    config.forceParameter("scene", base_scene + "o");

  Scene scene(config);
  PoseEstim estimator(scene, "energy"); // generate dest path

  const auto [cMos, cdMo] = scene.vsPoses(); {}
  const auto cdMo2(cdMo);
  //const vpHomogeneousMatrix cdMo(cdMop);
  const auto cdRc(cdMo.getRotationMatrix());

  Task task(scene.points, cdMo);

  Grid energy(scene.singular_point);
  vpColVector v(6);
  vpMatrix A(1, 6);
  const vpColVector b(1,0);

  vpRobotCamera robot;
/*
    // bad temptative to perform exhaustive VS within the square
    // nice figures though, but impossible to understand
  const auto fillFrom = [&](vpHomogeneousMatrix cMo)
  {
    robot.setPosition(cMo);

    while(true)
    {
      robot.getPosition(cMo);
      try {
        task.computeFrom(cMo);
      } catch (...) {
        return;
      }


      // write energy at this position
      const auto oTs(-(cMo.getRotationMatrix().t()*cMo.getTranslationVector())-scene.singular_point);

      energy.update(oTs[0], oTs[1], task.e.frobeniusNorm(), cMo);

      // keep 0 velocity on zo
      for(int i = 0; i < 3; ++i)
        A[0][i] = cMo[i][2];
      vpQuadProg::solveQPe(task.L, -task.e, A, b, v);
      const double norm(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
      //std::cout << " norm = " << norm << std::endl;
      if(norm > 2*step)
        v *= 2*step/norm;
      else if(norm < .1*step)
        break;

      robot.setVelocity(vpRobot::CAMERA_FRAME, v);
    }
    energy.printPercentage();
  };

  // start from border
  for(size_t xi = 0; xi < steps; xi++)
  {
    const bool border(xi == 0 || xi ==steps-1);
    double x(energy.value(xi));

    for(size_t yi = 0; yi < steps; yi += border ? 1 : steps-1)
    {
      double y(energy.value(yi));
      // init a VS from here
      fillFrom(scene.cMoFrom(scene.singular_point + vpTranslationVector{x, y, Zoffset}, cdMo));
  }
  }

  // fill up empty cells
  auto empty = std::find_if(energy.grid.begin(), energy.grid.end(), Cell::isEmpty);
  while(empty != energy.grid.end())
  {
    fillFrom(energy.averageNeighboors(empty));
    empty = std::find_if(energy.grid.begin(), energy.grid.end(), Cell::isEmpty);
  }

  */

  vpRotationMatrix Res(0, 0, 0.4);

  const auto fillFrom = [&](vpTranslationVector oTc)
  {
    vpHomogeneousMatrix cMo(-(Res*cdRc*oTc), Res*cdRc);
    //vpHomogeneousMatrix cMo(scene.cMoFrom(oTc, cdMo2));

    //cMo.insert(cdRc);
    task.computeFrom(cMo);
    const auto oTs(-(cMo.getRotationMatrix().t()*cMo.getTranslationVector())-scene.singular_point);
    energy.update(oTs[0], oTs[1], task.e.frobeniusNorm(), cMo);
  };
  // do all cells
  for(size_t xi = 0; xi < steps; xi++)
  {
    double x(energy.value(xi));
    energy.printPercentage();
    for(size_t yi = 0; yi < steps; yi++)
    {
      double y(energy.value(yi));
      // init a VS from here
      fillFrom(scene.singular_point + vpTranslationVector{x, y, Zoffset});
    }
  }

  // energy is visual error at (x,y)
  log2plot::Logger logger(config.fullName());
  std::array<double, steps> e;

  std::string dest("eR_" + std::to_string((int)side) + "_" + config.read<std::string>("vs_dist"));

  logger.save(e, dest, "[]", "");
  // desired position wrt singular point
  vpTranslationVector sTcd(cdMo.inverse());
  sTcd = sTcd - scene.singular_point;
  logger.showFixedObject({{sTcd[0], sTcd[1]}},
                         "[[0,0]]", "gD");
  logger.setPlotArgs(std::to_string(side));

  for(int xi = 0; xi < steps; ++xi)
  {
    for(int yi = 0; yi < steps; ++yi)
    {
      e[yi] = energy.at(xi, yi).err;
      if(e[yi] == std::numeric_limits<double>::max())
        log2plot::setNaN(e, yi, yi+1);
    }
    logger.update();
  }
  std::cout << "Saving to " << config.fullName() +
               dest << ".yaml\n";
}
