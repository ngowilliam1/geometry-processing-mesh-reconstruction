#include "fd_interpolate.h"
#include <cmath>

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////

  W.resize(P.rows(), nx*ny*nz);
  // best way according to Eigen Tutorial to fill a sparse matrix (W) is to build a list of triplets then convert it to a Sparsematrix
  // https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(P.rows()*8);

  // Obtaining the bottom-left-front corner position of grid in coordinate form so that we can later center at the origin
  double xCorner = corner(0);
  double yCorner = corner(1);
  double zCorner = corner(2);

  for  (int row=0; row<P.rows(); row++){
    // Obtain shifted coordinates of the points (shifted so leftmost point is centered at the origin)
    // Dividing by h makes the coordinates in units of h
    double x = (P(row, 0) - xCorner)/h;
    double y = (P(row, 1) - yCorner)/h;
    double z = (P(row, 2) - zCorner)/h;

    // x,y,z will be doubles, it is important to know which grid node they belong to (or rather are between), we'll take the floor for that
    int floorX = floor(x);
    int floorY = floor(y);
    int floorZ = floor(z);

    // Distance between the point and gridnode is also required to directly calculate weight
    double remainderX = x - floorX;
    double remainderY = y - floorY;
    double remainderZ = z - floorZ;


    // Computer the Weights for the 8 points (can look at this as all possible combination of a 3 variable 2 bit truth table)
    tripletList.push_back(T(row, floorX + floorY * nx + floorZ * (nx*ny), (1-remainderX)*(1-remainderY)*(1-remainderZ)));
    tripletList.push_back(T(row, floorX + floorY * nx + (floorZ+1) * (nx*ny), (1-remainderX)*(1-remainderY)*remainderZ));

    tripletList.push_back(T(row, floorX + (floorY+1) * nx + floorZ * (nx*ny), (1-remainderX)*remainderY*(1-remainderZ)));
    tripletList.push_back(T(row, floorX + (floorY+1) * nx + (floorZ+1) * (nx*ny), (1-remainderX)*remainderY*remainderZ));

    tripletList.push_back(T(row, (floorX+1) + floorY * nx + floorZ * (nx*ny), remainderX*(1-remainderY)*(1-remainderZ)));
    tripletList.push_back(T(row, (floorX+1) + floorY * nx + (floorZ+1) * (nx*ny), remainderX*(1-remainderY)*remainderZ));

    tripletList.push_back(T(row, (floorX+1) + (floorY+1) * nx + floorZ * (nx*ny), remainderX*remainderY*(1-remainderZ)));
    tripletList.push_back(T(row, (floorX+1) + (floorY+1) * nx + (floorZ+1) * (nx*ny), remainderX*remainderY*remainderZ));
  }

  W.setFromTriplets(tripletList.begin(),tripletList.end());
}
