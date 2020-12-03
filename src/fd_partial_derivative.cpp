#include "fd_partial_derivative.h"

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////

  // finding m (amount of rows of D)
  int m;
  int adjX = nx;
  int adjY = ny;
  int adjZ = nz;

  if (dir == 0){
    adjX--;
  } else if (dir == 1){
    adjY--;
  } else {
    adjZ--;
  }
  m = adjX*adjY*adjZ;
  D.resize(m, nx*ny*nz);

  // best way according to Eigen Tutorial to fill a sparse matrix (W) is to build a list of triplets then convert it to a Sparsematrix
  // https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletList;
  tripletList.reserve(D.rows()*2);

  for (int i = 0; i < adjX; i++){
    for (int j = 0; j < adjY; j++){
      for (int k = 0; k < adjZ; k++){
        int row = i + j*adjX + k*(adjX*adjY);
        int col1 = i + j*nx + k*(nx*ny);
        int col2 = i + j*nx + k*(nx*ny);
        if (dir == 0){
          col2 += 0;
        } else if (dir == 1){
          col2 += nx;
        } else {
          col2 += (nx*ny);
        }
        // Here I divided by 1/h to be consistent with course notes (making G a gradient matrix)
        // I think removing h will not change the look of the surface
        tripletList.push_back(T(row, col1, -1/h));
        tripletList.push_back(T(row, col2, 1/h));
      }
    }
  }
  D.setFromTriplets(tripletList.begin(),tripletList.end());
}
