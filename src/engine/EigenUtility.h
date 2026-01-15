#ifndef EIGEN_UTILITY_H
#define EIGEN_UTILITY_H

#include <Eigen/Dense>

template <typename Derived, typename ftype>
auto ClipEigen(const Eigen::ArrayBase<Derived>& array,
          const ftype& min,
          const ftype& max)
{
  return array.cwiseMin(max).cwiseMax(min);
}

template <typename Derived>
auto SafeSetEigen(Eigen::ArrayBase<Derived>& array,
          const Eigen::ArrayBase<Derived>& array2)
{
  int minDim = std::min(array.size(), array2.size());
  array.head(minDim) = array2.head(minDim);
}

#endif // EIGEN_UTILITY_H