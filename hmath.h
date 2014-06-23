
#include <assert.h>
#include <numeric>
#include <algorithm>
#include <valarray>
#include <vector>

#include "matrix.h"

#ifndef _HMAP2_MATH
#define _HMAP2_MATH

using namespace std;

template <typename T>
inline T square(const T& d) { return d*d; };

template <typename T>
T dot_product (valarray<T>& v1, valarray<T>& v2)
{
  assert (v1.size()==v2.size());
  valarray<T> v3(v1.size());
  transform (&v1[0],&v1[v1.size()],&v2[0],&v3[0],multiplies<T>());
  return v3.sum();
}

template <typename T>
T norm_dot_product (valarray<T>& v1, valarray<T>& v2)
{
  assert (v1.size()==v2.size());
  valarray<T> v3(v1.size()), v4(v2.size());
  transform (&v1[0],&v1[v1.size()],&v2[0],&v3[0],multiplies<T>());
  T res = v3.sum();

  v3 = v1.apply(&square<T>);
  v4 = v2.apply(&square<T>);
  T norm = sqrt(v3.sum()) * sqrt(v4.sum());

  return res/norm;
}

template <typename T>
void norm_elements (valarray<T>& res, valarray<T>& v1)
{
  T sum = v1.sum();
  res = v1.apply(&square<T>);
  T sumsq = res.sum();

  T avg = sum/float(v1.size());
  T var = sumsq/float(v1.size())-avg*avg;
  T std = std::sqrt (var);

  assert (std!=0);
  res = v1;
  res -= avg;
  res /= std;
}

template <typename T>
void norm_elements (matrix<T>& res, matrix<T>& m, int i0, int i1, int j0, int j1)
{
  assert (res.rows()==m.rows() && res.cols()==m.cols());
  if (i0>=i1 || j0>=j1)
    { i0 = 0; j0 = 0; i1 = m.rows(); j1 = m.cols(); }
  valarray<T> v1((i1-i0)*(j1-j0)),v2((i1-i0)*(j1-j0));
  
  int c = 0;
  for (int i=i0; i<i1; ++i) 
    for (int j=j0; j<j1; ++j) v1[c++] = m(i,j);
  
  norm_elements (v2,v1);
  
  c=0;
  for (int i=i0; i<i1; ++i) 
    for (int j=j0; j<j1; ++j) res(i,j) = v2[c++];
}

template <typename T>
void shift_elements (matrix<T>& res, matrix<T>& m, int i0, int i1, int j0, int j1,
		     T shift)
{
  assert (res.rows()==m.rows() && res.cols()==m.cols());
  if (i0>=i1 || j0>=j1)
    { i0 = 0; j0 = 0; i1 = m.rows(); j1 = m.cols(); }

  for (int i=i0; i<i1; ++i) 
    for (int j=j0; j<j1; ++j) res(i,j) = m(i,j) + shift;
}

template <typename T>
T pearson_corr (valarray<T>& v1, valarray<T>& v2)
{
  assert (v1.size()==v2.size());
  valarray<T> n1(v1.size()), n2(v1.size());

  norm_elements<T> (n1,v1);
  norm_elements<T> (n2,v2);

  return dot_product (n1,n2) / v1.size();
}

#endif  // _HMAP2_MATH
