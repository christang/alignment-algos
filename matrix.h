/**
 *  File: matrix.h
 *  Desc: Based on Bjarne Stroustrup's 2-dim array using valarray + slice.
 *        (See also 22.4.6 in C++PL).
 *  
 *  Note: matrix add/sub appear to work properly, though w/o error checking
 *
 *  Created on 11/11/06
 *  Author: cltang @ honig lab { derived from BS's C++PL ed.3 }
 *
 */

#ifndef _MATRIX_BS
#define _MATRIX_BS

#include <valarray>

using namespace std;

// forward declarations to allow friend declarations:
template<class T> class Slice_iter;
template<class T> bool operator==(const Slice_iter<T>&, const Slice_iter<T>&);
template<class T> bool operator!=(const Slice_iter<T>&, const Slice_iter<T>&);
template<class T> bool operator< (const Slice_iter<T>&, const Slice_iter<T>&);

template<class T> class Slice_iter {
    valarray<T>* v;
    slice s;
    size_t curr;    // index of current element

    T& ref(size_t i) const { return (*v)[s.start()+i*s.stride()]; }
public:
    Slice_iter(valarray<T>* vv, slice ss) :v(vv), s(ss), curr(0) { }

    Slice_iter end() const
    {
        Slice_iter t = *this;
        t.curr = s.size();  // index of last-plus-one element
        return t;
    }

    Slice_iter& operator++() { curr++; return *this; }
    Slice_iter operator++(int) { Slice_iter t = *this; curr++; return t; }

    T& operator[](size_t i) { return ref(i); }      // C style subscript
    T& operator()(size_t i) { return ref(i); }      // Fortran-style subscript
    T& operator*() { return ref(curr); }            // current element

    friend bool operator==<>(const Slice_iter& p, const Slice_iter& q);
    friend bool operator!=<>(const Slice_iter& p, const Slice_iter& q);
    friend bool operator< <>(const Slice_iter& p, const Slice_iter& q);

};


template<class T>
bool operator==(const Slice_iter<T>& p, const Slice_iter<T>& q)
{
    return p.curr==q.curr
        && p.s.stride()==q.s.stride()
        && p.s.start()==q.s.start();
}

template<class T>
bool operator!=(const Slice_iter<T>& p, const Slice_iter<T>& q)
{
    return !(p==q);
}

template<class T>
bool operator<(const Slice_iter<T>& p, const Slice_iter<T>& q)
{
    return p.curr<q.curr
        && p.s.stride()==q.s.stride()
        && p.s.start()==q.s.start();
}


//-------------------------------------------------------------


// forward declarations to allow friend declarations:
template<class T> class Cslice_iter;
template<class T> bool operator==(const Cslice_iter<T>&, 
				  const Cslice_iter<T>&);
template<class T> bool operator!=(const Cslice_iter<T>&, 
				  const Cslice_iter<T>&);
template<class T> bool operator< (const Cslice_iter<T>&, 
				  const Cslice_iter<T>&);


template<class T> class Cslice_iter 
{
    valarray<T>* v;
    slice s;
    size_t curr; // index of current element
    const T& ref(size_t i) const { return (*v)[s.start()+i*s.stride()]; }
public:
    Cslice_iter(valarray<T>* vv, slice ss): v(vv), s(ss), curr(0){}
    Cslice_iter end() const
    {
        Cslice_iter t = *this;
        t.curr = s.size(); // index of one plus last element
        return t;
    }
    Cslice_iter& operator++() { curr++; return *this; }
    Cslice_iter operator++(int) { Cslice_iter t = *this; curr++; return t; }
    
    const T& operator[](size_t i) const { return ref(i); }
    const T& operator()(size_t i) const { return ref(i); }
    const T& operator*() const { return ref(curr); }

    friend bool operator==<>(const Cslice_iter& p, const Cslice_iter& q);
    friend bool operator!=<>(const Cslice_iter& p, const Cslice_iter& q);
    friend bool operator< <>(const Cslice_iter& p, const Cslice_iter& q);

};

template<class T>
bool operator==(const Cslice_iter<T>& p, const Cslice_iter<T>& q)
{
    return p.curr==q.curr
        && p.s.stride()==q.s.stride()
        && p.s.start()==q.s.start();
}

template<class T>
bool operator!=(const Cslice_iter<T>& p, const Cslice_iter<T>& q)
{
    return !(p==q);
}

template<class T>
bool operator<(const Cslice_iter<T>& p, const Cslice_iter<T>& q)
{
    return p.curr<q.curr
        && p.s.stride()==q.s.stride()
        && p.s.start()==q.s.start();
}


//-------------------------------------------------------------


template <class val_t>
class matrix : public valarray<val_t>
{

  typedef valarray<val_t> Base;

public:
  
  matrix (int nrows, int ncols);
  matrix (const matrix& m);
  matrix& operator= (const matrix& m);

  int size() const { return nrows*ncols; }
  int rows() const { return nrows; }
  int cols() const { return ncols; }

  val_t operator() (int r, int c) const;
  val_t& operator() (int r, int c);

  Slice_iter<val_t> row(int i);
  Cslice_iter<val_t> row(int i) const;
  Slice_iter<val_t> col(int i);
  Cslice_iter<val_t> col(int i) const;

  // Overriding valarray []...
  Slice_iter<val_t> operator[] (int i);
  Cslice_iter<val_t> operator[] (int i) const;

private:

  int nrows;
  int ncols;

};


template <class val_t>
matrix<val_t>::matrix(int nr, int nc)
  : valarray<val_t> (nr*nc),
    nrows           (nr),
    ncols           (nc) {}

template <class val_t>
matrix<val_t>::matrix(const matrix<val_t>& m)
  : valarray<val_t> (m),
    nrows           (m.nrows),
    ncols           (m.ncols) {}

template <class val_t>
matrix<val_t>& matrix<val_t>::operator= (const matrix<val_t>& m)
{ 
  if (&m == this) return *this;
  Base::operator=(m);
  nrows = m.nrows;
  ncols = m.ncols;
  return *this;
}

template <class val_t>
inline val_t matrix<val_t>::operator() (int r, int c) const 
{ return Base::operator[](r*ncols+c); }

template <class val_t>
inline val_t& matrix<val_t>::operator() (int r, int c)
{ return Base::operator[](r*ncols+c); }

template <class val_t>
inline Slice_iter<val_t> matrix<val_t>::row (int i) 
{ return Slice_iter<val_t>(this,slice(i*ncols,ncols,1)); }

template <class val_t>
inline Cslice_iter<val_t> matrix<val_t>::row (int i) const
{ return Cslice_iter<val_t>(this,slice(i*ncols,ncols,1)); }

template <class val_t>
inline Slice_iter<val_t> matrix<val_t>::col (int i) 
{ return Slice_iter<val_t>(this,slice(i,nrows,ncols)); }

template <class val_t>
inline Cslice_iter<val_t> matrix<val_t>::col (int i) const 
{ return Cslice_iter<val_t>(this,slice(i,nrows,ncols)); }

template <class val_t>
inline Slice_iter<val_t> matrix<val_t>::operator[] (int i)
{ return row(i); }

template <class val_t>
inline Cslice_iter<val_t> matrix<val_t>::operator[] (int i) const
{ return row(i); }

#endif  // _MATRIX_BS
