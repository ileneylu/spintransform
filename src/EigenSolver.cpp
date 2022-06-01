// =============================================================================
// SpinXForm -- EigenSolver.cpp
// Keenan Crane
// August 16, 2011
//

#include "EigenSolver.h"
#include "LinearSolver.h"
#include <cmath>

void EigenSolver :: solve( QuaternionMatrix& A,
                           vector<Quaternion>& x )
// solves the eigenvalue problem Ax = cx for the
// eigenvector x with the smallest eigenvalue c
{
   // set the initial guess to the identity
   vector<Quaternion> b( x.size(), 1. );

   // perform a fixed number of inverse power iterations
   const int nIter = 3;
   for( int i = 0; i < nIter; i++ )
   {
      normalize( b );
      LinearSolver::solve( A, x, b, false );
      b = x;
   }

   // normalize the final solution
   normalize( x );
}

void EigenSolver :: solve( SparseMatrixd& A,
                           vector<Quaternion>& x )
// solves the eigenvalue problem Ax = cx for the
// eigenvector x with the smallest eigenvalue c
{
   // set the initial guess to the identity
   vector<Quaternion> b( x.size(), 1. );

   // perform a fixed number of inverse power iterations
   const int nIter = 3;
   for( int i = 0; i < nIter; i++ )
   {
      normalize( b );
      LinearSolver::solve( A, x, b, false );
      b = x;
   }

   // normalize the final solution
   normalize( x );
}

void EigenSolver :: normalize( vector<Quaternion>& x )
// rescales x to have unit length
{
   // compute length
   double norm = 0.;
   for( size_t i = 0; i < x.size(); i++ )
   {
      norm += x[i].norm2();
   }
   norm = sqrt( norm );

   // normalize
   for( size_t i = 0; i < x.size(); i++ )
   {
      x[i] /= norm;
   }
}

void EigenSolver :: toReal( const vector<Quaternion>& uQuat,
                             vector<double>& uReal )
// converts vector from quaternion- to real-valued entries
{
   for( size_t i = 0; i < uQuat.size(); i++ )
   {
      uReal[i*4+0] = uQuat[i].re();   // real
      uReal[i*4+1] = uQuat[i].im().x; // i
      uReal[i*4+2] = uQuat[i].im().y; // j
      uReal[i*4+3] = uQuat[i].im().z; // k
   }
}

