// ============================================================================
// SpinXForm -- Mesh.h
// Keenan Crane
// August 16, 2011
//
// Mesh is a very simple mesh data structure consisting of a list of vertices
// and a collection of triangles specified via indices into the vertex list.
// It also stores the data necessary to compute a spin transformation of the
// surface.
//
// To specify a deformation, the user should set a value of "rho" on each
// face corresponding to the desired change in curvature.  The deformation
// is computed by calling updateDeformation(), which puts the transformed
// vertices in the list "newVertices."
//

#ifndef SPINXFORM_MESH_H
#define SPINXFORM_MESH_H

#include <vector>
#include <string>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "Quaternion.h"
#include "QuaternionMatrix.h"
#include "Image.h"

using namespace std;

class Face
{
   public:
      int vertex[3]; // indices into vertex list
      Vector uv[3]; // texture coordinates (for visualization only)
};

class Mesh
{
   public:
      void read( const string& filename );
      // loads a triangle mesh in Wavefront OBJ format

      void write( const string& filename );
      // saves a triangle mesh in Wavefront OBJ format

      void readCurvatureChange( const string& filename);

      void setCurvatureChange( const Image& image, const double scale );
      // sets rho values by interpreting "image" as a square image
      // in the range [0,1] x [0,1] and mapping values to the
      // surface via vertex texture coordinates -- grayscale
      // values in the  range [0,1] get mapped (linearly) to values
      // in the range [-scale,scale]

      void setCurvatureChange( );

      void setBoundaryCondition( void );
      // setup matrix U that handles boundary condition

      void updateDeformation( void );
      // computes a conformal deformation using the current rho

      void resetDeformation( void );
      // restores surface to its original configuration

      double area( int i );
      // returns area of triangle i in the original mesh

      bool hasBoundary;
      // indicates whether the mesh has boundaries

      vector<Face> faces;
      // list of triangles as indices into vertex list

      vector<Quaternion> vertices, newVertices;
      // original and deformed vertex coordinates

      vector<double> rho;
      // controls change in curvature (one value per face)

      vector<vector<int>> boundary_loop;
      int numBoundaryVertices;
      Eigen::MatrixXd t;
      Eigen::MatrixXd tTilda;

   protected:

      vector<Quaternion> lambda;
      // local similarity transformation (one value per vertex)
      vector<double> a;
      vector<double> b;
      // local tangent directions (one pair per boundary vertex)

      vector<Quaternion> omega;
      // divergence of target edge vectors

      QuaternionMatrix L; // Laplace matrix
      QuaternionMatrix E; // matrix for eigenvalue problem

      Eigen::SparseMatrix<double> U; // Boundary matrix

      map<int,int> vertex_old_new_index_map;
      vector<double> rhoV;

      void buildEigenvalueProblem( void );
      void buildPoissonProblem( void );
      void buildLaplacian( void );
      void buildOmega( void );
      void normalizeSolution( void );

      // temporary function to handle boundary case
      void initBoundaryTangent();

      void toReal( const vector<Quaternion>& uQuat,
                             vector<double>& uReal );

      void toQuat( const vector<double>& uReal,
                             vector<Quaternion>& uQuat );
};

#endif
