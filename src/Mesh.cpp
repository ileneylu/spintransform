// ============================================================================
// SpinXForm -- Mesh.cpp
// Keenan Crane
// August 16, 2011
//

#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <set>
#include <map>
#include "Mesh.h"
#include "LinearSolver.h"
#include "EigenSolver.h"
#include "Utility.h"

#include <igl/readOBJ.h>
#include <igl/boundary_loop.h>
#include <Eigen/Core>

void Mesh :: updateDeformation( void )
{
   int t0 = clock();

   // solve eigenvalue problem for local similarity transformation lambda
   buildEigenvalueProblem();
   EigenSolver::solve( E, lambda );

   // solve Poisson problem for new vertex positions
   buildPoissonProblem();
   LinearSolver::solve( L, newVertices, omega );
   normalizeSolution();

   int t1 = clock();
   cout << "time: " << (t1-t0)/(double) CLOCKS_PER_SEC << "s" << endl;
}

void Mesh :: resetDeformation( void )
{
   // copy original mesh vertices to current mesh
   for( size_t i = 0; i < vertices.size(); i++ )
   {
      newVertices[i] = vertices[i];
   }

   normalizeSolution();
}

void Mesh :: setCurvatureChange( const Image& image, const double scale )
// sets rho values by interpreting "image" as a square image
// in the range [0,1] x [0,1] and mapping values to the
// surface via vertex texture coordinates -- grayscale
// values in the  range [0,1] get mapped (linearly) to values
// in the range [-scale,scale]
{
   double w = (double) image.width();
   double h = (double) image.height();

   for( size_t i = 0; i < faces.size(); i++ )
   {
      // get handle to current value of rho
      rho[i] = 0.;

      // compute average value over the face
      for( int j = 0; j < 3; j++ )
      {
         Vector uv = faces[i].uv[j];
         rho[i] += image.sample( uv.x*w, uv.y*h ) / 3.;
      }

      // map value to range [-scale,scale]
      rho[i] = (2.*(rho[i]-.5)) * scale;
   }
}

void Mesh :: initBoundaryTangent()
{
   t.resize(numBoundaryVertices, 3);
   tTilda.resize(numBoundaryVertices, 3);
   int curRow = 0;
   for (int i = 0; i < boundary_loop.size(); i++) 
   {
      int boundaryLength = boundary_loop[i].size();
      for (int j = 0; j < boundaryLength; j++) 
      {
         int vNext = boundary_loop[i][(j+1)%boundaryLength];
         int vPrev = boundary_loop[i][(j-1)%boundaryLength];
         Quaternion diff = vertices[vNext] - vertices[vPrev];
         Eigen::RowVectorXd row(3); 
         row << (diff.im())[0],(diff.im())[1],(diff.im())[2];
         t.row(curRow) = row;
         tTilda.row(curRow) = row;
         curRow++;
      }
   }
}

void Mesh :: setBoundaryTangent()
{
   initBoundaryTangent();
   int nV = vertices.size();
   int nB = t.rows();
   int nI = nV - nB;

   // Build W
   QuaternionMatrix W;
   W.resize(nV,nV);

   for (int i = 0; i < nI; i++) 
   {
      W(i,i) = Quaternion(1.,0.,0.,0.);
   }
   for (int i = 0; i < nB; i++) 
   {
      double theta;
      Quaternion w;
      if (((tTilda.row(i) + t.row(i)).norm()< 0.00001) ||
          ((tTilda.row(i) - t.row(i)).norm()< 0.00001)) 
      {
         // tTilda = +- t
         theta = 0;
         w = Quaternion(0.,1.,0.,0.);
      } else {
         double theta = acos(t.row(i).dot(tTilda.row(i)));
         Eigen::RowVector3d ti = t.row(i);
         Eigen::RowVector3d tTildai = tTilda.row(i);
         Eigen::Vector3d wVec = ti.cross(tTildai/sin(theta));
         w = Quaternion(0.,wVec[0],wVec[1],wVec[2]);
      }

      W(i+nI,i+nI) = Quaternion(cos(theta/2),0.,0.,0.) + sin(theta/2) * w;
   }

   // Build C
   Eigen::MatrixXd C = Eigen::MatrixXd::Zero(4*numBoundaryVertices, 2*numBoundaryVertices);
   for (int i = 0; i < tTilda.rows(); i++) {
      C(4*i,2*i) = 1;
      C(4*i+1,2*i+1) = tTilda(i,0);
      C(4*i+2,2*i+1) = tTilda(i,1);
      C(4*i+3,2*i+1) = tTilda(i,2);
   }

   // Build IC
   int numInternalVertices = vertices.size()-numBoundaryVertices;
   Eigen::MatrixXd I = Eigen::MatrixXd::Identity(4*numInternalVertices,4*numInternalVertices);
   Eigen::MatrixXd IC = Eigen::MatrixXd::Zero(4*vertices.size(), 4*numInternalVertices+2*numBoundaryVertices);
   IC.topLeftCorner(4*numInternalVertices,4*numInternalVertices) = I;
   IC.bottomRightCorner(4*numBoundaryVertices, 2*numBoundaryVertices) = C;

   // build U
   Eigen::MatrixXd U = W.toEigenReal() * IC;
}

double Mesh :: area( int i )
// returns area of triangle i in the original mesh
{
   Vector& p1 = vertices[ faces[i].vertex[0] ].im();
   Vector& p2 = vertices[ faces[i].vertex[1] ].im();
   Vector& p3 = vertices[ faces[i].vertex[2] ].im();

   return .5 * (( p2-p1 ) ^ ( p3-p1 )).norm();
}

void Mesh :: buildEigenvalueProblem( void )
{
   // allocate a sparse |V|x|V| matrix
   int nV = vertices.size();
   E.resize( nV, nV );

   // visit each face
   for( size_t k = 0; k < faces.size(); k++ )
   {
      double A = area(k);
      double a = -1. / (4.*A);
      double b = rho[k] / 6.;
      double c = A*rho[k]*rho[k] / 9.;

      // get vertex indices
      int I[3] =
      {
         faces[k].vertex[0],
         faces[k].vertex[1],
         faces[k].vertex[2]
      };

      // compute edges across from each vertex
      Quaternion e[3];
      for( int i = 0; i < 3; i++ )
      {
         e[i] = vertices[ I[ (i+2) % 3 ]] -
                vertices[ I[ (i+1) % 3 ]] ;
      }

      // increment matrix entry for each ordered pair of vertices
      for( int i = 0; i < 3; i++ )
      for( int j = 0; j < 3; j++ )
      {
         E(I[i],I[j]) += a*e[i]*e[j] + b*(e[j]-e[i]) + c;
      }
   }
}

void Mesh :: buildPoissonProblem( void )
{
   buildLaplacian();
   buildOmega();
}

void Mesh :: buildLaplacian( void )
// builds the cotan-Laplace operator
{
   // allocate a sparse |V|x|V| matrix
   int nV = vertices.size();
   L.resize( nV, nV );

   // visit each face
   for( size_t i = 0; i < faces.size(); i++ )
   {
      // visit each triangle corner
      for( int j = 0; j < 3; j++ )
      {
         // get vertex indices
         int k0 = faces[i].vertex[ (j+0) % 3 ];
         int k1 = faces[i].vertex[ (j+1) % 3 ];
         int k2 = faces[i].vertex[ (j+2) % 3 ];

         // get vertex positions
         Vector f0 = vertices[k0].im();
         Vector f1 = vertices[k1].im();
         Vector f2 = vertices[k2].im();

         // compute cotangent of the angle at the current vertex
         // (equal to cosine over sine, which equals the dot
         // product over the norm of the cross product)
         Vector u1 = f1 - f0;
         Vector u2 = f2 - f0;
         double cotAlpha = (u1*u2)/(u1^u2).norm();

         // add contribution of this cotangent to the matrix
         L( k1, k2 ) -= cotAlpha / 2.;
         L( k2, k1 ) -= cotAlpha / 2.;
         L( k1, k1 ) += cotAlpha / 2.;
         L( k2, k2 ) += cotAlpha / 2.;
      }
   }
}

void Mesh :: buildOmega( void )
{
   // clear omega
   for( size_t i = 0; i < omega.size(); i++ )
   {
      omega[i] = 0.;
   }

   // visit each face
   for( size_t i = 0; i < faces.size(); i++ )
   {
      // get indices of the vertices of this face
      int v[3] = { faces[i].vertex[0],
                   faces[i].vertex[1],
                   faces[i].vertex[2] };

      // visit each edge
      for( int j = 0; j < 3; j++ )
      {
         // get vertices
         Quaternion f0 = vertices[ v[ (j+0) % 3 ]];
         Quaternion f1 = vertices[ v[ (j+1) % 3 ]];
         Quaternion f2 = vertices[ v[ (j+2) % 3 ]];

         // determine orientation of this edge
         int a = v[ (j+1) % 3 ];
         int b = v[ (j+2) % 3 ];
         if( a > b )
         {
            swap( a, b );
         }

         // compute transformed edge vector
         Quaternion lambda1 = lambda[a];
         Quaternion lambda2 = lambda[b];
         Quaternion e = vertices[b] - vertices[a];
         Quaternion eTilde = (1./3.) * (~lambda1) * e * lambda1 +
                             (1./6.) * (~lambda1) * e * lambda2 +
                             (1./6.) * (~lambda2) * e * lambda1 +
                             (1./3.) * (~lambda2) * e * lambda2 ;

         // compute cotangent of the angle opposite the current edge
         Vector u1 = ( f1 - f0 ).im();
         Vector u2 = ( f2 - f0 ).im();
         double cotAlpha = (u1*u2)/(u1^u2).norm();

         // add contribution of this edge to the divergence at its vertices
         omega[a] -= cotAlpha * eTilde / 2.;
         omega[b] += cotAlpha * eTilde / 2.;
      }
   }

   removeMean( omega );
}

void Mesh :: normalizeSolution( void )
{
   // center vertices around the origin
   removeMean( newVertices );

   // find the vertex with the largest norm
   double r = 0.;
   for( size_t i = 0; i < vertices.size(); i++ )
   {
      r = max( r, newVertices[i].norm2() );
   }
   r = sqrt(r);

   // rescale so that vertices have norm at most one
   for( size_t i = 0; i < vertices.size(); i++ )
   {
      newVertices[i] /= r;
   }
}


// FILE I/O --------------------------------------------------------------------

void Mesh :: read( const string& filename )
// loads a triangle mesh in Wavefront OBJ format
{
   // set default to boundary free
   hasBoundary = false;
   
   // open mesh file
   ifstream in( filename.c_str() );
   if( !in.is_open() )
   {
      cerr << "Error: couldn't open file ";
      cerr << filename;
      cerr << " for input!" << endl;
      exit( 1 );
   }

   // temporary list of vertex coordinates
   vector<vector<double>> V, TC;
   vector<vector<int>> F, FTC;

   // parse mesh file
   string s;
   while( getline( in, s ))
   {
      stringstream line( s );
      string token;

      line >> token;

      if( token == "v" ) // vertex
      {
         double x, y, z;

         line >> x >> y >> z;

         vector<double> v = {x,y,z};
         V.emplace_back(v);
      }
      if( token == "vt" ) // texture coordinate
      {
         double u, v;

         line >> u >> v;

         vector<double> tc = {u,v};
         TC.emplace_back(tc);
      }
      else if( token == "f" ) // face
      {
         // Face triangle;
         vector<int> fv,ftc;

         // iterate over vertices
         for( int i = 0; i < 3; i++ )
         {
            line >> s;
            stringstream item( s );

            int I[3] = { -1, -1, -1 };

            // iterate over v, vt, and vn indices
            for( int j = 0; getline( item, s, '/' ) && j < 3; j++ )
            {
               stringstream index( s );
               index >> I[j];
            }

            fv.emplace_back(I[0]-1);

            if( I[1] != -1 )
            {
               ftc.emplace_back(I[1]-1);
            }
         }
         F.emplace_back(fv);
         FTC.emplace_back(ftc);
      }
   }

   // Construct F as Eigen Matrix
   Eigen::MatrixXi FMat(F.size(),F[0].size());
   for (int i = 0; i < F.size(); i++) {
      Eigen::RowVectorXi f(3);
      f << F[i][0],F[i][1],F[i][2];
      FMat.row(i) = f;
   }

   // boundary vertices
   set<int> boundary_vertices;
   numBoundaryVertices = 0;
   igl::boundary_loop(FMat, boundary_loop);
   for (int i = 0; i < boundary_loop.size(); i++) {
      numBoundaryVertices += boundary_loop[i].size();
      for (int j = 0; j < boundary_loop[i].size(); j++) {
         boundary_vertices.emplace(boundary_loop[i][j]);
      }
   }

   if (!boundary_vertices.empty()) {
      hasBoundary = true;
      // Re-indexing V with internal vertices followed by boundary vertices
      map<int,int> vertex_old_new_index_map;
      vector<vector<double>> internal_v, boundary_v;
      for (int i = 0; i < V.size(); i++) {
         if (boundary_vertices.find(i) != boundary_vertices.end()) {
            // boundary vertex
            vertex_old_new_index_map[i] = boundary_v.size();
            boundary_v.emplace_back(V[i]);
         } else {
            // internal vertex
            vertex_old_new_index_map[i] = internal_v.size();
            internal_v.emplace_back(V[i]);
         }
      }
      V.clear();
      V.reserve(internal_v.size() + boundary_v.size());
      V.insert(V.end(), internal_v.begin(), internal_v.end());
      V.insert(V.end(), boundary_v.begin(), boundary_v.end() );
      
      // update vertex_old_new_index_map
      int num_internal_vertices = internal_v.size();
      for (int vi : boundary_vertices) {
         vertex_old_new_index_map[vi] += num_internal_vertices;
      }

      // update F
      for (int i = 0; i < F.size(); i++) {
         for (int j = 0; j < F[i].size(); j++) {
            F[i][j] = vertex_old_new_index_map[F[i][j]];
         }
      }
   }

   for (int i = 0; i < V.size(); i++) {
      vertices.emplace_back(Quaternion(0.,V[i][0],V[i][1],V[i][2]));
      newVertices.emplace_back(Quaternion(0.,V[i][0],V[i][1],V[i][2]));
   }

   for (int i = 0; i < F.size(); i++) {
      Face triangle;
      for (int j = 0; j < 3; j++) {
         triangle.vertex[j] = F[i][j];
         vector<double> tc = TC[FTC[i][j]];
         triangle.uv[j] = Vector(tc[0],tc[1],0.);
      }
      faces.emplace_back(triangle);
   }

   // allocate space for mesh attributes
   lambda.resize( vertices.size() );
   omega.resize( vertices.size() );
   rho.resize( faces.size() );
   normalizeSolution();
}

void Mesh :: write( const string& filename )
// saves a triangle mesh in Wavefront OBJ format
{
   ofstream out( filename.c_str() );

   if( !out.is_open() )
   {
      cerr << "Error: couldn't open file ";
      cerr << filename;
      cerr << " for output!" << endl;
      return;
   }

   for( size_t i = 0; i < vertices.size(); i++ )
   {
      out << "v " << newVertices[i].im().x << " "
                  << newVertices[i].im().y << " "
                  << newVertices[i].im().z << endl;
   }

   for( size_t i = 0; i < faces.size(); i++ )
   {
      out << "f " << 1+faces[i].vertex[0] << " "
                  << 1+faces[i].vertex[1] << " "
                  << 1+faces[i].vertex[2] << endl;
   }
}

