// =============================================================================
// SpinXForm -- main.cpp
// Keenan Crane
// August 16, 2011
//

#include <iostream>
#include <string.h>
#include "Viewer.h"

using namespace std;

int main( int argc, char **argv )
{

   if( argc < 3 || argc > 5 )
   {
      cerr << "usage: " << argv[0] << " mesh.obj mode image.tga [result.obj]" << endl;
      return 1;
   }

   if( argc == 5 ) // batch mode
   {
      // load mesh
      Mesh mesh;
      mesh.read( argv[2] );

      int mode = std::stoi(argv[1]);
      if (mode == 0) {
      // load image
         Image image;
         image.read( argv[3] );
         const double scale = 5.;
         mesh.setCurvatureChange( image, scale );
      } else if (mode == 1) {
         mesh.readCurvatureChange( argv[3] );
         mesh.setCurvatureChange();
      }

      // apply transformation
      mesh.updateDeformation();

      // write result
      mesh.write( argv[4] );
   }
   else // interactive mode
   {
      Viewer viewer;

      // load mesh
      viewer.mesh.read( argv[2] );

      int mode = std::stoi(argv[1]);
      if (mode == 0) {
      // load image
         // Image image;
         // image.read( argv[3] );
         viewer.image.read( argv[3] );
         const double scale = 5.;
         viewer.mesh.setCurvatureChange( viewer.image, scale );
      } else if (mode == 1) {
         viewer.mesh.readCurvatureChange( argv[3] );
         viewer.mesh.setCurvatureChange();
      }

      // load image
      // viewer.image.read( argv[2] );

      // start viewer
      viewer.init( &argc, argv );
   }

   return 0;
}

