/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file testSurfaceHelper.cpp
 * @ingroup Tests
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2014/04/07
 *
 * Functions for testing class SurfaceHelper.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include "DGtal/base/Common.h"
#include "ConfigTest.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/io/boards/Board2D.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class SurfaceHelper.
///////////////////////////////////////////////////////////////////////////////
/**
 * Example of a test. To be completed.
 *
 */


bool testComputeInterior()
{
  unsigned int nbok = 0;
  unsigned int nb = 0;
  typedef  KhalimskySpaceND<2, int>::Cell Cell;
  typedef  KhalimskySpaceND<2, int>::SCell SCell;

  trace.beginBlock ( "Testing computation of interior from a freeman chain ..." );
  std::string freemanChainFilename = testPath + "samples/contourS.fc";
  std::string freemanChainFilename2 = testPath + "samples/SmallBall.fc";
  fstream fst;
  fst.open (freemanChainFilename.c_str(), ios::in);
  FreemanChain<Z2i::Space::Integer> fc(fst);
  fst.close();
  fst.open (freemanChainFilename2.c_str(), ios::in);
  FreemanChain<Z2i::Space::Integer> fc2(fst);
  fst.close();
  std::set<SCell> boundaryCell;
  std::set<SCell> boundaryCell2;
  DGtal::KhalimskySpaceND< 2, int > K; 
  int minx, miny, maxx, maxy; 
  fc.computeBoundingBox(minx, miny, maxx, maxy); 
  K.init(Z2i::Point(minx-5,miny-5), Z2i::Point(maxx+5, maxy+5), false);
  FreemanChain<Z2i::Space::Integer>::getInterPixelLinels(K, fc, boundaryCell ); 
  std::set<Cell> interiorCell;
  Surfaces<DGtal::KhalimskySpaceND< 2, int > >::uComputeInterior(K, boundaryCell, interiorCell, false);
  trace.info() << "Interior size: " << interiorCell.size() << " (awaited:  3014)" <<  std::endl;
  

  DGtal::KhalimskySpaceND< 2, int > K2; 
  int minx2, miny2, maxx2, maxy2; 
  fc2.computeBoundingBox(minx2, miny2, maxx2, maxy2); 
  K2.init(Z2i::Point(minx2-5,miny2-5), Z2i::Point(maxx2+5, maxy2+5), false);
  FreemanChain<Z2i::Space::Integer>::getInterPixelLinels(K2, fc2, boundaryCell2 ); 
  std::set<Cell> interiorCell2;
  Surfaces<DGtal::KhalimskySpaceND< 2, int > >::uComputeInterior(K2, boundaryCell2, interiorCell2, false);
  trace.info() << "Interior size2: " << interiorCell2.size() << " (awaited:  60196)" <<  std::endl;
  
  // Displaying interiorCell
  Board2D aBoard, aBoard2;  
  for( std::set<Cell>::const_iterator it=interiorCell.begin();  it!= interiorCell.end(); it++){
    aBoard << K.uCoords (*it); 
  }
  for( std::set<Cell>::const_iterator it=interiorCell2.begin();  it!= interiorCell2.end(); it++){
    aBoard2 << K.uCoords (*it); 
  }

  aBoard<< CustomStyle( (*(boundaryCell.begin())).className(), 
                        new CustomColors(  Color::Red,  Color::None ) );          
  aBoard2<< CustomStyle( (*(boundaryCell.begin())).className(), 
                        new CustomColors(  Color::Red,  Color::None ) );          
  for( std::set<SCell>::const_iterator it= boundaryCell.begin();  it!= boundaryCell.end(); it++){
    aBoard << *it;
  }
  for( std::set<SCell>::const_iterator it= boundaryCell2.begin();  it!= boundaryCell2.end(); it++){
    aBoard2 << *it;
  }
  aBoard << CustomStyle( fc.className(), 
                        new CustomColors(  Color::Red,  Color::None ) );        
  aBoard2 << CustomStyle( fc.className(), 
                        new CustomColors(  Color::Red,  Color::None ) );        
  aBoard << SetMode( fc.className(), "InterGrid" );
  aBoard2 << SetMode( fc.className(), "InterGrid" );
  aBoard << fc;
  aBoard2 << fc2;    
  
  aBoard.saveEPS("testSurfaceHelperComputeInterior.eps");
  aBoard2.saveEPS("testSurfaceHelperComputeInterior2.eps");
  nbok += (interiorCell.size()== 3014 && interiorCell2.size()== 60196)? 1 : 0; 
  nb++;
  trace.info() << "(" << nbok << "/" << nb << ") "
	       << interiorCell.size() << " (interiorCell.size()) == 3014 and "
               << interiorCell2.size() << " (interiorCell2.size()) == 60196 "
               << std::endl;
  trace.endBlock();
  
  return nbok == nb;
}

/**
* Checks that method Surfaces::findABel can take in argument any pair
* of points (one inside, one outside) to determine a boundary surfel
* by dichotomy. 
*/
template <typename KSpace3D>
bool testFindABel()
{
  typedef KSpace3D                   KSpace;
  typedef typename KSpace::Space     Space;
  typedef typename KSpace::Point     Point;
  typedef typename KSpace::SCell     SCell;
  typedef HyperRectDomain<Space>     Domain;
  typedef DigitalSetBySTLSet<Domain> DigitalSet;
  unsigned int nbok = 0;
  unsigned int nb = 0;
  trace.beginBlock ( "Testing Surfaces::findABel." );
  Point pts[ 8 ];  
  Point p1( -5, -5, -5 );
  Point p2(  5,  5,  5 );
  pts[ 0 ] = p1;
  pts[ 1 ] = p2;
  pts[ 2 ] = Point(  5, -5, -5 );
  pts[ 3 ] = Point( -5,  5, -5 );
  pts[ 4 ] = Point(  5,  5, -5 );
  pts[ 5 ] = Point( -5, -5,  5 );
  pts[ 6 ] = Point(  5, -5,  5 );
  pts[ 7 ] = Point( -5,  5,  5 );
  KSpace K; K.init( p1, p2, true );
  Domain domain( p1, p2 );
  DigitalSet aSet( domain );
  Shapes<Domain>::addNorm2Ball( aSet, Point::zero, 2 );
  for ( unsigned int i = 0; i < 8; ++i )
    {
      SCell bel = Surfaces<KSpace>::findABel( K, aSet, pts[ i ], Point::zero );
      trace.info() << "- Exterior point is " << pts[ i ] << std::endl;
      trace.info() << "  - Found bel = " << bel << std::endl;
      ++nb, nbok += K.sDim( bel ) == 2;
      trace.info() << "(" << nbok << "/" << nb << ") "
                   << " K.sDim( bel ) == " <<  K.sDim( bel ) << " (should be 2)" << std::endl;
      SCell vox_in = K.sDirectIncident( bel, K.sOrthDir( bel ) ); 
      SCell vox_out = K.sIndirectIncident( bel, K.sOrthDir( bel ) ); 
      ++nb, nbok += aSet( K.sCoords( vox_in ) ) ? 1 : 0;
      trace.info() << "(" << nbok << "/" << nb << ") "
                   << " vox_in should be inside the shape." << std::endl;
      ++nb, nbok += aSet( K.sCoords( vox_out ) ) ? 0 : 1;
      trace.info() << "(" << nbok << "/" << nb << ") "
                   << " vox_out should be outside the shape." << std::endl;
    }
  return nbok == nb;
}

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int argc, char** argv )
{
  trace.beginBlock ( "Testing class SurfaceHelper" );
  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;

  bool res = testComputeInterior()
    && testFindABel< KhalimskySpaceND<3,int> >();
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
  return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
