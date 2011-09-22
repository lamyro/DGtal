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

#pragma once

/**
 * @file Display3D.h
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2011/08/08
 *
 * Header file for module Display3D.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(Display3D_RECURSES)
#error Recursive header files inclusion detected in Display3D.h
#else // defined(Display3D_RECURSES)
/** Prevents recursive inclusion of headers. */
#define Display3D_RECURSES

#if !defined Display3D_h
/** Prevents repeated inclusion of headers. */
#define Display3D_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

#include "DGtal/base/Common.h"
#include "DGtal/base/CountedPtr.h"
#include "DGtal/io/Color.h"



//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class Display3D
  /**
   * Description of class 'Display3D' <p>
   * \brief Aim:   This semi abstract class  defines the stream mechanism to
   display 3d primitive (like PointVector, DigitalSetBySTLSet, Object
   ...). The class Viewer3D and Board3DTo2D implement two different
   ways to display 3D objects. The first one (Viewer3D), permits an
   interactive visualisation (based on @i OpenGL ) and the second one
   (Board3DTo2D) provides 3D visualisation from 2D vectorial display
   (based on the CAIRO library)
 
   @see Viewer3D, Board2Dto3D
  
  */
  class Display3D
  {
    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor.
     */
    virtual ~Display3D(){};


  protected:  
    Display3D(){};

    // ----------------------- Interface --------------------------------------
  public:
    enum StreamKey {addNewList, updateDisplay, shiftSurfelVisu};
  
 

    /**
     * Used to set the current fill color
     * @param aColor the fill color.
     **/ 
    virtual void setFillColor(DGtal::Color aColor);


    /**
     * Used to set the line fill color
     * @param aColor the line color.
     **/ 
    virtual void setLineColor(DGtal::Color aColor);

    
    /**
     * Used to get the fill color
     * @return the current fill color.
     **/ 
    
    virtual DGtal::Color getFillColor();
   
    /**
     * Used to get the line color
     * @return the current line color.
     **/ 
    
    virtual DGtal::Color getLineColor();

   
  
    /**
     * Add a new 3D Clipping plane represented by ax+by+cz+d = 0 
     * A maximal of five clipping plane can be added.
     *
     * @param a, b, c, d : plane equation.
     **/
  
    virtual void addClippingPlane(double a, double b, double c, double d, bool drawPlane);
  
  
    /**
     * Set camera up-vector.
     * @param x x coordinate of up-vector.
     * @param y y coordinate of up-vector.
     * @param z z coordinate of up-vector.
     */
    virtual void setCameraUpVector(double , double , double ){}; 
  
    /**
     * Set camera position.
     * @param x x position.
     * @param y y position.
     * @param z z position.
     */
    virtual void setCameraPosition(double , double , double ) {  };
  
    /**
     * Set near and far distance.
     * @param near near distance.
     * @param far far distance.
     */
    virtual void setNearFar(double , double ){};
  

    
    /**
     * Set camera direction.
     * @param x x direction.
     * @param y y direction.
     * @param z z direction.
     */
    virtual void setCameraDirection(double , double , double ) { };

  
  
    /**
     * @param objectName the name of the object (generally obtained
     * with a 'object.styleName()').
     *
     * @return the current mode for the given object name or "" if no
     * specific mode has been set.
     */
    virtual std::string getMode( const std::string & objectName ) const;
  
    /**
     * Used to create a new list containing new 3D objects
     * (useful to use transparency between different objects).
     * 
     **/  

    virtual void createNewLineList();
  

    /**
     * Used to create a new list containing new 3D objects
     * (useful to use transparency between different objects).
     * 
     **/  
  
    virtual void createNewPointList();


    /**
     * Used to create a new list containing new 3D objects
     * (useful to use transparency between different objects).
     * 
     **/  

    virtual void createNewVoxelList(bool depthTest=true);

  
    /**
     * Method to add a specific quad (used by @a addClippingPlane). The normal is computed from the vertex order.
     * @param x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4 the four coordinates of the quad.
     * @param aColor the quad color.
     */
    
    virtual void addQuad(double x1, double y1, double z1,  double x2, double y2, double z2,
			 double x3, double y3, double z3,  double x4, double y4, double z4, 
			 DGtal::Color aColor);
  

    /**
     * Method to add a line to the current display.
     * @param x1, y1, z1, x2, y2, z2  the two extremty line points.
     * @param color the line color.
     * @param with the line width
     *
     */
    
    virtual void addLine(double x1, double y1, double z1, double x2, double y2, double z2, 
			 const DGtal::Color &color=DGtal::Color(20,20,20,200), double width=1.5);
  

    /**
     * Method to add specific voxel. It includes several modes to
     * display the voxel with and without the wire visualisation.
     *
     * @param x, y, z  the voxel center.
     * @param color the voxel color.
     * @param width the voxel width.
     * @param widthWire if true add the wire representation.
     */

    virtual void addVoxel(int x, int y, int z, DGtal::Color color= DGtal::Color(220, 220, 220),
			  double width=0.5,bool withWire=false);
    

    /**
     * Method to add a point to the current display.
     * @param x, y, z  the point.
     * @param color the point color.
     * @param with the point width
     *
     */
    
    virtual void addPoint(double x, double y, double z ,const DGtal::Color &color=DGtal::Color(200,20,20),
			  double size=40);
   
  
    
    /**
     * Specific to display a KSSurfel from Kahlimsky space. The display can
     * take into accounts the sign of the cell.
     *
     * @param x,y,z the surfel center.
     * @param xSurfel, ySurfel , zSurfel  specify if the surfel has its main face in the direction of
     *                                     the x-axis, y-axis or z-axis.
     * @param sizeShiftFactor set the distance between the display of the surfel and potential KSVoxel.
     * @param isSigned to specify if we want to display an signed or unsigned Cell.
     * @param aSign if @ref isSigned is true it will be used to apply a different displays 
     *                             according this boolean  parameter 
     *                              (if @param aSign=true oriented in the direct axis orientation).
     * @param basicMode if true, a basic mode to display KSSurfel are used (i.e just a simple surfel face).  
     * 
     */
    
    virtual void addKSSurfel(double x, double y, double z, 
			     bool xSurfel, bool ySurfel, bool zSurfel, double sizeShiftFactor, 
			     bool isSigned= false, bool aSign=true, bool basicMode=false);
    
    /**
     * Add a KSVoxel from the Kahlimsky space.
     * 
     * @param x, y, z the center of the KSVoxel.
     * 
     */
    
    virtual void addKSVoxel(int x, int y, int z);
  
    
    /**
     * Add a KSPoint from the Kahlimsky space.
     * 
     * @param x, y, z the center of the KSVoxel.
     * @param size the point size (default= 0.1)
     * @param isSigned to specify if we want to display an signed or unsigned Cell point.
     * @param aSign if @ref isSigned is true it will be used to apply a different displays 
     *                             according this boolean  parameter 
     *                              (if @param aSign=true display a cross else, display a small cylinder.).     
     */
    
    virtual void addKSPointel(double x, double y, double z, double size=0.1,
			      bool isSigned=false, bool aSign=true);
  

    
    /**
     * Add a KSLinel from the Kahlimsky space. If the KSlinel is
     * signed its display will difffers acoording its sign (display as
     * a cone) else it will be displayed as simple cylinder.
     * 
     * @param x1, y1, z1 first point of the extremity point of the KSLinel.
     * @param x2, y2, z2 second point of the extremity point of the KSLinel.

     * @param width the width of the KSLinel representation (of its associated cylinder (default= 0.02))
     * @param isSigned to specify if we want to display an signed or unsigned Cell Linel.
     * @param aSign if @ref isSigned is true it will add the KSLinel reprensented by a cone oriented in 
     *              the direct axis orientation.   
     * 
     */
    
    virtual void addKSLinel(double x1, double y1, double z1,
			    double x2, double y2, double z2,
			    double width=0.02, bool isSigned=false, bool aSign=true);
  

    /**
     * Used to update the scene bounding box when objects are added. 
     *
     * @param x, y, z the coordinates to be taken into accounts. 
     *
     **/
  
    void updateBoundingBox(double x, double y, double z);
  

  
    /**
     * Draws the drawable [object] in this board. It should satisfy
     * the concept CDrawableWithViewer3D, which requires for instance a
     * method selfDraw( Viewer3D & ).
     *
     * @param object any drawable object.
     * @return a reference on 'this'.
     */
    template <typename TDrawableWithDisplay3D>
    Display3D & operator<<( const  TDrawableWithDisplay3D & object );
  


  
    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;
  
    /**
     * The associated map type for storing possible modes used for
     * displaying for digital objects.
     */
    typedef std::map< std::string, std::string > ModeMapping;
  
  
  
    /**
     * The associated map type for storing the default styles of
     * digital objects.
     */
    typedef std::map< std::string,CountedPtr<DrawableWithDisplay3D> > StyleMapping;
  



    // ------------------------- Protected Datas ------------------------------
  public:
  
    
    ModeMapping myModes;  
    /**
     * For instance, may associate a new style object T1 to the class
     * "HyperRectDomain": myStyles[ "HyperRectDomain" ] = T1.
     *
     * One can also store a new style T2 for a specific mode used for
     * drawing a class:  myStyles[ "HyperRectDomain/Paving" ] = T2.
     *
     * Modes may only be used in objects implementing the concept
     * CDrawableWithBoard2D.
     */
    StyleMapping myStyles;
  
  
  
    double  myBoundingPtUp [3];
    double  myBoundingPtLow [3];

  
  

 
    // ------------------------- Private Datas --------------------------------
  private:


  protected:


    /// Structure used to display KSPoint in 3D
    /// @see addKSPointel 
    ///
    
    struct pointD3D{
      double  x, y,z;
      unsigned int R,G,B,T;
      bool isSigned;
      bool signPos;
      double size;
    };


    /// Structure used to display KSLine in 3D
    /// @see addKSLinel 
    ///
    
  
    struct lineD3D{
      double x1, y1, z1;
      double x2, y2, z2;
      double width;
      unsigned int R,G,B,T;
      bool isSigned;
      bool signPos;
    };
    


    /**
     * Defines the 3D voxel.
     */

    struct voxelD3D{      
      ///  The center coordinate of the voxel.
      ///
      int x, y,z;
      
      ///  The display color of the voxel.
      ///
      unsigned int R,G,B,T;
      
      /// The width of a voxel face 
      ///
      double width;
    };
  
  
    /**
     * Used to define clipping planes (it uses the quadD3D structure)
     * @see Display3D, Viewer3D, Board3DTo2D, quadD3D
     **/ 
  
    struct clippingPlaneD3D{
      double a,b,c,d;
    };

  
    /**
     * This structure is used to display clipping planes and the
     * components of the myKSSurfelList (allowing to set normal and
     * color).
     * @see Display3D, Viewer3D, Board3DTo2D
     **/ 

    struct  quadD3D{
      double x1,y1,z1;
      double x2,y2,z2;
      double x3,y3,z3;
      double x4,y4,z4;    
      double nx, ny, nz;
      unsigned int R,G,B,T;
    };


  protected:


    DGtal::Color myCurrentFillColor;

    DGtal::Color myCurrentLineColor;


    /// Used to specialized visualisation with KS surfels/voxels.
    ///

    double myCurrentfShiftVisuKSSurfels;

  
    /// Used to represent all the list used in the display.
    /// 

    std::vector< std::vector<voxelD3D> > myVoxelSetList;
  
  
    /// Used to represent all the list of line primitive
    ///
  
    std::vector< std::vector<lineD3D> > myLineSetList;
  
  
    /// Used to represent all the list of point primitive
    /// 
  
    std::vector< std::vector<pointD3D> > myPointSetList;

  
    /// Represent all the clipping planes added to the scene (of maxSize=5).
    ///

    std::vector< clippingPlaneD3D > myClippingPlaneList;
 
  
 
    /// For saving all surfels of Khalimsky space (used to display Khalimsky Space Cell)
    ///

    std::vector< quadD3D > myKSSurfelList;
  
  
    /// For saving all pointels of Khalimsky space (used to display Khalimsky Space Cell)
    ///

    std::vector< pointD3D > myKSPointelList;


    /// For saving all linels of Khalimsky space (used to display Khalimsky Space Cell)
    ///

    std::vector< lineD3D > myKSLinelList;
  
  
    /// Represent all the drawed planes
    ///

    std::vector< quadD3D > myQuadList;


    /// Used to define if GL_TEST_DEPTH is used. 
    ///
  
    std::vector<bool> myListVoxelDepthTest;



    // ------------------------- Hidden services ------------------------------

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */


  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    Display3D ( const Display3D & other );
  
    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    Display3D & operator= ( const Display3D & other );
  
    // ------------------------- Internals ------------------------------------
  private:



  
  }; // end of class Display3D

  /**
   * Calculate the cross product of two 3d vectors and return it.
   * @param dst destination vector.
   * @param srcA source vector A.
   * @param srcB source vector B.
   */
  static void cross (double dst[3], double srcA[3], double srcB[3]);

  /**
   * Normalize the input 3d vector.
   * @param vec source & destination vector.
   */
  static void normalize (double vec[3]);


  /**
   * Overloads 'operator<<' for displaying objects of class 'Display3D'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'Display3D' to write.
   * @return the output stream after the writing.
   */
  std::ostream&
  operator<< ( std::ostream & out, const Display3D & object );





} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/io/Display3D.ih"


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined Display3D_h

#undef Display3D_RECURSES
#endif // else defined(Display3D_RECURSES)