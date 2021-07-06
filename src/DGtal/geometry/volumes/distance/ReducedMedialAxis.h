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
 * @file ReducedMedialAxis.h
 * @brief Linear in time distance transformation
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/12/08
 *
 * Header file for module ReducedMedialAxis.cpp
 *
 * This file is part of the DGtal library.
 *
 * @see testReducedMedialAxis.cpp, testReducedMedialAxisND.cpp, testReverseDT.cpp
 */

#if defined(ReducedMedialAxis_RECURSES)
#error Recursive header files inclusion detected in ReducedMedialAxis.h
#else // defined(ReducedMedialAxis_RECURSES)
/** Prevents recursive inclusion of headers. */
#define ReducedMedialAxis_RECURSES

#if !defined ReducedMedialAxis_h
/** Prevents repeated inclusion of headers. */
#define ReducedMedialAxis_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/PointHashFunctions.h"
#include "DGtal/kernel/NumberTraits.h"
#include "DGtal/geometry/volumes/distance/CPowerSeparableMetric.h"
#include "DGtal/geometry/volumes/distance/PowerMap.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class ReducedMedialAxis
  /**
   * Description of template class 'ReducedMedialAxis' <p>
   * \brief Aim: Implementation of the separable medial axis
   * extraction for the l_2 metric.
   *
   * This utility struct extract medial axis balls from a
   * PowerMap. Basically, each (weighted) site of the PowerMap defines
   * a digital maximal ball if its digital power cell restricted to
   * the input shape is not empty @cite dcoeurjo_pami_RDMA .
   *
   *        Optimal Separable Algorithms to Compute the Reverse
   *        Euclidean Distance Transformation and Discrete Medial Axis in
   *        Arbitrary Dimension, D. Coeurjolly and A. Montanvert, IEEE
   *        Transactions on Pattern Analysis and Machine Intelligence,
   *        29(3):437-448, 2007.
   *
   * The output is a vector of digital ball, with digital center and an
   * integer @e squared @e radius.
   *
   *
   * @note Following ReverseDistanceTransformation, the input shape is
   * defined as points with negative power distance.
   *
   * @tparam TPowerMap any specialized PowerMap type 
   *
   * @see testReducedMedialAxis.cpp
   */
  template <typename TPowerMap>
  struct ReducedMedialAxis
  {
    ///Point type
    typedef typename TPowerMap::Point Point;
    
    ///Weight type
    typedef typename TPowerMap::Weight SquaredRadius;
    
    ///Digital medial ball type (digital center, integer squared radius)
    typedef std::pair<Point,SquaredRadius> MABall;
    
    ///Medial Axis Constainer
    typedef std::vector<MABall> MABalls;
    
    /**
     * Extract reduced medial axis from a power map.
     * This methods is in @f$ O(|powerMap|)@f$.
     *
     * @param aPowerMap the input powerMap
     *
     * @return a vector of digital balls ({center, squared radius}).
     */
    static
    MABalls computeReducedMedialAxisFromPowerMap(const TPowerMap &aPowerMap)
    {
      //Temporary container, a set of MA ball center
      std::unordered_set<Point> centers;
      
      for (typename TPowerMap::Domain::ConstIterator it = aPowerMap.domain().begin(),
             itend = aPowerMap.domain().end(); it != itend; ++it)
        {
          const auto v  = aPowerMap( *it );
          const auto pv = aPowerMap.projectPoint( v );
          
          if ( aPowerMap.metricPtr()->powerDistance( *it, v, aPowerMap.weightImagePtr()->operator()( pv ) )
                      < NumberTraits<typename TPowerMap::PowerSeparableMetric::Value>::ZERO )
            centers.insert( v  );
        }

      //Final copy
      MABalls computedMA;
      for(auto &center: centers)
        computedMA.push_back({ center, aPowerMap.weightImagePtr()->operator()( aPowerMap.projectPoint( center ) ) });
      
      return computedMA;
    }
  }; // end of class ReducedMedialAxis



} // namespace DGtal

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ReducedMedialAxis_h

#undef ReducedMedialAxis_RECURSES
#endif // else defined(ReducedMedialAxisdesign pa_RECURSES)
