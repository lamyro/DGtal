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
 * @file LambdaMedialAxis.h
 * @brief Linear in time distance transformation
 * @author David Coeurjolly (\c david.coeurjolly@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2012/12/08
 *
 * Header file for module 
 *
 * This file is part of the DGtal library.
 *
 * @see testReducedMedialAxis.cpp, testReducedMedialAxisND.cpp, testReverseDT.cpp
 */

#if defined(LambdaMedialAxis_RECURSES)
#error Recursive header files inclusion detected in ReducedMedialAxis.h
#else // defined(ReducedMedialAxis_RECURSES)
/** Prevents recursive inclusion of headers. */
#define LambdaMedialAxis_RECURSES

#if !defined LambdaMedialAxis_h
/** Prevents repeated inclusion of headers. */
#define LambdaMedialAxis_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <vector>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/kernel/NumberTraits.h"
#include "DGtal/geometry/volumes/distance/PowerMap.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/images/ImageContainerBySTLMap.h"
#include "DGtal/images/Image.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/geometry/volumes/distance/VoronoiMap.h"
#include "DGtal/kernel/sets/CDigitalSet.h"
#include "DGtal/geometry/tools/ImplicitBallTools.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class LambdaMedialAxis
  /**
   * Description of template class 'LambdaMedialAxis' <p>
   * \brief Aim: Implementation of the separable medial axis
   * extraction.
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
   * The output is an image associating ball radii (weight of the
   * power map site) to maximal ball centers. Most methods output a
   * lightweight proxy to an image container (of type ImageContainer,
   * see below).
   *
   * @note Following ReverseDistanceTransformation, the input shape is
   * defined as points with negative power distance.
   *
   * @tparam TPowerMap any specialized PowerMap type @tparam
   * TImageContainer any model of CImage to store the medial axis
   * points (default: ImageContainerBySTLVector).
   *
   * @see testReducedMedialAxis.cpp
   */

  

  template <typename TSet,
            typename TDistanceTransformation,
            typename TPowerSeparableMetric,
            typename TImageContainer =  ImageContainerBySTLMap<typename TDistanceTransformation::Domain, double> >
  struct LambdaMedialAxis
  {
    //MA Container
    typedef
      typename DigitalSetSelector < typename TDistanceTransformation::Domain,
	       SMALL_DS + HIGH_ITER_DS >::Type SmallSet;
    typedef Image<TImageContainer> Type;
    typedef TDistanceTransformation DT;
    typedef VoronoiMap<typename TSet::Space, SmallSet, typename DT::SeparableMetric> VMap;
    typedef PowerMap<TImageContainer, TPowerSeparableMetric> PMap;
    //typedef ImplicitBall<typename TSet::Space> ImplicitBall;

    /**
     * Extract reduced medial axis from a power map.
     * This methods is in @f$ O(|powerMap|)@f$.
     *
     * @param aPowerMap the input powerMap
     *
     * @return a lightweight proxy to the ImageContainer specified in
     * template arguments.
     */
    static
    Type getLambdaMedialAxis(TSet & set, typename DT::SeparableMetric & metric, double lambda) {

      SmallSet smallSet(set.domain());
      for (auto it : set) {
        if (set(it)) {
          smallSet.insert(it);
        }
      }
      TImageContainer *computedMA = new TImageContainer( smallSet.domain() );
      TImageContainer *squaredDT = new TImageContainer( smallSet.domain() );
      Z2i::SmallObject8_4 setSmallObject(Z2i::dt8_4, smallSet);
      VMap vmap(smallSet.domain(), smallSet, metric);
      DT dt(smallSet.domain(), smallSet, metric);
      for (typename DT::Domain::ConstIterator it = dt.domain().begin(); it != dt.domain().end(); it++) {
        squaredDT->setValue((*it), pow(dt(*it), 2));
      }
      PMap aPowerMap(squaredDT->domain(), squaredDT, Z2i::l2PowerMetric);
      for (typename PMap::Domain::ConstIterator it = aPowerMap.domain().begin(),
              itend = aPowerMap.domain().end(); it != itend; ++it)
        {
          const auto v  = aPowerMap( *it );
          const auto pv = aPowerMap.projectPoint( v );
          
          if ( aPowerMap.metricPtr()->powerDistance( *it, v, aPowerMap.weightImagePtr()->operator()( pv ) )
                      < NumberTraits<typename PMap::PowerSeparableMetric::Value>::ZERO ) {
            
            Z2i::Point vmapPoint(vmap(*it));
            double vMapDistance = pow(vmapPoint[0] - (*it)[0], 2) + pow(vmapPoint[1] - (*it)[1], 2);
            std::vector<Z2i::Point> setForPotential;
            SmallSet neighborhood(setSmallObject.neighborhood(*it).pointSet());
            for (auto itbis : neighborhood) {
              Z2i::Point vMapDuVoisin(vmap(itbis));
              double distanceAuVmapDuVoisin = pow(vMapDuVoisin[0] - (*it)[0], 2) + pow(vMapDuVoisin[1] - (*it)[1], 2);
              if (distanceAuVmapDuVoisin == vMapDistance) {
                setForPotential.push_back(vMapDuVoisin);
              }
            }
            SmallSet R(set.domain());
            ImplicitBall<typename TSet::Space> enclosing(smallestEuclideanEnclosingBall(setForPotential, R));
            if (enclosing.radius() > lambda) {
              computedMA->setValue( v, aPowerMap.weightImagePtr()->operator()( pv ) );
            }
          
        }
      }
      return Type( computedMA );
    }
  }; // end of class LambdaMedialAxis

} // namespace DGtal

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined LambdaMedialAxis_h

#undef LambdaMedialAxis_RECURSES
#endif // else defined(LambdaMedialAxisdesign pa_RECURSES)