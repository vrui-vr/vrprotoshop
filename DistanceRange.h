/***********************************************************************
DistanceRange - Class to represent ranges of distances and check point
distances.
Copyright (c) 2002-2020 Oliver Kreylos

This file is part of VR ProtoShop.

VR ProtoShop is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

VR ProtoShop is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along
with VR ProtoShop; if not, write to the Free Software Foundation, Inc.,
59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
***********************************************************************/

#ifndef DISTANCERANGE_INCLUDED
#define DISTANCERANGE_INCLUDED

#include "MDGeometry.h"

namespace MD {

class DistanceRange
	{
	/* Elements: */
	Scalar min,max; // Original distance range
	Scalar min2,max2; // Squared distance range
	
	/* Constructors and destructors: */
	public:
	DistanceRange(Scalar sMin,Scalar sMax)
		:min(sMin),max(sMax),min2(min*min),max2(max*max)
		{
		};
	
	/* Methods: */
	Scalar getMin(void) const
		{
		return min;
		};
	Scalar getMax(void) const
		{
		return max;
		};
	bool isInRange(Scalar distance) const // Checks if a distance is inside the range
		{
		return min<=distance&&distance<=max;
		};
	bool isSqrInRange(Scalar distance2) const // Checks if a squared distance is inside the range
		{
		return min2<=distance2&&distance2<=max2;
		};
	bool areInRange(const Point& p1,const Point& p2) const // Checks if two points are inside distance range
		{
		Scalar distance2=Geometry::sqrDist(p1,p2);
		return min2<=distance2&&distance2<=max2;
		};
	};

}

#endif
